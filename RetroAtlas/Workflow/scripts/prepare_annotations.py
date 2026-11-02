"""Fetch fixed assembly tracks and atomically publish validated BED files.

Gene BEDs retain the original UCSC columns (MANE BED12+13, RefSeq BED12+8).
"""
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

SOURCES = {
    'hg38': {
        'genes': 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb',
        'repeats': 'https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed',
    },
    'T2T': {
        'genes': 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/ncbiRefSeq/ncbiRefSeqSelectCurated.bb',
        'repeats': 'https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed',
    },
}


def validate_bed(path, kind, assembly):
    count = 0
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.rstrip('\n').split('\t')
            minimum = (25 if assembly == 'hg38' else 20) if kind == 'genes' else 3
            if len(fields) < minimum or (kind == 'genes' and len(fields) != minimum):
                raise ValueError(f'{path}:{number}: unexpected annotation schema; expected {minimum} BED fields')
            start, end = map(int, fields[1:3])
            if not 0 <= start < end:
                raise ValueError(f'{path}:{number}: invalid BED interval')
            if kind == 'genes':
                sizes = [int(v) for v in fields[10].rstrip(',').split(',')]
                starts = [int(v) for v in fields[11].rstrip(',').split(',')]
                if (fields[5] not in {'+', '-'} or len(sizes) != int(fields[9])
                        or len(starts) != len(sizes) or not sizes
                        or any(n <= 0 or p < 0 or start+p+n > end for p, n in zip(starts, sizes))):
                    raise ValueError(f'{path}:{number}: invalid BED exon blocks')
            count += 1
    if not count:
        raise ValueError(f'{path}: empty annotation BED')
    return count


def prepare(assembly, kind, destination, bigbed=None):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = SOURCES[assembly][kind]
    with tempfile.TemporaryDirectory(prefix='.annotation-', dir=destination.parent) as tmp:
        raw, bed = Path(tmp)/'source', Path(tmp)/'annotation.bed'
        if bigbed:
            if kind != 'genes':
                raise ValueError('--bigbed is only valid for gene annotations')
            shutil.copyfile(bigbed, raw)
        else:
            # A failed/truncated download never becomes a valid cached output.
            with urllib.request.urlopen(url, timeout=120) as response, raw.open('wb') as handle:
                shutil.copyfileobj(response, handle)
        if kind == 'genes':
            subprocess.run(['bigBedToBed', str(raw), str(bed)], check=True)
        else:
            raw.replace(bed)
        records = validate_bed(bed, kind, assembly)
        digest = hashlib.sha256(bed.read_bytes()).hexdigest()
        metadata = dict(assembly=assembly, kind=kind, source=str(bigbed) if bigbed else url,
                        records=records, sha256=digest)
        meta = Path(tmp)/'source.json'
        meta.write_text(json.dumps(metadata, indent=2)+'\n')
        meta.replace(str(destination)+'.source.json')
        bed.replace(destination)
    return metadata


if 'snakemake' in globals():
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.log[0], 'w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print(prepare(snakemake.params.assembly, snakemake.params.kind, snakemake.output.bed))
elif __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly', choices=SOURCES, required=True)
    parser.add_argument('--kind', choices=['genes', 'repeats'], required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--bigbed', help='Optional already downloaded gene .bb file')
    args = parser.parse_args()
    print(prepare(args.assembly, args.kind, args.output, args.bigbed))
