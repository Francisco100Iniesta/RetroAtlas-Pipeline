"""Validate reference names AND lengths before running any sample analysis."""
import csv
import contextlib
import traceback
import json
import subprocess
import warnings
from pathlib import Path


def load_catalog(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    aliases = {}
    for row in rows:
        row['length'] = int(row['length'])
        for alias in row['aliases'].split(';'):
            if alias in aliases and aliases[alias]['FINAL_CHR'] != row['FINAL_CHR']:
                raise ValueError(f'Ambiguous chromosome alias: {alias}')
            aliases[alias] = row
    return rows, aliases


def read_fai(path):
    seqs = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2 or fields[0] in seqs or int(fields[1]) <= 0:
                raise ValueError(f'Invalid/duplicate FAI entry: {line!r}')
            seqs[fields[0]] = int(fields[1])
    if not seqs:
        raise ValueError('Empty reference index')
    return seqs


def parse_sam_header(header):
    seqs = {}
    for line in header.splitlines():
        if not line.startswith('@SQ\t'):
            continue
        fields = dict(x.split(':', 1) for x in line.split('\t')[1:] if ':' in x)
        name = fields.get('SN')
        if not name or name in seqs or int(fields.get('LN', 0)) <= 0:
            raise ValueError(f'Invalid/duplicate BAM/CRAM @SQ: {line}')
        seqs[name] = int(fields['LN'])
    if not seqs:
        raise ValueError('No @SQ sequences in BAM/CRAM')
    return seqs


def validate_sequences(fai, bam, catalog, sex='U'):
    rows, aliases = load_catalog(catalog)
    mapped, errors, used = [], [], {}
    for chrom in dict.fromkeys([*fai, *bam]):
        entry = aliases.get(chrom)
        expected = entry['length'] if entry else ''
        status = 'NON_PRIMARY_OR_UNKNOWN'
        if entry:
            canonical = entry['FINAL_CHR']
            if chrom not in fai or chrom not in bam:
                status = 'MISSING_IN_REFERENCE' if chrom not in fai else 'MISSING_IN_ALIGNMENT'
            elif fai[chrom] != expected or bam[chrom] != expected:
                status = 'LENGTH_MISMATCH'
            elif canonical in used:
                status = 'DUPLICATE_PRIMARY_ALIAS'
            else:
                status = 'SEX_F_EXCLUDES_Y' if sex == 'F' and canonical == 'chrY' else 'PASS'
                used[canonical] = chrom
            if status not in {'PASS', 'SEX_F_EXCLUDES_Y'}:
                errors.append(f'{chrom}: {status}; expected={expected}, FASTA={fai.get(chrom)}, BAM/CRAM={bam.get(chrom)}')
        mapped.append(dict(CHROM=chrom, FINAL_CHR=entry['FINAL_CHR'] if entry else '',
                           ASSEMBLY=rows[0]['assembly'], expected_length=expected,
                           fasta_length=fai.get(chrom, ''), alignment_length=bam.get(chrom, ''),
                           status=status, included=status == 'PASS'))
    # Require autosomes + X. Y and mitochondrial sequence are optional.
    required = {f'chr{i}' for i in range(1, 23)} | {'chrX'}
    missing = required - set(used)
    if missing:
        errors.append('Missing validated primary chromosomes: ' + ', '.join(sorted(missing)))
    return mapped, errors


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def normalize_bed(source, destination, catalog, mapping):
    """Rename aliases within the selected assembly, never project coordinates."""
    _, aliases = load_catalog(catalog)
    actual = {x['FINAL_CHR']: x for x in mapping if x['included']}
    kept = skipped = 0
    with open(destination, 'w') as out:
        if source:
            with open(source) as src:
                for line in src:
                    if not line.strip() or line.startswith(('#', 'track', 'browser')):
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 3:
                        raise ValueError(f'Malformed BED row: {line.strip()}')
                    entry = aliases.get(fields[0])
                    target = actual.get(entry['FINAL_CHR']) if entry else None
                    if target is None:
                        skipped += 1
                        continue
                    start, end = int(fields[1]), int(fields[2])
                    if not 0 <= start < end <= target['expected_length']:
                        raise ValueError(f'BED outside selected assembly: {line.strip()}')
                    fields[0] = target['CHROM']
                    out.write('\t'.join(fields) + '\n')
                    kept += 1
    if source and not kept:
        raise ValueError(f'No BED records match validated primary chromosomes: {source}')
    return dict(kept=kept, skipped=skipped)


def run(smk):
    header = subprocess.run(['samtools', 'view', '-H', '-T', smk.input.ref, smk.input.alignment],
                            check=True, capture_output=True, text=True).stdout
    mapping, errors = validate_sequences(read_fai(smk.input.fai), parse_sam_header(header), smk.input.catalog, smk.params.sex)
    write_tsv(smk.output.mapping, mapping, list(mapping[0]))
    included = [r for r in mapping if r['included']]
    order = {r['FINAL_CHR']: i for i, r in enumerate(load_catalog(smk.input.catalog)[0])}
    included.sort(key=lambda r: order[r['FINAL_CHR']])
    write_tsv(smk.output.lengths, [dict(FINAL_CHR=r['FINAL_CHR'], length=r['expected_length']) for r in included], ['FINAL_CHR', 'length'])
    issues = []
    if not smk.input.repeats:
        issues.append('No repeats BED for this assembly; tandem-repeat context unavailable.')
    status = dict(assembly=mapping[0]['ASSEMBLY'], sex=smk.params.sex, errors=errors, warnings=issues,
                  mismatch_policy=smk.params.policy, reference=str(smk.input.ref), alignment=str(smk.input.alignment))
    Path(smk.output.validation).write_text(json.dumps(status, indent=2) + '\n')
    for error in errors + issues:
        warnings.warn('REFERENCE ALARM: ' + error)
    if errors and smk.params.policy == 'error':
        raise ValueError('Reference validation failed. Check hg38/T2T selection and BAM/FASTA assembly. ' + '; '.join(errors))
    if not included:
        raise ValueError('No valid primary chromosomes remain')
    with open(smk.output.windows, 'w') as handle:
        for row in included:
            for start in range(0, row['expected_length'], smk.params.window):
                handle.write(f"{row['CHROM']}\t{start}\t{min(start + smk.params.window, row['expected_length'])}\n")
    status['repeat_annotation'] = normalize_bed(str(smk.input.repeats[0]) if smk.input.repeats else None, smk.output.repeats, smk.input.catalog, mapping)
    if getattr(smk.input, 'genes', None):
        status['gene_annotation'] = normalize_bed(smk.input.genes, smk.output.genes, smk.input.catalog, mapping)
    Path(smk.output.validation).write_text(json.dumps(status, indent=2) + '\n')


if 'snakemake' in globals():
    with open(snakemake.log[0], 'w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            run(snakemake)
        except Exception:
            traceback.print_exc()
            raise

