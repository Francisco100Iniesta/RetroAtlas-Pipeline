"""Parallel BLAST/TSD events; publish a sample TSV only after ALL events succeed."""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

FIELDS_EMPTY = ['event_id', 'tsd_detected', 'tsd_length', 'tsd_sequence']


def run(smk):
    with open(smk.input.bed) as handle:
        events = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
    if any(len(row) != 3 for row in events) or len({r[0] for r in events}) != len(events):
        raise ValueError('TSD BED requires unique IDs and three columns')
    output = Path(smk.output.tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    cfg = dict(smk.params.settings)
    with open(smk.log[0], 'w') as log, tempfile.TemporaryDirectory(prefix='tsd.', dir=output.parent) as root:
        root = Path(root)
        if events:
            # Copy contigs into the temporary workspace; faidx runs exactly once.
            contigs = root / 'contigs.fa'
            contigs.symlink_to(Path(smk.input.contigs).resolve())
            subprocess.run(['samtools', 'faidx', str(contigs)], check=True, stdout=log, stderr=log)
        def process_worker(item):
            idx, (event, region, size) = item
            work = root / str(idx)
            work.mkdir()
            query, subject = work / 'query.fa', work / 'subject.fa'
            for fasta, name, dest in [(smk.input.ref, region, query), (str(contigs), event, subject)]:
                proc = subprocess.run(['samtools', 'faidx', fasta, name], check=True, capture_output=True, text=True)
                lines = proc.stdout.splitlines()
                if len(lines) < 2 or not ''.join(lines[1:]).strip():
                    raise ValueError(f'No sequence for {name}: {proc.stderr}')
                # Ensure that BLAST query/subject share the event ID used by the detector.
                dest.write_text('>' + event + '\n' + '\n'.join(lines[1:]) + '\n')
            for fmt, filename in [('6', 'blast.tsv'), ('0', 'alignment.txt')]:
                subprocess.run(['blastn', '-query', str(query), '-subject', str(subject), '-num_threads', '1',
                                '-task', str(cfg['blast_task']), '-evalue', str(cfg['blast_evalue']),
                                '-word_size', str(cfg['blast_word_size']), '-outfmt', fmt,
                                '-out', str(work / filename)], check=True, stdout=log, stderr=log)
            result = work / 'result.tsv'
            command = [sys.executable, smk.input.detector, '--blast', str(work / 'blast.tsv'),
                       '--alignment', str(work / 'alignment.txt'), '--ins-size', size,
                       '--event-id', event, '--contig-side', 'subject', '--output', str(result)]
            for key, val in cfg.items():
                if key != 'flank_bp' and not key.startswith('blast_'):
                    command += ['--' + key.replace('_', '-'), str(val)]
            subprocess.run(command, check=True, stdout=log, stderr=log)
            with result.open() as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                rows = list(reader)
                if len(rows) != 1 or rows[0].get('event_id') != event:
                    raise ValueError(f'Invalid TSD result for {event}')
                return reader.fieldnames, rows[0]
        def worker(item):
            try:
                return process_worker(item)
            finally:
                shutil.rmtree(root / str(item[0]), ignore_errors=True)
        with ThreadPoolExecutor(max_workers=max(1, smk.threads)) as executor:
            results = list(executor.map(worker, enumerate(events)))
        merged = root / 'merged.tsv'
        with merged.open('w') as handle:
            fields = results[0][0] if results else FIELDS_EMPTY
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for header, row in results:
                if header != fields:
                    raise ValueError('Inconsistent TSD headers')
                writer.writerow(row)
        os.replace(merged, output)


if 'snakemake' in globals():
    run(snakemake)
