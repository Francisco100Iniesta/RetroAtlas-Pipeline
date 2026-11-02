"""BED0 target panels: translate chromosome aliases, never genomic coordinates."""
import bisect
import csv
from collections import defaultdict
from pathlib import Path


def read_bed(path):
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except (IndexError, ValueError) as exc:
                raise ValueError(f'Invalid BED row {number}: {path}') from exc
            if start < 0 or end <= start:
                raise ValueError(f'Invalid BED interval on row {number}: {path}')
            yield chrom, start, end


def merge_intervals(rows):
    merged = []
    for chrom, start, end in sorted(rows):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = chrom, merged[-1][1], max(end, merged[-1][2])
        else:
            merged.append((chrom, start, end))
    return merged


class PanelIndex:
    def __init__(self, path):
        self.intervals = defaultdict(list)
        for chrom, start, end in merge_intervals(read_bed(path)):
            self.intervals[chrom].append((start, end))
        self.starts = {chrom: [a for a, _ in rows] for chrom, rows in self.intervals.items()}

    def contains(self, chrom, pos):
        """VCF positions are 1-based; BED intervals are [start,end)."""
        i = bisect.bisect_left(self.starts.get(chrom, []), pos) - 1
        return i >= 0 and pos <= self.intervals[chrom][i][1]


def normalize_panel(source, mapping_path, catalog_path, output, audit_path):
    with open(mapping_path) as handle:
        mapping = {r['FINAL_CHR']: r for r in csv.DictReader(handle, delimiter='\t')}
    with open(catalog_path) as handle:
        aliases = {a: r['FINAL_CHR'] for r in csv.DictReader(handle, delimiter='\t')
                   for a in r['aliases'].split(';') if a}
    rows = []
    excluded = []
    for chrom, start, end in read_bed(source):
        canonical = aliases.get(chrom)
        if canonical is None or canonical not in mapping:
            raise ValueError(f'Panel chromosome not validated in this assembly: {chrom}')
        m = mapping[canonical]
        if end > int(m['expected_length']):
            raise ValueError(f'Panel interval outside reference: {chrom}:{start}-{end}')
        if m['included'].lower() not in {'true', '1', 't'}:
            excluded.append((chrom, start, end, m['status']))
        else:
            rows.append((m['CHROM'], start, end))
    if not rows:
        raise ValueError('Panel has no intervals on chromosomes included for this sample')
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as handle:
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(merge_intervals(rows))
    with open(audit_path, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['CHROM', 'start', 'end', 'EXCLUSION_REASON'])
        writer.writerows(excluded)


if 'snakemake' in globals():
    normalize_panel(snakemake.input.panel, snakemake.input.mapping, snakemake.input.catalog,
                    snakemake.output.bed, snakemake.output.audit)
