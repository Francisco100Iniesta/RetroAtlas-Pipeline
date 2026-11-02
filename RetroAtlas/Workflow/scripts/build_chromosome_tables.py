"""Build primary chromosome aliases from the checked-in NCBI assembly reports."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / 'data' / 'reference'
SOURCES = {
    'hg38': ('GRCh38.p14', 'GCF_000001405.40', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt'),
    'T2T': ('T2T-CHM13v2.0', 'GCF_009914755.1', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt'),
}


def build():
    for key, (assembly, accession, url) in SOURCES.items():
        with (ROOT / f'{key}.assembly_report.txt').open() as src, (ROOT / f'{key}.chromosomes.tsv').open('w') as out:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(['assembly', 'assembly_accession', 'FINAL_CHR', 'length', 'sequence_name', 'genbank', 'refseq', 'ucsc', 'aliases', 'source'])
            for row in csv.reader((line for line in src if not line.startswith('#')), delimiter='\t'):
                if row[1] != 'assembled-molecule':
                    continue
                name = row[2]
                canonical = 'chrM' if name == 'MT' else f'chr{name}'
                aliases = [row[0], row[4], row[6], row[9], canonical]
                if name == 'MT':
                    aliases += ['M', 'MT', 'chrMT']
                # Unversioned accessions are accepted only with an exact length match.
                aliases += [x.split('.')[0] for x in (row[4], row[6]) if x != 'na']
                aliases = list(dict.fromkeys(x for x in aliases if x != 'na'))
                writer.writerow([assembly, accession, canonical, row[8], row[0], row[4], row[6], row[9], ';'.join(aliases), url])


if __name__ == '__main__':
    build()
