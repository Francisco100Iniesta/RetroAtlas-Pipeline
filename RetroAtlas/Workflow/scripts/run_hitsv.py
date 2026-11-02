"""Un proceso por chunk; ventanas secuenciales dentro de cada proceso."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def bam_contigs(path):
    """HitSV's -S/-E are zero-based BAM target IDs, not FASTA FAI indices."""
    header = subprocess.run(['samtools', 'view', '-H', str(path)], check=True,
                            capture_output=True, text=True).stdout
    contigs = {}
    for line in header.splitlines():
        if not line.startswith('@SQ\t'):
            continue
        fields = dict(field.split(':', 1) for field in line.split('\t')[1:] if ':' in field)
        name = fields.get('SN')
        length = int(fields.get('LN', 0))
        if not name or name in contigs or length <= 0:
            raise ValueError(f'Invalid/duplicate BAM @SQ: {line}')
        contigs[name] = (len(contigs), length)
    if not contigs:
        raise ValueError(f'No @SQ sequences in BAM: {path}')
    return contigs


def validate_insertions(vcf, chrom, region):
    # BND records may describe another breakend; reconstructed INS must belong
    # to the requested chromosome. Never silently accept a misindexed run.
    with vcf.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                raise RuntimeError(f'Malformed HitSV VCF: {vcf}')
            if 'SVTYPE=INS' in fields[7].split(';') and fields[0] != chrom:
                raise RuntimeError(f'HitSV chromosome mismatch for {region}: '
                                   f'expected {chrom}, got {fields[0]} in {vcf}')


def run(smk):
    results = Path(smk.output.results)
    results.mkdir(parents=True, exist_ok=True)
    chunks = sorted(Path(smk.input.chunks).glob("windows_*.bed"))
    with open(smk.input.fai) as handle:
        contigs = {fields[0]: int(fields[1])
                   for line in handle if (fields := line.split())}
    jobs, seen = [], set()
    for chunk in chunks:
        windows = []
        for line in chunk.read_text().splitlines():
            if not line.strip():
                continue
            chrom, start, end, region = line.split("\t")
            if not re.fullmatch(r"Atlas_[0-9]+", region) or region in seen:
                raise ValueError(f"ID inválido o repetido: {region}")
            seen.add(region)
            length = contigs[chrom]
            start, end = int(start), min(int(end), length)
            if not 0 <= start < end:
                raise ValueError(f"Ventana inválida: {line}")
            windows.append((chrom, start, end, region))
        jobs.append((chunk.stem, windows))

    targets = bam_contigs(smk.input.bam) if seen else {}
    used_chroms = dict.fromkeys(chrom for _, windows in jobs for chrom, _, _, _ in windows)
    for chrom in used_chroms:
        if chrom not in targets or targets[chrom][1] != contigs[chrom]:
            raise ValueError(f'HitSV BAM/FASTA chromosome missing or length mismatch: {chrom}')

    # Logs fuera de los chunks temporales: se conservan para diagnosticar fallos.
    logdir = Path(smk.log[0]).parent / "hitsv_windows"
    logdir.mkdir(parents=True, exist_ok=True)
    with (logdir.parent / 'hitsv_contigs.tsv').open('w') as audit:
        audit.write('CHROM\tBAM_TARGET_ID\tLENGTH\n')
        for chrom in used_chroms:
            idx, length = targets[chrom]
            audit.write(f'{chrom}\t{idx}\t{length}\n')
    scratch = Path(tempfile.mkdtemp(prefix="chunks_", dir=results.parent))

    def process(job):
        chunk_name, windows = job
        folder = scratch / chunk_name
        folder.mkdir()
        for chrom, start, end, region in windows:
            idx = targets[chrom][0]
            vcf = folder / f"{region}.vcf"
            command = [str(smk.input.binary), "call", "-S", str(idx), "-E", str(idx),
                       "-s", str(start), "-F", str(end), "-l", str(smk.input.bam),
                       "-r", str(smk.input.ref), "-p", "ERR_PRONE", "-P", "asm5",
                       "-f", "0", "-H", "-v", "0", "-m", str(smk.params.min_svlen),
                       "-o", str(vcf)]
            with (logdir / f"{region}.log").open("w") as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            if not vcf.is_file():
                raise RuntimeError(f"HitSV no produjo {vcf}")
            validate_insertions(vcf, chrom, region)
        return folder

    with open(smk.log[0], "w") as log:
        log.write(f"Chunks: {len(jobs)}; workers: {smk.threads}\n")
        try:
            with ThreadPoolExecutor(max_workers=max(1, int(smk.threads))) as pool:
                folders = list(pool.map(process, jobs))
            for folder in folders:
                for vcf in folder.glob("*.vcf"):
                    shutil.copy2(vcf, results / vcf.name)
            if len(list(results.glob("*.vcf"))) != len(seen):
                raise RuntimeError("La colección de VCF no coincide con las ventanas procesadas")
            log.write(f"Completado: {len(seen)} ventanas.\n")
        except Exception:
            log.write(f"Fallo; chunks conservados para diagnóstico: {scratch}\n")
            raise
        else:
            shutil.rmtree(scratch)


if "snakemake" in globals():
    run(snakemake)
