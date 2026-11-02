"""Export a private explorer and de-identified catalogue + full-VCF window backup."""
import csv
import json
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
import subprocess
import sys
from pathlib import Path


def normalize_sample(task):
    """Each worker streams one sample's VCFs and owns its audit file."""
    from normalize_evidence import mapping, normalize_vcf
    sid, prepared, norm, originals = task
    m = mapping(prepared / (sid + '.mapping.tsv'))
    audit_path = norm / (sid + '.excluded_chromosomes.tsv')
    with audit_path.open('w') as handle:
        audit = csv.writer(handle)
        for caller, source in originals:
            normalize_vcf(source, norm / (sid + '_' + caller + '.vcf'), m, audit)
    return audit_path


def run(smk):
    sys.path.insert(0, str(Path(smk.input.normalizer).parent))
    private = Path(smk.input.internal).resolve().parent
    public = Path(smk.output.catalog).resolve().parent
    norm = private / 'WindowEvidence'
    norm.mkdir(parents=True, exist_ok=True)
    prepared = Path(smk.input.prepared).resolve()
    with open(prepared / 'samples.tsv') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = list(smk.params.sample_ids)
    if [r['sample_id'] for r in rows] != ids:
        raise ValueError('Window manifest sample order mismatch')
    cores = max(1, int(getattr(smk, 'threads', 1)))
    originals = {caller: getattr(smk.input, 'original_' + caller)
                 for caller in ('cutesv', 'sniffles2', 'severus')}
    if any(len(files) != len(ids) for files in originals.values()):
        raise ValueError('One original VCF per caller and sample is required')
    tasks = [(sid, prepared, norm, [(c, files[i]) for c, files in originals.items()])
             for i, sid in enumerate(ids)]
    workers = min(cores, len(tasks))
    if workers > 1 and 'fork' in multiprocessing.get_all_start_methods():
        with ProcessPoolExecutor(max_workers=workers,
                                 mp_context=multiprocessing.get_context('fork')) as pool:
            audit_paths = list(pool.map(normalize_sample, tasks))
    else:
        audit_paths = [normalize_sample(task) for task in tasks]
    sources = {}
    with open(norm / 'excluded_chromosomes.tsv', 'w') as h:
        audit = csv.writer(h, delimiter='\t')
        audit.writerow(['FILE', 'ORIGINAL_CHROM', 'POS', 'REASON'])
        for path in audit_paths:
            with path.open() as records:
                audit.writerows(csv.reader(records))
        for i, sid in enumerate(ids):
            sources[sid] = []
            for caller in ('cutesv', 'sniffles2', 'severus'):
                dest = norm / (sid + '_' + caller + '.vcf')
                # Original VCFs retain all validated coordinates, even outside the
                # rescue BED or target panel. Coverage will mark those as unassayed.
                sources[sid].append(dict(path=str(dest), caller=caller, source='original', min_support=smk.params.primary_support))
            for caller in ('cutesv', 'sniffles2'):
                sources[sid].append(dict(path=str(Path(smk.input.evidence).resolve() / (sid + '_' + caller + '.vcf')),
                                         caller=caller, source='rescue', min_support=smk.params.rescue_support))
    spec = dict(internal=str(Path(smk.input.internal).resolve()), catalog=str(Path(smk.input.catalog).resolve()),
                helpers=str(Path(smk.input.helpers).resolve()), engine=str(Path(smk.input.engine).resolve()),
                sex_module=str(Path(smk.input.sex_module).resolve()), mapping=str(prepared / 'chromosomes.tsv'),
                module=str(Path(smk.input.module).resolve()), template=str(Path(smk.input.template).resolve()),
                assembly=smk.params.assembly, settings=dict(smk.params.settings), sources=sources, cores=cores,
                min_depth=smk.params.min_depth, upper_ratio=smk.params.upper_ratio, min_svlen=smk.params.min_svlen,
                primary_support=smk.params.primary_support, rescue_support=smk.params.rescue_support,
                page_size=smk.params.page_size, private=str(private), public=str(public),
                public_link=os.path.relpath(public, private))
    request = norm / 'request_PRIVATE.json'
    request.write_text(json.dumps(spec, indent=2))
    Path(smk.log[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(smk.log[0], 'w') as log:
        subprocess.run(['Rscript', str(smk.input.cli), str(request)], check=True, stdout=log, stderr=subprocess.STDOUT)


if 'snakemake' in globals():
    run(snakemake)
