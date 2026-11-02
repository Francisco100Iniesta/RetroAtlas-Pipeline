"""Configuration shared by the DAG and its command wrappers."""
import csv
import re
from pathlib import Path


def resource_plan(cores, ram, allocated_cores=None, allocated_ram=None):
    """Discrete caller thread tiers; RAM is the total job budget, in MB."""
    if any(isinstance(v, bool) or not isinstance(v, int) for v in (cores, ram)):
        raise ValueError('--cores and --resources mem_mb must be integers')
    # Newer Snakemake versions wrap global limits in Resource(value=...).
    allocated_cores = getattr(allocated_cores, 'value', allocated_cores)
    allocated_ram = getattr(allocated_ram, 'value', allocated_ram)
    if allocated_cores is not None:
        cores = min(cores, int(allocated_cores))
    if allocated_ram is not None:
        ram = min(ram, int(allocated_ram))
    if cores < 8 or ram < 32000:
        raise ValueError('RetroAtlas requires --cores >= 8 and --resources mem_mb >= 32000 (MB)')
    tier = 16 if cores >= 16 else 8
    return dict(cores=cores, ram=ram, severus_threads=tier, severus_mem=32000,
                cutesv_threads=tier, cutesv_mem=16000,
                sniffles_threads=8 if cores >= 16 else 4, sniffles_mem=8000,
                samtools_threads=tier, samtools_mem=8000,
                repeatmasker_threads=tier, repeatmasker_mem=32000,
                depth_threads=4, depth_mem=5000, tsd_threads=8, tsd_mem=8000)


def cli_resource_plan(cores, resources):
    """Use the same explicit budgets as Snakemake's global scheduler."""
    cores = getattr(cores, 'value', cores)
    ram = resources.get('mem_mb')
    ram = getattr(ram, 'value', ram)
    if ram is None:
        raise ValueError('Specify total memory with --resources mem_mb=64000 (MB; 64 GB)')
    return resource_plan(cores, ram)


BASE_DEFAULTS = {
    'Outpath': 'Results', 'Single': True, 'Compare': False, 'RetroAtlas': False,
    'VariantCalling': {'Rescue_supp': 1, 'supp': 3, 'min_svlen': 30, 'cluster_size': 500},
    'Callabillity': {'MinDepth': 10, 'UpperDepthRatio': 3},
    'AdaptiveSampling': False, 'AdaptiveSamplingBED': None,
    'T2T': False, 'SEED': 1231239, 'CleanAll': False,
}


def normalize_settings(settings, defaults):
    """Translate public sentence-case keys to the tools' existing internal keys.

    Case-insensitive legacy spellings are accepted, but unknown/duplicate keys
    are rejected instead of silently leaving an option unused.
    """
    if not isinstance(settings, dict):
        raise ValueError('Configuration must be a mapping')
    known = {key.casefold(): key for key in defaults}
    result = {}
    for name, value in settings.items():
        key = known.get(str(name).casefold())
        if key is None:
            raise ValueError(f'Unknown configuration parameter: {name}')
        if key in result:
            raise ValueError(f'Duplicate configuration parameter: {name}')
        result[key] = normalize_settings(value, defaults[key]) if isinstance(defaults[key], dict) else value
    return result


def load_samples(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        first = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter='\t' if '\t' in first else ',')
        reader.fieldnames = [name.strip().upper() for name in (reader.fieldnames or [])]
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate Samples columns after case normalization')
        if not {'ID', 'PATH', 'FILE', 'REF'}.issubset(reader.fieldnames or []):
            raise ValueError('Samples requires ID,PATH,FILE,REF; optional SEX,PROCEDENCE.')
        samples = {}
        for raw in reader:
            if None in raw:
                raise ValueError('Extra fields in Samples row')
            row = {k: (v or '').strip() for k, v in raw.items()}
            sid = row['ID']
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in samples:
                raise ValueError(f'Invalid/duplicate sample ID: {sid!r}')
            row['FILE'] = row['FILE'].upper()
            if row['FILE'] not in {'BAM', 'CRAM'} or Path(row['PATH']).suffix.upper() != '.' + row['FILE']:
                raise ValueError(f'Wrong file extension: sample {sid}; FILE={row["FILE"]!r}, PATH={row["PATH"]!r}. Expected matching BAM or CRAM.')
            if not row['REF']:
                raise ValueError(f'Invalid REF for sample {sid}')
            row['SEX'] = row.get('SEX', 'U').upper()
            row['PROCEDENCE'] = row.get('PROCEDENCE', 'UNKNOW').upper()
            row['SEX'] = {'MALE': 'M', 'FEMALE': 'F', 'UNKNOW': 'U'}.get(row['SEX'], row['SEX'])
            if row['PROCEDENCE'] == 'U':
                row['PROCEDENCE'] = 'UNKNOW'
            if row['SEX'] not in {'M', 'F', 'U'}:
                raise ValueError(f'{sid}: SEX must be M/Male, F/Female or U/Unknow (case insensitive)')
            if row['PROCEDENCE'] not in {'EUROPEAN', 'AFRICAN', 'ASIAN', 'UNKNOW'}:
                raise ValueError(f'{sid}: invalid PROCEDENCE; expected EUROPEAN, AFRICAN, ASIAN or U/Unknow (case insensitive)')
            samples[sid] = row
    if not samples:
        raise ValueError('Samples is empty')
    return samples


def merge_settings(defaults, overrides):
    overrides = normalize_settings(overrides, defaults)
    result = dict(defaults)
    for key, value in overrides.items():
        if key not in defaults:
            raise ValueError(f'Unknown advanced parameter: {key}')
        if isinstance(defaults[key], dict):
            if not isinstance(value, dict):
                raise ValueError(f'{key} must be a mapping')
            result[key] = merge_settings(defaults[key], value)
        else:
            result[key] = value
    return result


def cli_parameters(settings):
    args = []
    for key, value in settings.items():
        if value is None:
            continue
        args += ['--' + key.replace('_', '-'), str(value).lower() if isinstance(value, bool) else str(value)]
    return args


def validate_advanced(settings, defaults, prefix=''):
    """Fail early for malformed numeric settings, before expensive sample jobs."""
    import math
    for key, default in defaults.items():
        val = settings[key]
        name = f'{prefix}.{key}' if prefix else key
        if isinstance(default, dict):
            validate_advanced(val, default, name)
        elif isinstance(default, bool):
            if not isinstance(val, bool):
                raise ValueError(f'{name} must be a YAML boolean')
        elif isinstance(default, (int, float)):
            if isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val) or (val < 0 if default >= 0 else val > 0):
                raise ValueError(f'{name} must be a finite {"nonnegative" if default >= 0 else "nonpositive"} number')
            if isinstance(default, int) and not isinstance(val, int):
                raise ValueError(f'{name} must be an integer')
            if isinstance(default, float) and 0 < default <= 1 and val > 1:
                raise ValueError(f'{name} must lie in [0,1]')
            if key in {'query_threshold', 'model_threshold', 'min_anchor_identity'} and val > 100:
                raise ValueError(f'{name} must lie in [0,100]')
            if key in {'window_bp', 'min_svlen', 'min_tsd', 'max_tsd', 'min_anchor', 'blast_word_size', 'max_repeat_elements'} and val == 0:
                raise ValueError(f'{name} must be positive')
        elif isinstance(default, list):
            if not isinstance(val, list) or not val or not all(isinstance(v, str) and v for v in val):
                raise ValueError(f'{name} must be a nonempty list of strings')
        elif not isinstance(val, str) or not val:
            raise ValueError(f'{name} must be a nonempty string')


def load_comparisons(path, samples):
    """Only explicit pairs; never expand a sample list into all-vs-all."""
    with open(path, encoding='utf-8-sig', newline='') as handle:
        first = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter='\t' if '\t' in first else ',')
        if not {'SampleA', 'SampleB', 'ID_Comparison'}.issubset(reader.fieldnames or []):
            raise ValueError('CompareMatrix requires SampleA,SampleB,ID_Comparison')
        pairs = {}
        seen = set()
        for raw in reader:
            if None in raw:
                raise ValueError('Extra fields in comparison row')
            a, b, cid = [(raw.get(k) or '').strip() for k in ('SampleA', 'SampleB', 'ID_Comparison')]
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', cid) or cid in pairs or cid in samples:
                raise ValueError(f'Invalid/duplicate comparison ID or collision with sample: {cid!r}')
            if a == b or a not in samples or b not in samples:
                raise ValueError(f'{cid}: two distinct samples present in Samples are required: {a}, {b}')
            if frozenset((a, b)) in seen:
                raise ValueError(f'Duplicate comparison pair: {a}, {b}')
            seen.add(frozenset((a, b)))
            pairs[cid] = (a, b)
    if not pairs:
        raise ValueError('Compare requires at least one explicit comparison')
    return pairs


def validate_variant_calling(settings):
    """Do not silently truncate thresholds or interpret YAML booleans as counts."""
    for key in ('Rescue_supp', 'supp', 'min_svlen', 'cluster_size'):
        value = settings[key]
        minimum = 0 if key == 'cluster_size' else 1
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f'VariantCalling.{key} must be an integer >= {minimum}')


def validate_compare(settings):
    """Reject the R command's invalid modes and bounds before sample processing."""
    p = settings['Comparison']
    for key, choices in {'primary_seq_mode': {'off', 'conflicts', 'all'},
                         'clip_input_mode': {'preload', 'region'}}.items():
        if p[key] not in choices:
            raise ValueError(f'Comparison.{key} must be one of {sorted(choices)}')
    if not 0 < p['primary_window'] <= 5000 or p['rescue_window'] < p['primary_window']:
        raise ValueError('Invalid Comparison primary/rescue window')
    if not 0 <= p['duplicate_geometry_window'] <= p['duplicate_window'] <= 100000:
        raise ValueError('Invalid Comparison duplicate windows')
    if not 0 <= p['duplicate_max_alignments_per_source'] <= 1000:
        raise ValueError('Comparison.duplicate_max_alignments_per_source must lie in [0,1000]')
    for key in ('primary_min_size_sim', 'rescue_size_sim', 'duplicate_geometry_size_sim', 'duplicate_align_min_size_sim'):
        if not 0 < p[key] <= 1:
            raise ValueError(f'Comparison.{key} must lie in (0,1]')
    for key in ('primary_seq_identity', 'primary_seq_coverage', 'duplicate_seq_identity', 'duplicate_seq_coverage'):
        if not 0 <= p[key] <= 100:
            raise ValueError(f'Comparison.{key} must lie in [0,100]')
