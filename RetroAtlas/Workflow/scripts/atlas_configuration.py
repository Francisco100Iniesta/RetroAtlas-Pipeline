"""Early validation of cohort settings, before scheduling sample processing."""
import re


def validate_atlas(settings, samples, final_window):
    p = settings['Population']
    if p['min_samples'] < 2:
        raise ValueError('Population.min_samples must be at least 2')
    if len(samples) < p['min_samples']:
        raise ValueError(f"RetroAtlas requires at least {p['min_samples']} samples; configure Population.min_samples")
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', p['catalog_id']):
        raise ValueError('Population.catalog_id must be a safe identifier (1–80 characters)')
    for key, choices in dict(public_detail={'summary', 'vectors'},
                             primary_seq_mode={'off', 'conflicts', 'all'},
                             clip_input_mode={'preload', 'region'},
                             clip_evaluation={'unresolved_only', 'all'},
                             clip_error_policy={'abort', 'report'},
                             evidence_scheduler={'dynamic', 'static'}).items():
        if p[key] not in choices:
            raise ValueError(f'Population.{key} must be one of {sorted(choices)}')
    if not p['write_tsv']:
        raise ValueError('RetroAtlas requires Population.write_tsv: true')
    if p['source_clip_cache_enabled'] and p['clip_input_mode'] != 'preload':
        raise ValueError('source_clip_cache_enabled requires clip_input_mode: preload')
    if not 0 < p['primary_window'] <= 500:
        raise ValueError('Population.primary_window must lie in (0, 500]')
    if p['rescue_window'] < p['primary_window'] or not p['primary_window'] <= p['reconcile_window'] <= 100000:
        raise ValueError('Invalid Population rescue/reconcile window')
    if p['duplicate_window'] > 100000 or p['duplicate_geometry_window'] > p['duplicate_window']:
        raise ValueError('Invalid Population duplicate window')
    for key in ('min_public_stratum_n', 'max_candidate_edges', 'progress_batch_loci', 'worker_chunk_loci', 'io_cores'):
        if p[key] < 1:
            raise ValueError(f'Population.{key} must be positive')
    for key, value in p.items():
        if key.endswith('size_sim') or re.fullmatch(r'q\d+_size', key):
            if not 0 < value <= 1:
                raise ValueError(f'Population.{key} must lie in (0, 1]')
        if key.endswith(('seq_identity', 'seq_coverage', 'min_identity', 'min_coverage')) and value > 100:
            raise ValueError(f'Population.{key} must lie in [0, 100]')
    for key in ('min_softclip', 'breakpoint_cluster_bp'):
        if settings['Megadepth'][key] < 1:
            raise ValueError(f'Megadepth.{key} must be positive')
    # A representative can be at either end of the complete-linkage cluster.
    search = max(p['vcf_window'], p['rescue_window'], p['fragment_window']) + p['primary_window']
    if max(final_window, search) + settings['Megadepth']['flank_bp'] < p['clip_window']:
        raise ValueError('Atlas search padding + Megadepth flank must cover Population.clip_window')

    w = settings['WindowBackup']
    if w['workers'] < 1:
        raise ValueError('WindowBackup.workers must be positive')
    if w['window_bp'] < 1 or not 0 < w['min_callable_fraction'] <= 1 or not 0 < w['min_size_similarity'] <= 1:
        raise ValueError('Invalid WindowBackup window/fraction/size similarity')
    if settings['Report']['page_size'] < 1:
        raise ValueError('Report.page_size must be positive')
