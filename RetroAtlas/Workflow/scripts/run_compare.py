import csv
import json
import subprocess
import sys
from pathlib import Path

def run(smk):
    sys.path.insert(0,str(Path(smk.input.normalizer).parent))
    sys.path.insert(0,str(Path(smk.input.config_helper).parent))
    from normalize_evidence import mapping,normalize_vcf,normalize_clips
    from configuration import cli_parameters
    from adaptive_sampling import PanelIndex
    prepared=Path(smk.input.prepared).resolve()
    out=Path(smk.output.rds).resolve().parent
    norm=out/'Evidence';norm.mkdir(parents=True,exist_ok=True)
    with open(prepared/'samples.tsv') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    with open(norm/'excluded_evidence.tsv','w') as handle:
        audit=csv.writer(handle,delimiter='\t');audit.writerow(['FILE','ORIGINAL_CHROM','POS','REASON'])
        for i,row in enumerate(rows):
            sid=row['sample_id'];m=mapping(prepared/(sid+'.mapping.tsv'))
            panel=PanelIndex(row['panel_bed']) if row.get('panel_bed') else None
            for caller in ('cutesv','sniffles2','severus'):
                src=getattr(smk.input,caller)[i]
                dest=norm/(sid+'_'+caller+'.vcf')
                normalize_vcf(src,dest,m,audit,panel);row['vcf_'+caller]=str(dest)
            dest=norm/(sid+'.softclips.detail.tsv.gz')
            normalize_clips(Path(smk.input.clips[i])/(sid+'.softclips.detail.tsv.gz'),dest,m,audit,panel)
            row['softclips']=str(dest)
    manifest=norm/'samples.tsv'
    with open(manifest,'w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    softclip=norm/'softclip_settings.json'
    softclip.write_text(json.dumps(dict(smk.params.softclip),indent=2))
    cfg=dict(smk.params.settings)
    cfg['softclip_config']=str(softclip)
    cfg.update(manifest=str(manifest),outdir=str(out),compare_id=smk.wildcards.comparison,a=rows[0]['sample_id'],b=rows[1]['sample_id'],reference_build=smk.params.assembly,cores=smk.threads,
               chromosome_map=str(prepared/'chromosomes.tsv'),reference_helpers=str(Path(smk.input.helpers).resolve()),
               primary_min_support=smk.params.primary_support,min_svlen=smk.params.min_svlen,require_pass_cover=True,require_pass_region=True,
               callable_min_depth=smk.params.min_depth,callable_max_depth_ratio=smk.params.upper_ratio,vcf_min_support=smk.params.support)
    for name in ('io_cores','clip_cores'):
        cfg[name]=min(cfg[name] or smk.threads,smk.threads)
    cmd=['Rscript',smk.input.cli]+cli_parameters(cfg)
    (out/'resolved_compare.json').write_text(json.dumps(dict(parameters=cfg,softclip=dict(smk.params.softclip),command=cmd),indent=2))
    with open(smk.log[0],'w') as log:
        subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
if 'snakemake' in globals():
    run(snakemake)
