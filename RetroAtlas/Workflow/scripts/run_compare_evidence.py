import sys
from pathlib import Path

def run(smk):
    sys.path.insert(0,str(Path(smk.input.module).parent))
    from evidence import callers, megadepth
    prepared=Path(smk.input.prepared)
    sid=smk.wildcards.sample
    bed=str(prepared/(sid+'.regions.native.bed'))
    with open(smk.log[0],'w') as log:
        if smk.params.kind=='megadepth':
            megadepth(smk.input.bam,bed,str(prepared/(sid+'.expanded.native.bed')),smk.output[0],sid,smk.threads,dict(smk.params.settings),log)
        else:
            callers(smk.params.kind,smk.input.bam,smk.input.ref,smk.input.fai,bed,smk.input.tandem,smk.output[0],sid,smk.params.support,smk.threads,dict(smk.params.settings),log)
if 'snakemake' in globals():
    run(snakemake)
