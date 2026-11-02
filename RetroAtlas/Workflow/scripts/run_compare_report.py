import json
import subprocess
from pathlib import Path

def run(smk):
    a,b=smk.params.ids
    p=Path(smk.input.prepared)
    cmd=['Rscript',smk.input.cli,'--rds',smk.input.rds,'--parameters',smk.input.parameters,
         '--glossary',smk.input.glossary,'--template',smk.input.template,'--out',smk.output[0],
         '--regions-a',str(p/(a+'.regions.bed')),'--regions-b',str(p/(b+'.regions.bed')),
         '--mapping',str(p/'chromosomes.tsv'),'--helpers',smk.input.helpers,'--catalog',smk.input.catalog,
         '--promoter-bp',str(smk.params.promoter),
         '--run-config',json.dumps(smk.params.run_config)]
    if smk.input.mane:
        cmd+=['--mane',smk.input.mane[0],'--mane-assembly',smk.params.assembly]
    with open(smk.log[0],'w') as log:
        subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
if 'snakemake' in globals():
    run(snakemake)
