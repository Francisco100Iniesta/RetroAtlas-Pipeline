import json
import subprocess
from pathlib import Path

def run(smk):
    spec = dict(smk.params.spec)
    spec['out'] = str(Path(smk.output[0]).resolve())
    spec['helpers'] = str(Path(smk.input.helpers).resolve())
    Path(spec['out']).mkdir(parents=True, exist_ok=True)
    request = Path(spec['out']) / 'request.json'
    request.write_text(json.dumps(spec, indent=2))
    with open(smk.log[0], 'w') as log:
        subprocess.run(['Rscript', smk.input.cli, str(request)], check=True, stdout=log, stderr=subprocess.STDOUT)
if 'snakemake' in globals():
    run(snakemake)
