import json
import subprocess
import sys
from pathlib import Path
import yaml



def run(smk):
    sys.path.insert(0, str(Path(smk.input.config_helper).parent))
    from configuration import cli_parameters
    out = Path(smk.output.rds).parent
    out.mkdir(parents=True, exist_ok=True)
    classifier_config = out / f'{smk.wildcards.sample}_classifier_config.yaml'
    classifier_config.write_text(yaml.safe_dump(dict(smk.params.classifier)))
    command = ['Rscript', smk.input.aggregate,
               '--predata', smk.input.predata, '--repeatmasker', smk.input.rm,
               '--sample', smk.wildcards.sample, '--regions', smk.input.regions,
               '--out', str(out), '--classifier-script', smk.input.classifier,
               '--classifier-config', str(classifier_config),
               '--chromosome-map', smk.input.mapping, '--reference-helpers', smk.input.helpers,
               '--assembly', smk.params.assembly, '--sex', smk.params.sex,
               '--procedence', smk.params.procedence,
               '--primary-min-support', str(smk.params.primary_support), '--min-svlen', str(smk.params.min_svlen),
               '--min-depth', str(smk.params.min_depth), '--max-depth-ratio', str(smk.params.upper_ratio),
               '--tsd-flank-bp', str(smk.params.flank), '--simple-repeats-bed', smk.input.repeats,
               '--vcf-cutesv', smk.input.cutesv, '--vcf-severus', smk.input.severus,
               '--vcf-sniffles2', smk.input.sniffles]
    if getattr(smk.input, 'panel', []):
        command += ['--panel-bed', smk.input.panel[0]]
    command += cli_parameters(dict(smk.params.settings))
    (out / f'{smk.wildcards.sample}_parameters.json').write_text(json.dumps(dict(
        reconstruction=dict(smk.params.settings), classifier=dict(smk.params.classifier),
        min_depth=smk.params.min_depth, upper_ratio=smk.params.upper_ratio,
        assembly=smk.params.assembly, sex=smk.params.sex, procedence=smk.params.procedence,
        command=command), indent=2) + '\n')
    with open(smk.log[0], 'w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


if 'snakemake' in globals():
    run(snakemake)
