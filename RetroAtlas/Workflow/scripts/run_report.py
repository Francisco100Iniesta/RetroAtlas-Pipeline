import subprocess

command = ['Rscript', snakemake.input.cli, '--rds', snakemake.input.rds,
           '--regions-bed', snakemake.input.regions, '--out', snakemake.output.html,
           '--vcf', snakemake.output.vcf, '--tsv', snakemake.output.tsv,
           '--chromosome-map', snakemake.input.mapping, '--chromosome-lengths', snakemake.input.lengths,
           '--catalog', snakemake.input.catalog, '--helpers', snakemake.input.helpers,
           '--template', snakemake.input.template, '--reference', snakemake.params.assembly,
           '--primary-min-support', str(snakemake.params.primary_support), '--min-svlen', str(snakemake.params.min_svlen),
           '--min-depth', str(snakemake.params.min_depth), '--upper-depth-ratio', str(snakemake.params.upper_ratio),
           '--promoter-bp', str(snakemake.params.promoter)]
if snakemake.input.mane:
    command += ['--mane', snakemake.input.mane[0], '--mane-assembly', snakemake.params.assembly]
if snakemake.input.panel:
    command += ['--panel-bed', snakemake.input.panel[0]]
with open(snakemake.log[0], 'w') as log:
    subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
