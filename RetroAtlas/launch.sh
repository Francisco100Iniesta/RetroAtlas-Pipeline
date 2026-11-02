# Please activate your Snakemake environment before running this command.
# If your launch script is executed outside this path, please modify the current path (Workflow/Snakefile) to the new one.
snakemake -s Workflow/Snakefile --sdm conda -p "$@" --cores 24 --resources mem_mb=64000
