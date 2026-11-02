# Shared assembly-specific downloads; never removed by Cleanall.
rule prepare_gene_annotation:
    output: bed=MANE
    params: assembly=ASSEMBLY_KEY, kind="genes"
    threads: 1
    resources: mem_mb=2048
    conda: "../envs/annotations.yaml"
    log: OUTPATH + "/logs/gene_annotation.log"
    script: "../scripts/prepare_annotations.py"


rule prepare_repeat_annotation:
    output: bed=REPEATS
    params: assembly=ASSEMBLY_KEY, kind="repeats"
    threads: 1
    resources: mem_mb=1024
    conda: "../envs/annotations.yaml"
    log: OUTPATH + "/logs/repeat_annotation.log"
    script: "../scripts/prepare_annotations.py"
