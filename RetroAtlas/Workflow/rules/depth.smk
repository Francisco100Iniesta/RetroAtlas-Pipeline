rule depth:
    input:
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        bam=bam_for_sample,
        index=index_for_sample,
        ref=ref_for_sample,
        fai=fai_for_sample,
        bed=OUTPATH + "/{sample}/Reference/{sample}_windows.bed"
    output:
        regions=OUTPATH + "/{sample}/Depth/{sample}.regions.bed.gz"
    params:
        prefix=OUTPATH + "/{sample}/Depth/{sample}"
    threads: RESOURCES["depth_threads"]
    resources:
        mem_mb=RESOURCES["depth_mem"]
    conda: "../envs/depth.yaml"
    log: OUTPATH + "/{sample}/logs/depth.log"
    shell:
        "mosdepth -n -t {threads} -f {input.ref:q} -m --by {input.bed:q} "
        "{params.prefix:q} {input.bam:q} > {log:q} 2>&1"
