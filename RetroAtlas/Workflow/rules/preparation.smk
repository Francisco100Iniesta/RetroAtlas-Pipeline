rule validate_reference:
    input:
        alignment=lambda wc: SAMPLES[wc.sample]["PATH"],
        ref=ref_for_sample,
        fai=fai_for_sample,
        catalog=CATALOG,
        repeats=[REPEATS],
        genes=MANE
    output:
        mapping=OUTPATH + "/{sample}/Reference/{sample}_chromosomes.tsv",
        lengths=OUTPATH + "/{sample}/Reference/{sample}_chromosome_lengths.tsv",
        validation=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        windows=OUTPATH + "/{sample}/Reference/{sample}_windows.bed",
        repeats=OUTPATH + "/{sample}/Reference/{sample}_repeats.bed",
        genes=OUTPATH + "/{sample}/Reference/{sample}_genes.bed"
    params:
        sex=lambda wc: SAMPLES[wc.sample]["SEX"],
        policy=ADV["ReferenceValidation"]["mismatch_policy"],
        window=ADV["Depth"]["window_bp"]
    threads: 1
    resources:
        mem_mb=2048
    conda: "../envs/preparation.yaml"
    log: OUTPATH + "/{sample}/logs/reference.log"
    script: "../scripts/reference_validation.py"


if ADAPTIVE:
    rule prepare_panel:
        input:
            panel=PANEL,
            mapping=OUTPATH + "/{sample}/Reference/{sample}_chromosomes.tsv",
            catalog=CATALOG
        output:
            bed=OUTPATH + "/{sample}/Reference/{sample}_panel.bed",
            audit=OUTPATH + "/{sample}/Reference/{sample}_panel_excluded.tsv"
        threads: 1
        resources: mem_mb=1024
        conda: "../envs/preparation.yaml"
        script: "../scripts/adaptive_sampling.py"


rule cram_to_bam:
    input:
        cram=lambda wc: SAMPLES[wc.sample]["PATH"],
        ref=ref_for_sample,
        fai=fai_for_sample,
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json"
    output:
        bam=OUTPATH + "/{sample}/BAM/{sample}.bam",
        bai=OUTPATH + "/{sample}/BAM/{sample}.bam.bai"
    threads: RESOURCES["samtools_threads"]
    resources:
        mem_mb=RESOURCES["samtools_mem"]
    conda: "../envs/preparation.yaml"
    log: OUTPATH + "/{sample}/logs/cram_to_bam.log"
    shell:
        r"""
        if [ {threads} -lt 8 ]; then
            echo 'CRAM conversion requires at least 8 allocated cores; launch with --cores >= 8.' > {log:q}
            exit 1
        fi
        samtools view -@ $(( {threads} - 1 )) -T {input.ref:q} -b -o {output.bam:q} {input.cram:q} > {log:q} 2>&1
        samtools index -@ $(( {threads} - 1 )) {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """
