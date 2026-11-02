# Un calling por caller y muestra; el soporte procede de VariantCalling.supp.
rule sniffles2:
    input:
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        bam=bam_for_sample,
        index=index_for_sample,
        ref=ref_for_sample,
        fai=fai_for_sample,
        tandem=OUTPATH + "/{sample}/Reference/{sample}_repeats.bed"
    output:
        raw=OUTPATH + "/{sample}/VariantCalling/raw/{sample}_sniffles2.vcf"
    params:
        support=SUPPORT,
        min_svlen=MIN_SVLEN
    threads: RESOURCES["sniffles_threads"]
    resources:
        mem_mb=RESOURCES["sniffles_mem"]
    conda: "../envs/sv_calling_env_latest.yml"
    log: OUTPATH + "/{sample}/logs/sniffles2.log"
    shell:
        r"""
        tandem_args=()
        if [[ -s {input.tandem:q} ]]; then
            tandem_args=(--tandem-repeats {input.tandem:q})
        fi
        sniffles --threads {threads} --reference {input.ref:q} --output-rnames \
            --minsupport {params.support} --minsvlen {params.min_svlen} --phase \
            --phase-conflict-threshold 0.1 --input {input.bam:q} --vcf {output.raw:q} \
            "${{tandem_args[@]}}" > {log:q} 2>&1
        """


rule cutesv:
    input:
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        bam=bam_for_sample,
        index=index_for_sample,
        ref=ref_for_sample,
        fai=fai_for_sample
    output:
        raw=OUTPATH + "/{sample}/VariantCalling/raw/{sample}_cutesv.vcf"
    params:
        support=SUPPORT,
        min_svlen=MIN_SVLEN,
        work=OUTPATH + "/{sample}/VariantCalling/tmp"
    threads: RESOURCES["cutesv_threads"]
    resources:
        mem_mb=RESOURCES["cutesv_mem"]
    conda: "../envs/sv_calling_env_latest.yml"
    log: OUTPATH + "/{sample}/logs/cutesv.log"
    shell:
        r"""
        mkdir -p {params.work:q}
        work=$(mktemp -d {params.work:q}/cutesv.XXXXXX)
        trap 'rm -rf "$work"' EXIT
        cuteSV -t {threads} -s {params.support} -l {params.min_svlen} \
            --report_readid --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
            --genotype {input.bam:q} {input.ref:q} {output.raw:q} "$work" > {log:q} 2>&1
        """


rule severus:
    input:
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        bam=bam_for_sample,
        index=index_for_sample,
        ref=ref_for_sample,
        fai=fai_for_sample,
        tandem=OUTPATH + "/{sample}/Reference/{sample}_repeats.bed"
    output:
        raw=OUTPATH + "/{sample}/VariantCalling/raw/{sample}_severus.vcf"
    params:
        support=SUPPORT,
        min_svlen=MIN_SVLEN,
        work=OUTPATH + "/{sample}/VariantCalling/tmp"
    threads: RESOURCES["severus_threads"]
    resources:
        mem_mb=RESOURCES["severus_mem"]
    conda: "../envs/sv_calling_env_latest.yml"
    log: OUTPATH + "/{sample}/logs/severus.log"
    shell:
        r"""
        mkdir -p {params.work:q}
        work=$(mktemp -d {params.work:q}/severus.XXXXXX)
        trap 'rm -rf "$work"' EXIT
        export REF_PATH={input.ref:q}
        tandem_args=()
        if [[ -s {input.tandem:q} ]]; then
            tandem_args=(--vntr-bed {input.tandem:q})
        fi
        severus --target-bam {input.bam:q} --out-dir "$work" --threads {threads} \
            --min-sv-size {params.min_svlen} --min-support {params.support} \
            --output-read-ids "${{tandem_args[@]}}" > {log:q} 2>&1
        test -s "$work/all_SVs/severus_all.vcf"
        cp "$work/all_SVs/severus_all.vcf" {output.raw:q}
        """


rule filter_insertions:
    input:
        raw=OUTPATH + "/{sample}/VariantCalling/raw/{sample}_{caller}.vcf",
        panel=panel_for_sample
    output:
        vcf=OUTPATH + "/{sample}/VariantCalling/{sample}_{caller}_filtered.vcf.gz",
        tbi=OUTPATH + "/{sample}/VariantCalling/{sample}_{caller}_filtered.vcf.gz.tbi"
    threads: BCF_THREADS
    resources:
        mem_mb=min(RAM, 8000)
    conda: "../envs/sv_calling_env_latest.yml"
    log: OUTPATH + "/{sample}/logs/filter_{caller}.log"
    shell:
        r"""
        panel_args=()
        for panel in {input.panel:q}; do
            panel_args=(--targets-file "$panel" --targets-overlap 0)
        done
        bcftools view --threads {threads} -i 'INFO/SVTYPE="INS"' \
            "${{panel_args[@]}}" {input.raw:q} -Oz -o {output.vcf:q} > {log:q} 2>&1
        bcftools index --tbi --force --threads {threads} {output.vcf:q} >> {log:q} 2>&1
        """


rule merge_windows:
    input:
        sniffles=OUTPATH + "/{sample}/VariantCalling/{sample}_sniffles2_filtered.vcf.gz",
        cutesv=OUTPATH + "/{sample}/VariantCalling/{sample}_cutesv_filtered.vcf.gz",
        severus=OUTPATH + "/{sample}/VariantCalling/{sample}_severus_filtered.vcf.gz",
        indexes=lambda wc: expand(OUTPATH + "/{sample}/VariantCalling/{sample}_{caller}_filtered.vcf.gz.tbi",
                                  sample=wc.sample, caller=CALLERS)
    output:
        rds=OUTPATH + "/{sample}/VariantCalling/{sample}_vcf_info.rds",
        bed=OUTPATH + "/{sample}/VariantCalling/{sample}_regions.bed",
        chunks=directory(OUTPATH + "/{sample}/VariantCalling/windows_chunks")
    params:
        chunks=CHUNK_COUNT,
        window=MERGE_DISTANCE
    threads: 1
    resources:
        mem_mb=min(RAM, 8000)
    conda: "../envs/sv_calling_env_latest.yml"
    log: OUTPATH + "/{sample}/logs/merge.log"
    script: "../scripts/merge.R"
