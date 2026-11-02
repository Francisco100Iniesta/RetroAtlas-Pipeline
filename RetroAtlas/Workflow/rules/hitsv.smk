rule hitsv:
    input:
        chunks=OUTPATH + "/{sample}/VariantCalling/windows_chunks",
        validated=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        bam=bam_for_sample,
        index=index_for_sample,
        ref=ref_for_sample,
        fai=fai_for_sample,
        binary=HITSV_BIN
    output:
        results=directory(OUTPATH + "/{sample}/HitSV/results")
    params:
        # Se conserva el -m 50 del comando original de HitSV.
        min_svlen=ADV["HitSV"]["min_svlen"]
    threads: CORES
    resources:
        mem_mb=RAM
    conda: "../envs/preparation.yaml"
    log: OUTPATH + "/{sample}/logs/hitsv.log"
    script: "../scripts/run_hitsv.py"


rule insertion_fasta:
    input:
        results=OUTPATH + "/{sample}/HitSV/results"
    output:
        fasta=OUTPATH + "/{sample}/HitSV/{sample}.fasta",
        rds=OUTPATH + "/{sample}/HitSV/{sample}_pre-data_.rds"
    params:
        min_svlen=MIN_SVLEN
    threads: 1
    resources:
        mem_mb=min(RAM, 16000)
    conda: "../envs/R_rm_env.yaml"
    log: OUTPATH + "/{sample}/logs/fasta.log"
    script: "../scripts/Fasta.R"
