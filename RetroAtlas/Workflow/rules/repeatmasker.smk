# Descarga del archivo raiz obligatorio de Dfam 4.0.
# prepare_famdb.py descarga tambien URL + ".md5", valida y descomprime.
rule prepare_famdb_root:
    output:
        database=FAMDB_ROOT_H5
    params:
        url=FAMDB_ROOT_URL
    threads: 1
    resources:
        mem_mb=1024
    conda: "../envs/rm_env.yaml"
    log: OUTPATH + "/logs/famdb_root.log"
    script: "../scripts/prepare_famdb.py"


rule prepare_famdb:
    output:
        database=FAMDB_H5
    params:
        url=FAMDB_URL
    threads: 1
    resources:
        mem_mb=1024
    conda: "../envs/rm_env.yaml"
    log: OUTPATH + "/logs/famdb.log"
    script: "../scripts/prepare_famdb.py"


rule repeatmasker:
    input:
        fasta=OUTPATH + "/{sample}/HitSV/{sample}.fasta",
        database=FAMDB_H5,
        root=FAMDB_ROOT_H5
    output:
        out=OUTPATH + "/{sample}/RepeatMasker/{sample}.fasta.out",
        gff=OUTPATH + "/{sample}/RepeatMasker/{sample}.fasta.out.gff",
        masked=OUTPATH + "/{sample}/RepeatMasker/{sample}.fasta.masked"
    params:
        species="Homo sapiens",
        # RMBlast emplea hasta cuatro threads por búsqueda (-pa).
        parallel=lambda wildcards, threads: max(1, min(24, threads // 4))
    threads: RESOURCES["repeatmasker_threads"]
    resources:
        mem_mb=RESOURCES["repeatmasker_mem"]
    conda: "../envs/rm_env.yaml"
    log: OUTPATH + "/{sample}/logs/repeatmasker.log"
    script: "../scripts/run_repeatmasker.py"
