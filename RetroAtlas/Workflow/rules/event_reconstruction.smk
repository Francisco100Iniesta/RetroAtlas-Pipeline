rule event_reconstruction:
    input:
        panel=panel_for_sample,
        predata=OUTPATH + "/{sample}/HitSV/{sample}_pre-data_.rds",
        rm=OUTPATH + "/{sample}/RepeatMasker/{sample}.fasta.out",
        regions=OUTPATH + "/{sample}/Depth/{sample}.regions.bed.gz",
        mapping=OUTPATH + "/{sample}/Reference/{sample}_chromosomes.tsv",
        repeats=OUTPATH + "/{sample}/Reference/{sample}_repeats.bed",
        cutesv=OUTPATH + "/{sample}/VariantCalling/{sample}_cutesv_filtered.vcf.gz",
        severus=OUTPATH + "/{sample}/VariantCalling/{sample}_severus_filtered.vcf.gz",
        sniffles=OUTPATH + "/{sample}/VariantCalling/{sample}_sniffles2_filtered.vcf.gz",
        aggregate=str(WORKFLOW_DIR / "scripts/event_reconstruction/Retroposon_aggregate_v5_final_fragment_rescue.R"),
        classifier=str(WORKFLOW_DIR / "scripts/event_reconstruction/te_event_classifier_refined_v7.R"),
        helpers=str(WORKFLOW_DIR / "scripts/reference_helpers.R"),
        config_helper=str(WORKFLOW_DIR / "scripts/configuration.py")
    output:
        rds=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_retroatlas_intra_final.rds",
        bed=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_samtools.bed",
        contigs=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_contig.fasta",
        excluded=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_excluded_predata.tsv"
    params:
        settings=ADV["EventReconstruction"], classifier=ADV["TEClassifier"],
        primary_support=SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"],
        flank=ADV["TSD"]["flank_bp"], assembly=ASSEMBLY,
        sex=lambda wc: SAMPLES[wc.sample]["SEX"],
        procedence=lambda wc: SAMPLES[wc.sample]["PROCEDENCE"]
    threads: 1
    resources:
        mem_mb=min(RAM, 32000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH + "/{sample}/logs/event_reconstruction.log"
    script: "../scripts/run_reconstruction.py"


rule tsd:
    input:
        bed=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_samtools.bed",
        contigs=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_contig.fasta",
        ref=ref_for_sample, fai=fai_for_sample,
        detector=str(WORKFLOW_DIR / "scripts/event_reconstruction/tsd_detector_primary_only_float.py")
    output:
        tsv=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_TSD.tsv"
    params:
        settings=ADV["TSD"]
    threads: RESOURCES["tsd_threads"]
    resources:
        mem_mb=RESOURCES["tsd_mem"]
    conda: "../envs/tsd.yaml"
    log: OUTPATH + "/{sample}/logs/tsd.log"
    script: "../scripts/run_tsd.py"


rule finalize_sample:
    input:
        regions=OUTPATH + "/{sample}/Depth/{sample}.regions.bed.gz",
        panel=panel_for_sample,
        rds=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_retroatlas_intra_final.rds",
        tsd=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_TSD.tsv",
        mapping=OUTPATH + "/{sample}/Reference/{sample}_chromosomes.tsv",
        validation=OUTPATH + "/{sample}/Reference/{sample}_validation.json",
        excluded=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_excluded_predata.tsv",
        helpers=str(WORKFLOW_DIR / "scripts/reference_helpers.R"),
        cli=str(WORKFLOW_DIR / "scripts/procesar_intrasample_tsd_BED_todos_loci_CLI.R")
    output:
        rds=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_Final.rds",
        tsv=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_Final.tsv",
        candidates=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_rescue_candidates.tsv",
        nocallable=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_nocallablecandidates.tsv",
        excluded=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_excluded_loci.tsv",
        bed=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_regions.bed",
        events=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_events.tsv",
        membership=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_events_per_region.tsv"
    params:
        out=OUTPATH + "/{sample}/Event_Reconstruction", assembly=ASSEMBLY,
        primary_support=SUPPORT, min_svlen=MIN_SVLEN, min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"],
        panel_arg="--panel-bed" if ADAPTIVE else "",
        sex=lambda wc: SAMPLES[wc.sample]["SEX"], procedence=lambda wc: SAMPLES[wc.sample]["PROCEDENCE"],
        window=ADV["Finalization"]["window"], min_tsd=ADV["Finalization"]["min_tsd"]
    threads: 1
    resources:
        mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH + "/{sample}/logs/finalize.log"
    shell:
        "Rscript {input.cli:q} --rds {input.rds:q} --tsd {input.tsd:q} --mapping {input.mapping:q} "
        "--validation {input.validation:q} --excluded-predata {input.excluded:q} --helpers {input.helpers:q} "
        "--sample {wildcards.sample:q} --sex {params.sex:q} --procedence {params.procedence:q} "
        "--assembly {params.assembly:q} --window {params.window} --min-tsd {params.min_tsd} "
        "--regions {input.regions:q} --min-depth {params.min_depth} --upper-ratio {params.upper_ratio} "
        "--primary-min-support {params.primary_support} --min-svlen {params.min_svlen} "
        "--out {params.out:q} {params.panel_arg} {input.panel:q} > {log:q} 2>&1"


rule report_sample:
    input:
        panel=panel_for_sample,
        rds=OUTPATH + "/{sample}/Event_Reconstruction/{sample}_Final.rds",
        mapping=OUTPATH + "/{sample}/Reference/{sample}_chromosomes.tsv",
        lengths=OUTPATH + "/{sample}/Reference/{sample}_chromosome_lengths.tsv",
        regions=OUTPATH + "/{sample}/Depth/{sample}.regions.bed.gz",
        mane=[OUTPATH + "/{sample}/Reference/{sample}_genes.bed"], catalog=CATALOG,
        helpers=str(WORKFLOW_DIR / "scripts/reference_helpers.R"),
        cli=str(WORKFLOW_DIR / "scripts/report/report_final.R"),
        template=str(WORKFLOW_DIR / "scripts/report/report_template.html")
    output:
        html=OUTPATH + "/{sample}/Report/{sample}_report.html",
        vcf=OUTPATH + "/{sample}/Report/{sample}_insertions.vcf",
        tsv=OUTPATH + "/{sample}/Report/{sample}_report.tsv"
    params:
        assembly=ASSEMBLY, primary_support=SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"],
        promoter=ADV["Report"]["promoter_bp"],
        run_config={"configuration": config, "advance_config": ADV, "resources": RESOURCES}
    threads: 1
    resources:
        mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH + "/{sample}/logs/report.log"
    script: "../scripts/run_report.py"
