    log: OUTPATH+"/{comparison}/logs/{rescue_caller}_{sample}.log"
    script: "../scripts/run_compare_evidence.py"
rule classify_comparison:
    input:
        prepared=OUTPATH+"/{comparison}/Inputs",
        cutesv=lambda wc: pair_evidence(wc,"cutesv"),
        sniffles2=lambda wc: pair_evidence(wc,"sniffles2"),
        severus=lambda wc: pair_paths(wc,"VariantCalling/raw/{sample}_severus.vcf"),
        clips=lambda wc: pair_evidence(wc,"clips"),
        cli=COMPARE_DIR+"/RetroAtlas_compare.R",
        modules=[COMPARE_DIR+"/"+n for n in ("RetroAtlas_primary_ALT_pairwise.R","retroatlas_softclip_rescue_v2.R","RetroAtlas_sexaware_frequency.R")],
        par=str(WORKFLOW_DIR/"data/reference/T2T.PAR.bed"),
        normalizer=COMPARE_DIR+"/normalize_evidence.py",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R"),
        config_helper=str(WORKFLOW_DIR/"scripts/configuration.py"),
        panel_helper=str(WORKFLOW_DIR/"scripts/adaptive_sampling.py")
    output:
        rds=OUTPATH+"/{comparison}/Compare_{comparison}.rds",
        parameters=OUTPATH+"/{comparison}/00_parameters.tsv",
        classified=OUTPATH+"/{comparison}/02_loci_classified.tsv",
        summary=OUTPATH+"/{comparison}/01_category_summary.tsv",
        evidence=directory(OUTPATH+"/{comparison}/Evidence"),
        audits=[OUTPATH+"/{comparison}/"+n+".tsv" for n in ("03_primary_matches","04_rescues","05_review","06_not_callable","07_differentials","08_vcf_load_audit","09_primary_edges","10_softclip_decisions","11_sex_excluded","12_evidence_reuse_qc","13_plausible_internal_locus_links","14_duplicate_guard_summary")]
    params:
        settings=CADV["Comparison"], softclip=CADV["Softclip"], assembly=ASSEMBLY, support=RESCUE_SUPPORT, primary_support=SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"]
    threads: CORES
    resources: mem_mb=RAM
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH+"/{comparison}/logs/compare.log"
    script: "../scripts/run_compare.py"
rule report_comparison:
    input:
        rds=OUTPATH+"/{comparison}/Compare_{comparison}.rds",
        parameters=OUTPATH+"/{comparison}/00_parameters.tsv",
        prepared=OUTPATH+"/{comparison}/Inputs",
        mane=([MANE] if MANE else []),catalog=CATALOG,
        cli=COMPARE_DIR+"/report_compare.R",
        template=COMPARE_DIR+"/report_template.html",
        glossary=COMPARE_DIR+"/COMPARE_D1_GLOSSARY.md",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R")
    output:
        html=OUTPATH+"/{comparison}/{comparison}_report.html",
        annotation_audit=OUTPATH+"/{comparison}/{comparison}_report.html.mane_excluded.tsv"
    params:
        ids=lambda wc: COMPARISONS[wc.comparison], promoter=CADV["Report"]["promoter_bp"], assembly=ASSEMBLY,
        run_config={"configuration": config, "advance_config": ADV, "advance_compare": CADV, "resources": RESOURCES}
    threads: 1
    resources: mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH+"/{comparison}/logs/report.log"
    script: "../scripts/run_compare_report.py"
