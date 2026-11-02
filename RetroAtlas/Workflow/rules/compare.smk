# Only explicitly requested A/B pairs. The two directions are evaluated inside R.
wildcard_constraints:
    comparison="|".join(re.escape(c) for c in COMPARISONS)

COMPARE_DIR = str(WORKFLOW_DIR / "scripts/compare")

def pair_paths(wc, suffix):
    return [OUTPATH + f"/{sid}/" + suffix.format(sample=sid) for sid in COMPARISONS[wc.comparison]]

def pair_evidence(wc, kind):
    return [OUTPATH + f"/{wc.comparison}/" + (f"Megadepth/{sid}" if kind == "clips" else f"VariantCalling/{sid}_{kind}.vcf") for sid in COMPARISONS[wc.comparison]]

def compare_spec(wc):
    return dict(assembly=ASSEMBLY, primary_min_support=SUPPORT, min_svlen=MIN_SVLEN, min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"],
        flank_bp=CADV["Megadepth"]["flank_bp"], search_padding=max(0, math.ceil(CADV["Comparison"]["vcf_window"] - ADV["Finalization"]["window"])), samples=[dict(id=sid, sex=SAMPLES[sid]["SEX"], procedence=SAMPLES[sid]["PROCEDENCE"],
        rds=os.path.abspath(OUTPATH+f"/{sid}/Event_Reconstruction/{sid}_Final.rds"),
        bed=os.path.abspath(OUTPATH+f"/{sid}/Event_Reconstruction/{sid}_regions.bed"),
        coverage=os.path.abspath(OUTPATH+f"/{sid}/Depth/{sid}.regions.bed.gz"),
        panel=os.path.abspath(OUTPATH+f"/{sid}/Reference/{sid}_panel.bed") if ADAPTIVE else None,
        mapping=os.path.abspath(OUTPATH+f"/{sid}/Reference/{sid}_chromosomes.tsv")) for sid in COMPARISONS[wc.comparison]])

rule prepare_comparison:
    input:
        panels=lambda wc: pair_paths(wc,"Reference/{sample}_panel.bed") if ADAPTIVE else [],
        finals=lambda wc: pair_paths(wc,"Event_Reconstruction/{sample}_Final.rds"),
        beds=lambda wc: pair_paths(wc,"Event_Reconstruction/{sample}_regions.bed"),
        maps=lambda wc: pair_paths(wc,"Reference/{sample}_chromosomes.tsv"),
        coverage=lambda wc: pair_paths(wc,"Depth/{sample}.regions.bed.gz"),
        cli=COMPARE_DIR+"/prepare_compare.R",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R")
    output: directory(OUTPATH+"/{comparison}/Inputs")
    params: spec=compare_spec
    threads: 1
    resources: mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH+"/{comparison}/logs/prepare.log"
    script: "../scripts/run_compare_prepare.py"

rule comparison_megadepth:
    input:
        prepared=OUTPATH+"/{comparison}/Inputs",
        bam=bam_for_sample, index=index_for_sample,
        module=COMPARE_DIR+"/evidence.py"
    output: directory(OUTPATH+"/{comparison}/Megadepth/{sample}")
    params: kind="megadepth", settings=CADV["Megadepth"]
    threads: RESOURCES["depth_threads"]
    resources: mem_mb=RESOURCES["depth_mem"]
    conda: "../envs/megadepth.yaml"
    log: OUTPATH+"/{comparison}/logs/megadepth_{sample}.log"
    script: "../scripts/run_compare_evidence.py"

rule comparison_callers:
    input:
        prepared=OUTPATH+"/{comparison}/Inputs",
        bam=bam_for_sample, index=index_for_sample, ref=ref_for_sample, fai=fai_for_sample,
        tandem=OUTPATH+"/{sample}/Reference/{sample}_repeats.bed",
        module=COMPARE_DIR+"/evidence.py"
    output: OUTPATH+"/{comparison}/VariantCalling/{sample}_{rescue_caller}.vcf"
    wildcard_constraints: rescue_caller="sniffles2|cutesv"
    params: kind=lambda wc: wc.rescue_caller, support=RESCUE_SUPPORT, settings=CADV["Callers"]
    threads: lambda wc: RESOURCES["sniffles_threads" if wc.rescue_caller == "sniffles2" else "cutesv_threads"]
    resources:
        mem_mb=lambda wc: RESOURCES["sniffles_mem" if wc.rescue_caller == "sniffles2" else "cutesv_mem"]
    conda: "../envs/sv_calling_env_latest.yml"
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
    params: ids=lambda wc: COMPARISONS[wc.comparison], promoter=CADV["Report"]["promoter_bp"], assembly=ASSEMBLY
    threads: 1
    resources: mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: OUTPATH+"/{comparison}/logs/report.log"
    script: "../scripts/run_compare_report.py"
