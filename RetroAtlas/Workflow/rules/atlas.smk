# One cohort-wide union and one complete-linkage clustering, independent of CompareMatrix.
ATLAS_DIR = str(WORKFLOW_DIR / "scripts/atlas")
ATLAS_COMPARE_DIR = str(WORKFLOW_DIR / "scripts/compare")
ATLAS_SAMPLE_IDS = list(SAMPLES)


def atlas_paths(suffix):
    return [OUTPATH + f"/{sid}/" + suffix.format(sample=sid) for sid in ATLAS_SAMPLE_IDS]


def atlas_spec(wc):
    p = AADV["Population"]
    radius = max(p["vcf_window"], p["rescue_window"], p["fragment_window"]) + p["primary_window"]
    return dict(assembly=ASSEMBLY, primary_min_support=SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"],
        flank_bp=AADV["Megadepth"]["flank_bp"], search_padding=max(0, math.ceil(radius - ADV["Finalization"]["window"])),
        samples=[dict(id=sid, sex=SAMPLES[sid]["SEX"], procedence=SAMPLES[sid]["PROCEDENCE"],
        sex_karyotype=SAMPLES[sid].get("SEX_KARYOTYPE", "UNKNOWN") or "UNKNOWN",
        origin=SAMPLES[sid]["PROCEDENCE"],
        rds=os.path.abspath(OUTPATH+f"/{sid}/Event_Reconstruction/{sid}_Final.rds"),
        bed=os.path.abspath(OUTPATH+f"/{sid}/Event_Reconstruction/{sid}_regions.bed"),
        coverage=os.path.abspath(OUTPATH+f"/{sid}/Depth/{sid}.regions.bed.gz"),
        panel=os.path.abspath(OUTPATH+f"/{sid}/Reference/{sid}_panel.bed") if ADAPTIVE else None,
        mapping=os.path.abspath(OUTPATH+f"/{sid}/Reference/{sid}_chromosomes.tsv")) for sid in ATLAS_SAMPLE_IDS])


rule prepare_atlas:
    input:
        panels=atlas_paths("Reference/{sample}_panel.bed") if ADAPTIVE else [],
        finals=atlas_paths("Event_Reconstruction/{sample}_Final.rds"),
        beds=atlas_paths("Event_Reconstruction/{sample}_regions.bed"),
        maps=atlas_paths("Reference/{sample}_chromosomes.tsv"),
        coverage=atlas_paths("Depth/{sample}.regions.bed.gz"),
        cli=ATLAS_COMPARE_DIR+"/prepare_compare.R",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R")
    output: directory(ATLAS_OUT+"/Inputs")
    params: spec=atlas_spec
    threads: 1
    resources: mem_mb=min(RAM, 16000)
    conda: "../envs/event_reconstruction.yaml"
    log: ATLAS_OUT+"/logs/prepare.log"
    script: "../scripts/run_compare_prepare.py"


rule atlas_megadepth:
    input:
        prepared=ATLAS_OUT+"/Inputs", bam=bam_for_sample, index=index_for_sample,
        module=ATLAS_COMPARE_DIR+"/evidence.py"
    output: directory(ATLAS_OUT+"/Megadepth/{sample}")
    params: kind="megadepth", settings=AADV["Megadepth"]
    threads: RESOURCES["depth_threads"]
    resources: mem_mb=RESOURCES["depth_mem"]
    conda: "../envs/megadepth.yaml"
    log: ATLAS_OUT+"/logs/megadepth_{sample}.log"
    script: "../scripts/run_compare_evidence.py"


rule atlas_callers:
    input:
        prepared=ATLAS_OUT+"/Inputs", bam=bam_for_sample, index=index_for_sample,
        ref=ref_for_sample, fai=fai_for_sample,
        tandem=OUTPATH+"/{sample}/Reference/{sample}_repeats.bed", module=ATLAS_COMPARE_DIR+"/evidence.py"
    output: ATLAS_OUT+"/VariantCalling/{sample}_{rescue_caller}.vcf"
    wildcard_constraints: rescue_caller="sniffles2|cutesv"
    params: kind=lambda wc: wc.rescue_caller, support=RESCUE_SUPPORT, settings=AADV["Callers"]
    threads: lambda wc: RESOURCES["sniffles_threads" if wc.rescue_caller == "sniffles2" else "cutesv_threads"]
    resources: mem_mb=lambda wc: RESOURCES["sniffles_mem" if wc.rescue_caller == "sniffles2" else "cutesv_mem"]
    conda: "../envs/sv_calling_env_latest.yml"
    log: ATLAS_OUT+"/logs/{rescue_caller}_{sample}.log"
    script: "../scripts/run_compare_evidence.py"


rule merge_atlas:
    input:
        prepared=ATLAS_OUT+"/Inputs",
        cutesv=expand(ATLAS_OUT+"/VariantCalling/{sample}_cutesv.vcf", sample=ATLAS_SAMPLE_IDS),
        sniffles2=expand(ATLAS_OUT+"/VariantCalling/{sample}_sniffles2.vcf", sample=ATLAS_SAMPLE_IDS),
        severus=atlas_paths("VariantCalling/raw/{sample}_severus.vcf"),
        clips=expand(ATLAS_OUT+"/Megadepth/{sample}", sample=ATLAS_SAMPLE_IDS),
        cli=ATLAS_DIR+"/RetroAtlas_population_merge.R",
        modules=[ATLAS_COMPARE_DIR+"/"+n for n in ("RetroAtlas_compare.R", "RetroAtlas_primary_ALT_pairwise.R", "retroatlas_softclip_rescue_v2.R", "RetroAtlas_sexaware_frequency.R")],
        par=str(WORKFLOW_DIR/"data/reference/T2T.PAR.bed"),
        normalizer=ATLAS_COMPARE_DIR+"/normalize_evidence.py",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R"),
        config_helper=str(WORKFLOW_DIR/"scripts/configuration.py"),
        panel_helper=str(WORKFLOW_DIR/"scripts/adaptive_sampling.py")
    output:
        catalog=ATLAS_OUT+"/Atlas_merge.rds",
        internal=ATLAS_OUT+"/population_"+AADV["Population"]["catalog_id"]+"_internal.rds",
        summary=ATLAS_OUT+"/Atlas_merge_summary.tsv",
        states=ATLAS_OUT+"/02_locus_sample_states_PRIVATE.tsv",
        evidence=directory(ATLAS_OUT+"/Evidence"),
        resolved=ATLAS_OUT+"/resolved_atlas.json",
        audits=[ATLAS_OUT+"/"+n+".tsv" for n in ("00_parameters_PRIVATE", "00_stage_checkpoints", "00_evidence_performance_PRIVATE", "01_population_loci", "03_primary_members_PRIVATE", "04_nearby_locus_relationships", "05_primary_edge_audit", "06_vcf_load_audit", "07_evidence_reuse_audit", "08_competing_rescue_audit", "09_public_strata_suppressed")]
    params:
        settings=AADV["Population"], softclip=AADV["Softclip"], assembly=ASSEMBLY, seed=config.get("SEED",1231239),
        support=RESCUE_SUPPORT, primary_support=SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"], sample_ids=ATLAS_SAMPLE_IDS
    threads: CORES
    resources: mem_mb=RAM
    conda: "../envs/event_reconstruction.yaml"
    log: ATLAS_OUT+"/logs/merge.log"
    script: "../scripts/run_atlas.py"


rule validate_atlas:
    input:
        internal=ATLAS_OUT+"/population_"+AADV["Population"]["catalog_id"]+"_internal.rds",
        catalog=ATLAS_OUT+"/Atlas_merge.rds", cli=ATLAS_DIR+"/RetroAtlas_population_validate.R"
    output: ATLAS_OUT+"/validation.txt"
    threads: 1
    resources: mem_mb=RAM
    conda: "../envs/event_reconstruction.yaml"
    log: ATLAS_OUT+"/logs/validate.log"
    shell: "Rscript {input.cli:q} {input.internal:q} {input.catalog:q} > {output:q} 2> {log:q}"


rule export_atlas:
    input:
        validated=ATLAS_OUT+"/validation.txt",
        internal=ATLAS_OUT+"/population_"+AADV["Population"]["catalog_id"]+"_internal.rds",
        catalog=ATLAS_OUT+"/Atlas_merge.rds", evidence=ATLAS_OUT+"/Evidence", prepared=ATLAS_OUT+"/Inputs",
        original_cutesv=atlas_paths("VariantCalling/raw/{sample}_cutesv.vcf"),
        original_sniffles2=atlas_paths("VariantCalling/raw/{sample}_sniffles2.vcf"),
        original_severus=atlas_paths("VariantCalling/raw/{sample}_severus.vcf"),
        normalizer=ATLAS_COMPARE_DIR+"/normalize_evidence.py",
        helpers=str(WORKFLOW_DIR/"scripts/reference_helpers.R"),
        engine=ATLAS_COMPARE_DIR+"/RetroAtlas_compare.R",
        sex_module=ATLAS_COMPARE_DIR+"/RetroAtlas_sexaware_frequency.R",
        par=str(WORKFLOW_DIR/"data/reference/T2T.PAR.bed"),
        module=ATLAS_DIR+"/atlas_outputs.R", cli=ATLAS_DIR+"/export_atlas.R",
        template=ATLAS_DIR+"/report_template.html"
    output:
        catalog=ATLAS_PUBLIC+"/Atlas_merge.rds",
        summary=ATLAS_PUBLIC+"/population_summary.tsv", bed=ATLAS_PUBLIC+"/population_loci.bed",
        rescues=ATLAS_PUBLIC+"/rescues_summary.tsv",
        private_rescues=ATLAS_OUT+"/rescues_PRIVATE.tsv",
        window_records=ATLAS_OUT+"/window_evidence_PRIVATE.tsv",
        window_sources=ATLAS_OUT+"/window_sources_PRIVATE.tsv",
        evidence=directory(ATLAS_OUT+"/WindowEvidence"),
        html=ATLAS_OUT+"/atlas_report.html"
    params:
        sample_ids=ATLAS_SAMPLE_IDS, settings=AADV["WindowBackup"], assembly=ASSEMBLY,
        primary_support=SUPPORT, rescue_support=RESCUE_SUPPORT, min_svlen=MIN_SVLEN,
        min_depth=CALLABILITY["MinDepth"], upper_ratio=CALLABILITY["UpperDepthRatio"], page_size=AADV["Report"]["page_size"], clean_all=CLEAN_ALL
    threads: min(CORES, AADV["WindowBackup"]["workers"])
    resources: mem_mb=RAM
    conda: "../envs/event_reconstruction.yaml"
    log: ATLAS_OUT+"/logs/export.log"
    script: "../scripts/run_atlas_outputs.py"
