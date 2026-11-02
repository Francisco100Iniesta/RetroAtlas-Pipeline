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
