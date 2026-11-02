#!/usr/bin/env Rscript
# RetroAtlas Population N v1 — conservative N-sample PRIMARY clustering + evidence states.
# Scientific candidate release: cluster/evidence tests MUST pass on a real R installation.
# Uses the D1 Compare engines (normalization, VCF, softclip, callability, sequence, sex).
# The public summary never contains sample IDs or raw sequences by default.
suppressPackageStartupMessages(library(data.table))
PM_BUILD <- 'RETROATLAS_POPULATION_N_SNAKEMAKE_20261021'
pm_script_dir <- local({
  a <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(a)) dirname(normalizePath(sub('^--file=','',a[[1L]]),mustWork=TRUE)) else getwd()
})
pm_path <- function(x,base=pm_script_dir,must=TRUE) {
  if(length(x)!=1L || is.na(x) || !nzchar(x)) stop('Empty path')
  p <- if(startsWith(x,'/') || startsWith(x,'~')) path.expand(x) else file.path(base,x)
  if(must && !file.exists(p)) stop('Missing file: ',p)
  if(must) normalizePath(p,mustWork=TRUE) else p
}
pm_cfg_defaults <- function() {
  # Merge-specific defaults take priority over unused two-sample identifiers.
  list(manifest='',outdir='RetroAtlas_POPULATION',catalog_id='COHORT',
    reference_build='GRCh38',seed=20261013L,public_detail='summary',overwrite=FALSE,
    min_samples=10L,min_public_stratum_n=5L,cores=16L,io_cores=2L,clip_cores=4L,
    primary_min_support=3,min_svlen=50,
    require_pass_cover=TRUE,require_pass_region=TRUE,
    require_resolved_class=TRUE,require_resolved_event=TRUE,
    require_resolved_completeness=TRUE,
    primary_window=500,primary_min_size_sim=0.80,
    primary_seq_mode='conflicts',primary_seq_identity=95,
    primary_seq_coverage=90,primary_seq_veto=TRUE,
    primary_seq_max_alignments=5000L,primary_seq_max_cells=80000000,
    q10_pos=5,q10_size=0.98,q9_pos=10,q9_size=0.95,
    q8_pos=25,q8_size=0.90,q7_pos=50,q7_size=0.90,
    q6_pos=100,q6_size=0.90,q5_pos=250,q5_size=0.90,
    q4_pos=500,q4_size=0.85,q3_pos=500,q3_size=0.80,
    rescue_window=1000,rescue_size_sim=0.80,
    reconcile_window=1000,reconcile_min_size_sim=0.80,
    reconcile_align=TRUE,reconcile_max_alignments=300L,
    duplicate_guard_enabled=TRUE,duplicate_window=1000,
    duplicate_geometry_window=100,duplicate_geometry_size_sim=0.85,
    duplicate_align_enabled=TRUE,duplicate_align_min_size_sim=0.80,
    duplicate_seq_identity=95,duplicate_seq_coverage=90,
    duplicate_max_alignments_per_source=4L,
    undermerge_seq_identity=95,undermerge_seq_coverage=90,
    vcf_window=500,vcf_size_sim=0.80,vcf_min_support=2,
    vcf_require_pass=FALSE,vcf_allow_single_alt_caller=TRUE,
    clip_enabled=TRUE,clip_input_mode='preload',clip_evaluation='unresolved_only',
    clip_window=5000,clip_core_window=1000,clip_min_identity=90,
    clip_min_coverage=90,clip_min_clips=3,clip_min_cluster=3,
    clip_min_union_bp=200,
    # Progress is logged after each completed group of loci. Errors are fatal by default.
    progress_batch_loci=1000L,clip_error_policy='abort',
    evidence_scheduler='dynamic',worker_chunk_loci=4L,
    source_clip_cache_enabled=TRUE,source_clip_cache_save_rds=TRUE,
    fragment_enabled=TRUE,fragment_window=1000,
    fragment_max_components=4,fragment_max_span=1000,fragment_size_sim=0.80,
    callable_min_depth=0,callable_max_depth_ratio=3,
    chromosome_map='',reference_helpers='../reference_helpers.R',softclip_config='',
    # Reject secondary rescue from noncallable/unknown-quality loci before expensive IO.
    rescue_require_callable=TRUE,
    require_all_vcf_for_differential=TRUE,
    include_sex_chromosomes=TRUE,write_tsv=TRUE,
    max_candidate_edges=3000000L,
    sequence_module='../compare/RetroAtlas_primary_ALT_pairwise.R',
    sex_module='../compare/RetroAtlas_sexaware_frequency.R',
    softclip_module='../compare/retroatlas_softclip_rescue_v2.R',
    compare_engine='../compare/RetroAtlas_compare.R')
}
pm_parse <- function(argv=commandArgs(trailingOnly=TRUE)) {
  d <- pm_cfg_defaults(); x <- d
  if('--help' %in% argv || '-h' %in% argv) {
    cat('Rscript RetroAtlas_population_merge.R --manifest samples.tsv --outdir OUTPUT [--catalog-id COHORT] [--option value ...]\n')
    for(k in names(d)) cat(sprintf(' --%-42s %s\n',gsub('_','-',k),as.character(d[[k]])))
    quit(status=0L)
  }
  if(length(argv)%%2L) stop('Arguments must be --flag value pairs')
  if(length(argv)) for(i in seq(1L,length(argv),by=2L)) {
    flag<-argv[[i]];key<-gsub('-','_',sub('^--','',flag))
    if(!startsWith(flag,'--') || !key %in% names(d)) stop('Unknown option: ',flag)
    old<-d[[key]];z<-argv[[i+1L]]
    if(is.logical(old)) {
      if(!tolower(z) %in% c('true','false','yes','no','1','0','t','f')) stop('Invalid boolean ',flag)
      x[[key]]<-tolower(z) %in% c('true','yes','1','t')
    } else if(is.numeric(old)) {
      v<-suppressWarnings(as.numeric(z));if(!is.finite(v)) stop('Invalid number for ',flag)
      x[[key]]<-if(is.integer(old)) as.integer(v) else v
    } else x[[key]]<-z
  }
  if(!nzchar(x$manifest)) stop('--manifest is required')
  if(!grepl('^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$',x$catalog_id)) stop('Unsafe catalog ID')
  if(!x$public_detail %in% c('summary','vectors')) stop('--public-detail summary|vectors')
  if(!x$clip_input_mode %in% c('preload','region')) stop('--clip-input-mode preload|region')
  if(!x$clip_evaluation %in% c('unresolved_only','all')) stop('--clip-evaluation unresolved_only|all')
  if(!x$clip_error_policy %in% c('abort','report')) stop('--clip-error-policy abort|report')
  if(!x$evidence_scheduler %in% c('dynamic','static')) stop('--evidence-scheduler dynamic|static')
  if(x$worker_chunk_loci<1L) stop('--worker-chunk-loci must be positive')
  if(isTRUE(x$source_clip_cache_enabled) && x$clip_input_mode!='preload')
    stop('--source-clip-cache-enabled true requires --clip-input-mode preload')
  if(!is.finite(x$progress_batch_loci) || x$progress_batch_loci<1L ||
     x$progress_batch_loci!=as.integer(x$progress_batch_loci))
    stop('--progress-batch-loci must be a positive integer')
  if(!x$primary_seq_mode %in% c('off','conflicts','all')) stop('Invalid primary-seq-mode')
  if(x$primary_window<1 || x$primary_window>500 || x$rescue_window<x$primary_window ||
     x$reconcile_window<x$primary_window || x$reconcile_window>100000 ||
     x$duplicate_window<0 || x$duplicate_window>100000 ||
     x$primary_min_size_sim<=0 || x$primary_min_size_sim>1 ||
     x$reconcile_min_size_sim<=0 || x$reconcile_min_size_sim>1 ||
     x$cores<1 || x$io_cores<1 || x$clip_cores<1 ||
     x$min_public_stratum_n<1 || x$max_candidate_edges<1 || x$reconcile_max_alignments<0)
    stop('Invalid parameter range')
  x
}
pm_modules <- function(cfg) {
  e <- new.env(parent=globalenv())
  # Source an already-audited Compare D1 engine as a *library*.
  # Its Rscript entrypoint is guarded by if(sys.nframe()==0).
  sys.source(pm_path(cfg$compare_engine),envir=e)
  e$compare_script_dir <- dirname(pm_path(cfg$compare_engine))
  sys.source(pm_path(cfg$reference_helpers),envir=e)
  sys.source(pm_path(cfg$sequence_module),envir=e)
  sys.source(pm_path(cfg$sex_module),envir=e)
  if(cfg$clip_enabled) {
    cenv <- new.env(parent=e)
    sys.source(pm_path(cfg$softclip_module),envir=cenv)
    for(fn in c('retroatlas_clip_classify','load_softclip_table')) {
      if(!exists(fn,envir=cenv,inherits=FALSE,mode='function')) stop('Missing clip API ',fn)
      assign(fn,get(fn,envir=cenv,inherits=FALSE),envir=e)
    }
    # Exercise the NONEMPTY table-selection API (old preflight only tested empty clips).
    # The clip module's ra_select_table() uses base data.frame indexing.
    probe_dt<-data.table(CLIP_ID='POP_TEST_CLIP',CHR='chr1',POS0=999,POS1=1000,
      SOFTCLIP_LEN=50,SEQUENCE=paste(rep('A',50),collapse=''))
    probe<-cenv$ra_select_table(as.data.frame(probe_dt),'chr1',1000,10)
    if(!is.data.frame(probe) || nrow(probe)!=1L ||
       !identical(as.character(probe$CLIP_ID[[1L]]),'POP_TEST_CLIP'))
      stop('[POP][SOFTCLIP_PREFLIGHT_FAILED] Nonempty table selection failed')
    message('[POP][SOFTCLIP_PREFLIGHT] Nonempty softclip table selection OK')
  }
  for(fn in c('ra_pair_q','ra_normalize_intra','ra_read_vcf','ra_quality','ra_vcf_evaluate',
               'ra_fragment_evaluate','ra_softclip_evaluate','sa_manifest','sa_ploidy',
               'sa_region','ra_gt_info','ra_size_sim','ra_duplicate_guard')) {
    if(!exists(fn,envir=e,inherits=FALSE,mode='function')) stop('Missing Compare D1 API: ',fn)
  }
  e$ra_compare_preflight()
  e$ra_vcf_regression_preflight()
  e$ra_duplicate_regression_preflight()
  if(cfg$clip_enabled) e$ra_softclip_regression_preflight()
  message('[POP][MODULES] Compare D1 functions + ALT, sex and softclip validated at startup')
  e
}
pm_manifest <- function(path,e,cfg) {
  d<-fread(path,na.strings=c('','NA','N/A'),colClasses=c(sample_id='character'))
  req<-c('sample_id','final_rds','vcf_cutesv','vcf_severus','vcf_sniffles2',
    'quality_regions','softclips','assembly')
  miss<-setdiff(req,names(d));if(length(miss)) stop('Missing manifest fields: ',paste(miss,collapse=','))
  if(cfg$min_samples<2L || nrow(d)<cfg$min_samples || anyDuplicated(d$sample_id) || anyNA(d$sample_id))
    stop('Need at least ',cfg$min_samples,' unique samples (Population.min_samples)')
  if(!'sex' %in% names(d)) d[,sex:='U']
  if(!'sex_karyotype' %in% names(d)) d[,sex_karyotype:='UNKNOWN']
  d<-e$sa_manifest(d)
  # sa_manifest normalizes the sole origin field from PROCEDENCE.
  if(anyNA(d$assembly) || any(d$assembly!=cfg$reference_build)) stop('Mixed/reference assembly mismatch')
  basedir<-dirname(path)
  for(k in setdiff(req,c('sample_id','assembly'))) {
    d[,(k):=vapply(as.character(get(k)),function(z) {
      if(is.na(z)||!nzchar(z)) return(NA_character_)
      pm_path(z,base=basedir,must=FALSE)
    },character(1))]
  }
  for(k in c('final_rds','quality_regions'))
    if(any(!file.exists(d[[k]]))) stop('Required file missing [',k,']: ',paste(d[[k]][!file.exists(d[[k]])],collapse=';'))
  if(cfg$clip_enabled && any(!file.exists(d$softclips))) stop('Softclips files missing')
  d
}
pm_checkpoint <- function(outdir,stage,detail='',status='DONE',start=NULL) {
  elapsed <- if(is.null(start)) NA_real_ else as.numeric(difftime(Sys.time(),start,units='secs'))
  x<-data.table(stage=stage,status=status,elapsed_s=round(elapsed,2),detail=as.character(detail))
  path<-file.path(outdir,'00_stage_checkpoints.tsv')
  fwrite(x,path,sep='\t',append=file.exists(path),col.names=!file.exists(path))
  message('[POP][',stage,'][',status,'] ',detail,if(is.finite(elapsed)) sprintf(' | %.1fs',elapsed) else '')
}
pm_first <- function(d,aliases,default=NA) {
  a<-aliases[aliases %in% names(d)]
  if(length(a)) return(d[[a[[1L]]]])
  rep(default,nrow(d))
}
pm_load_sample <- function(row,e,cfg) {
  sid<-row$sample_id[[1L]]
  obj<-e$ra_normalize_intra(row$final_rds[[1L]],sid,cfg)
  d<-obj$loci
  src<-readRDS(row$final_rds[[1L]])
  raw<-as.data.table(if(is.list(src) && 'loci' %in% names(src)) src$loci else src)
  if(is.list(src) && is.data.frame(src$rescue_candidates) && nrow(src$rescue_candidates))
    raw<-rbindlist(list(raw,as.data.table(src$rescue_candidates)),fill=TRUE,use.names=TRUE)
  if(nrow(raw)!=nrow(d)) stop('Lost row order between normalized and raw RDS for ',sid)
  d[,`:=`(polyA=as.character(pm_first(raw,c('polyA','POLYA'),NA_character_)),
    polyA_length=suppressWarnings(as.numeric(pm_first(raw,c('polyA_length','POLYA_LENGTH'),NA_real_))),
    tsd_detected=as.character(pm_first(raw,c('tsd_detected','TSD_DETECTED'),NA_character_)),
    tsd_length=suppressWarnings(as.numeric(pm_first(raw,c('tsd_length','TSD_LENGTH'),NA_real_))),
    tsd_sequence=as.character(pm_first(raw,c('tsd_sequence','TSD_SEQUENCE'),NA_character_)),
    element_type=as.character(pm_first(raw,c('element_type','ELEMENT_TYPE',
      'TE_SUBFAMILY','SVA_SUBFAMILY','ALU_SUBFAMILY','element_subfamily','subfamily'),NA_character_)))]
  quality<-e$ra_quality(row$quality_regions[[1L]],cfg$chromosome_map,d,cfg,row$sex[[1L]])
  d[,source_callable:=quality$callable[match(uid,quality$uid)]]
  d[,primary_eligible:=primary_eligible & source_callable %in% TRUE]
  d[,sex_region:=e$sa_region(CHROM,POS,cfg$reference_build)]
  if(!cfg$include_sex_chromosomes) d[sex_region!='AUTOSOME',primary_eligible:=FALSE]
  d[,ploidy:=e$sa_ploidy(sex_region,row$sex_karyotype[[1L]])]
  d[is.na(ploidy)|ploidy==0L,primary_eligible:=FALSE]
  list(loci=d,fragments=obj$fragments)
}
# All candidate edges are auditable. Range restriction is physical breakpoint +/-500 bp.
pm_primary_edges <- function(primary,e,cfg) {
  p<-copy(primary);p[,chr_rank:=suppressWarnings(as.integer(sub('^chr','',CHROM)))]
  p[CHROM=='chrX',chr_rank:=23L];p[CHROM=='chrY',chr_rank:=24L]
  setorder(p,chr_rank,CHROM,POS,uid)
  entries<-list();k<-0L;tested<-0L
  for(chr in unique(p$CHROM)) {
    x<-p[CHROM==chr];n<-nrow(x);if(n<2L) next
    for(i in seq_len(n-1L)) {
      j<-i+1L
      while(j<=n && x$POS[[j]]-x$POS[[i]]<=cfg$primary_window) {
        if(x$sample_id[[i]]!=x$sample_id[[j]]) {
          tested<-tested+1L
          if(tested>cfg$max_candidate_edges)
            stop('Exceeded --max-candidate-edges; increase after reviewing cohort density')
          q<-e$ra_pair_q(as.list(x[i]),as.list(x[j]),cfg)
          if(!is.null(q)) {
            k<-k+1L
            entries[[k]]<-data.table(uid_a=x$uid[[i]],uid_b=x$uid[[j]],
              sample_a=x$sample_id[[i]],sample_b=x$sample_id[[j]],
              Q=as.integer(q$Q),score=q$score,dpos=q$dpos,size_sim=q$size_sim,
              seq_status=q$seq_status,seq_identity=q$seq_identity,
              seq_coverage=q$seq_coverage,reason=q$reason)
          }
        }
        j<-j+1L
      }
    }
  }
  ans<-if(length(entries)) rbindlist(entries) else data.table(uid_a=character(),uid_b=character(),
    sample_a=character(),sample_b=character(),Q=integer(),score=double(),dpos=double(),
    size_sim=double(),seq_status=character(),seq_identity=double(),seq_coverage=double(),reason=character())
  if(nrow(ans)) setorder(ans,-Q,-score,uid_a,uid_b)
  message('[POP][CANDIDATES] tested=',tested,' | accepted=',nrow(ans))
  ans
}
# Greedy, strongest-edge-first, strict COMPLETE-LINKAGE: every cross-cluster
# member pair must have an acceptable PRIMARY edge. Never chain A~B~C if A!~C.
# Exactly ONE reconstructed PRIMARY per sample per cluster. No haplotype collapse.
pm_cluster <- function(primary,edges,cfg) {
  p<-copy(primary)
  if(!nrow(p)) stop('No eligible PRIMARY observations')
  if(anyDuplicated(p$uid)) stop('Duplicate PRIMARY UID')
  u<-p$uid;sid<-p$sample_id;n<-nrow(p)
  parent<-seq_len(n);members<-lapply(seq_len(n),function(i)i)
  idx<-setNames(seq_len(n),u)
  pairkey<-function(a,b) paste(sort(c(a,b)),collapse='||')
  ok<-new.env(parent=emptyenv(),hash=TRUE)
  if(nrow(edges)) for(i in seq_len(nrow(edges))) assign(
    pairkey(edges$uid_a[[i]],edges$uid_b[[i]]),TRUE,envir=ok)
  root<-function(j) {while(parent[[j]]!=j) j<-parent[[j]];j}
  decisions<-vector('list',nrow(edges))
  for(k in seq_len(nrow(edges))) {
    a<-idx[[edges$uid_a[[k]]]];b<-idx[[edges$uid_b[[k]]]]
    x<-root(a);y<-root(b);accepted<-FALSE;why<-'ALREADY_TOGETHER'
    if(x!=y) {
      aa<-members[[x]];bb<-members[[y]]
      if(length(intersect(sid[aa],sid[bb]))) why<-'SAME_SAMPLE_COMPETITION'
      else if(max(p$POS[c(aa,bb)])-min(p$POS[c(aa,bb)])>cfg$primary_window)
        why<-'COMPLETE_LINKAGE_POSITION_SPAN'
      else {
        allcross<-all(vapply(aa,function(i) all(vapply(bb,function(j)
          exists(pairkey(u[[i]],u[[j]]),envir=ok,inherits=FALSE),logical(1))),logical(1)))
        if(!allcross) why<-'COMPLETE_LINKAGE_PAIR_INCOMPATIBLE'
        else {
          parent[[y]]<-x;members[[x]]<-c(aa,bb);members[[y]]<-integer()
          accepted<-TRUE;why<-'JOINED_COMPLETE_LINKAGE'
        }
      }
    }
    decisions[[k]]<-data.table(uid_a=edges$uid_a[[k]],uid_b=edges$uid_b[[k]],
      accepted=accepted,decision=why)
  }
  all_roots<-vapply(seq_len(n),root,integer(1))
  p[,root_id:=all_roots]
  loc<-p[,.(CHROM=CHROM[[1L]],POS=stats::median(POS),
      POS_MIN=min(POS),POS_MAX=max(POS),SVLEN=stats::median(SVLEN),
      N_PRIMARY=.N,N_PRIMARY_SAMPLES=uniqueN(sample_id),
      CLASS_MODE=names(sort(table(event_class),decreasing=TRUE))[[1L]],
      EVENT_MODE=names(sort(table(event),decreasing=TRUE))[[1L]],
      Q_MIN=NA_integer_),by=root_id]
  loc[,chr_sort:=suppressWarnings(as.integer(sub('^chr','',CHROM)))]
  loc[CHROM=='chrX',chr_sort:=23L];loc[CHROM=='chrY',chr_sort:=24L]
  loc[is.na(chr_sort),chr_sort:=999L]
  setorder(loc,chr_sort,CHROM,POS,root_id)
  loc[,locus_id:=sprintf('LocusAtlas_%07d',seq_len(.N))]
  if(nrow(edges)) {
    for(i in seq_len(nrow(loc))) {
      uu<-p[root_id==loc$root_id[[i]],uid]
      w<-edges[uid_a %in% uu & uid_b %in% uu]
      if(nrow(w)) loc$Q_MIN[[i]]<-min(w$Q)
    }
  }
  p[loc,on='root_id',locus_id:=i.locus_id]
  setcolorder(loc,c('locus_id',setdiff(names(loc),'locus_id')))
  loc[,chr_sort:=NULL]
  p[,root_id:=NULL];loc[,root_id:=NULL]
  if(anyDuplicated(p[,.(locus_id,sample_id)])) stop('COMPLETE LINKAGE failed: duplicate sample in locus')
  audit<-if(nrow(edges)) cbind(copy(edges),rbindlist(decisions)[,.(accepted,decision)]) else data.table()
  list(loci=loc,primary=p,edge_decisions=audit)
}
# Build an inter-locus graph AFTER freezing PRIMARY identity.
# An edge does not merge loci. It records candidate identity/allelic linkage.
pm_links <- function(loci,primary,e,cfg) {
  empty<-data.table(locus_a=character(),locus_b=character(),CHROM=character(),
    distance_bp=double(),size_sim=double(),same_class=logical(),same_event=logical(),
    shared_intrasample_source=logical(),alt_qc=character(),link_basis=character(),
    plausible=logical(),unresolved_nearby=logical())
  if(!nrow(loci)) return(empty)
  loc<-copy(loci);setorder(loc,CHROM,POS,locus_id)
  # Representative is used to SCREEN relationship; full original ALTs stay in internal RDS.
  r<-primary[, .SD[order(-fifelse(is.finite(RQ_SCORE),RQ_SCORE,0),-fifelse(is.finite(alt_support),alt_support,0),abs(POS-median(POS)),uid)][1L],by=locus_id]
  reps<-setNames(seq_len(nrow(r)),r$locus_id)
  out<-list();nlinks<-0L;align_n<-0L
  for(chr in unique(loc$CHROM)) {
    z<-loc[CHROM==chr];N<-nrow(z);if(N<2L) next
    for(i in seq_len(N-1L)) {
      j<-i+1L
      while(j<=N && z$POS[[j]]-z$POS[[i]]<=cfg$reconcile_window) {
        a<-z$locus_id[[i]];b<-z$locus_id[[j]]
        d<-z$POS[[j]]-z$POS[[i]]
        sim<-e$ra_size_sim(z$SVLEN[[i]],z$SVLEN[[j]])
        cl<-identical(z$CLASS_MODE[[i]],z$CLASS_MODE[[j]])
        ev<-identical(z$EVENT_MODE[[i]],z$EVENT_MODE[[j]])
        pa<-primary[locus_id==a];pb<-primary[locus_id==b]
        same_src<-FALSE
        shared_samples<-intersect(pa$sample_id,pb$sample_id)
        if(length(shared_samples)) for(sid in shared_samples) {
          aa<-pa[sample_id==sid];bb<-pb[sample_id==sid]
          if(length(intersect(e$ra_source_tokens(aa$source_ids[[1L]],aa$source_event_ids[[1L]]),
              e$ra_source_tokens(bb$source_ids[[1L]],bb$source_event_ids[[1L]])))) same_src<-TRUE
        }
        ar<-as.list(r[reps[[a]]]);br<-as.list(r[reps[[b]]])
        exact<-isTRUE(ar$seq_valid) && isTRUE(br$seq_valid) &&
          identical(toupper(ar$ALT_SEQUENCE),toupper(br$ALT_SEQUENCE))
        qc<-if(exact) 'EXACT_ALT' else 'NOT_EVALUATED'
        if(!exact && cfg$reconcile_align && is.finite(sim) &&
           sim>=cfg$reconcile_min_size_sim && cl &&
           isTRUE(ar$seq_valid) && isTRUE(br$seq_valid) &&
           align_n<cfg$reconcile_max_alignments) {
          align_n<-align_n+1L
          q<-tryCatch(e$retroatlas_alt_pairwise(ar$ALT_SEQUENCE,br$ALT_SEQUENCE,
            params=list(min_identity=cfg$duplicate_seq_identity,
              min_reciprocal_coverage=cfg$duplicate_seq_coverage,
              max_align_cells=cfg$primary_seq_max_cells)),error=function(err)NULL)
          qc<-if(is.null(q)) 'ALIGNMENT_UNAVAILABLE' else as.character(q$seq_qc[[1L]])
        }
        geometric<-d<=cfg$duplicate_geometry_window && cl && ev &&
          is.finite(sim) && sim>=cfg$duplicate_geometry_size_sim
        seq_ok<-qc %in% c('EXACT_ALT','SEQ_CONCORDANT')
        plausible<-same_src || seq_ok || geometric
        unresolved<-!plausible && cl && is.finite(sim) &&
          sim>=cfg$reconcile_min_size_sim &&
          !qc %in% c('SEQ_DISCORDANT','SEQ_PARTIAL')
        basis<-if(same_src) 'SHARED_EVENT_WITHIN_SAMPLE' else if(seq_ok)
          'ALT_SEQUENCE_CONCORDANT' else if(geometric) 'NEAR_GEOMETRY'
          else if(unresolved) 'CLOSE_STRUCTURAL_IDENTITY_UNRESOLVED' else 'DISTINCT_OR_WEAK'
        nlinks<-nlinks+1L
        out[[nlinks]]<-data.table(locus_a=a,locus_b=b,CHROM=chr,
          distance_bp=d,size_sim=sim,same_class=cl,same_event=ev,
          shared_intrasample_source=same_src,alt_qc=qc,link_basis=basis,
          plausible=plausible,unresolved_nearby=unresolved)
        j<-j+1L
      }
    }
  }
  if(length(out)) rbindlist(out) else empty
}
pm_neighbors <- function(locus_id,relations) {
  if(!nrow(relations)) return(relations[0L])
  relations[locus_a==locus_id | locus_b==locus_id]
}
pm_slice <- function(dt,chr,pos,win) {
  # Support per-chromosome caches to avoid O(VCF size * N loci) rescans.
  if(is.list(dt) && !is.data.frame(dt)) {
    s<-dt[[chr]]
    if(is.null(s)) return(data.table())
  } else {
    if(!nrow(dt)) return(dt)
    s<-dt[CHROM==chr]
  }
  if(!nrow(s)) return(s)
  lo<-findInterval(pos-win-1,s$POS)+1L
  hi<-findInterval(pos+win,s$POS)
  if(hi<lo) return(s[0L])
  s[seq.int(lo,hi)]
}
# CALL_OPT — cache softclips once per sample for BOTH source and target roles.
# Windows are derived solely from the original PRIMARY locus representatives.
# Callability controls whether a target locus is evaluated; it NEVER filters the
# clips retained from a source sample. We preserve original file order for
# source matching, and POS-sorted order for target matching (as CALL V3 did).
pm_clip_cache_windows <- function(reps,cfg) {
  w<-unique(reps[,.(CHROM=as.character(CHROM),
    start=pmax(0,round(as.numeric(POS))-as.numeric(cfg$clip_window)),
    end=round(as.numeric(POS))+as.numeric(cfg$clip_window))])
  w<-w[!is.na(CHROM)&nzchar(CHROM)&is.finite(start)&is.finite(end)]
  if(!nrow(w)) stop('[POP][CLIP_CACHE] No valid PRIMARY representative windows')
  setorder(w,CHROM,start,end)
  # Merge overlapping / immediately adjacent windows without shrinking the
  # inclusive boundaries. Keep the exact +/-clip_window used by D1.
  runs<-vector('list',length(unique(w$CHROM)));idx<-0L
  for(chr in unique(w$CHROM)) {
    z<-w[CHROM==chr]
    lo<-z$start[[1L]];hi<-z$end[[1L]]
    for(i in seq_len(nrow(z))) {
      if(i==1L) next
      a<-z$start[[i]];b<-z$end[[i]]
      if(a<=hi+1) hi<-max(hi,b) else {
        idx<-idx+1L;runs[[idx]]<-data.table(CHROM=chr,start=lo,end=hi)
        lo<-a;hi<-b
      }
    }
    idx<-idx+1L;runs[[idx]]<-data.table(CHROM=chr,start=lo,end=hi)
  }
  ans<-rbindlist(runs[seq_len(idx)])
  setorder(ans,CHROM,start)
  ans
}
pm_clip_cache_select <- function(dt,windows) {
  # Vectorized interval search: one binary-search index per clip, no cross
  # product between all clips and all representatives. No QC-based filtering.
  if(!nrow(dt)) return(dt)
  keep<-rep(FALSE,nrow(dt))
  for(chr in intersect(unique(dt$CHROM),unique(windows$CHROM))) {
    ii<-which(dt$CHROM==chr)
    z<-windows[CHROM==chr]
    k<-findInterval(dt$POS[ii],z$start)
    take<-k>0L
    valid<-which(take)
    if(length(valid)) take[valid]<-dt$POS[ii[valid]]<=z$end[k[valid]]
    keep[ii]<-take
  }
  dt[keep]
}
pm_clip_cache_prepare <- function(row,windows,e,cfg,cache_dir) {
  sid<-as.character(row$sample_id[[1L]])
  path<-as.character(row$softclips[[1L]])
  if(grepl('\\.(fa|fasta)(\\.gz)?$',path,ignore.case=TRUE)) {
    message('[POP][CLIP_CACHE_FALLBACK] sample=',sid,
      ' | format=FASTA | legacy region reader retained')
    return(NULL)
  }
  t0<-Sys.time()
  # load_softclip_table reads the gzip exactly once, preserving original rows.
  raw<-as.data.table(e$load_softclip_table(path))
  needed<-c('CLIP_ID','CHR','POS0','POS1','SOFTCLIP_LEN','SEQUENCE')
  if(!all(needed %in% names(raw))) stop('[POP][CLIP_CACHE_INPUT] Missing fields: ',
    paste(setdiff(needed,names(raw)),collapse=','),' sample=',sid)
  read_rows<-nrow(raw)
  raw[,`:=`(CHROM=e$ra_chr(CHR),POS=suppressWarnings(as.numeric(POS1)))]
  raw<-raw[is.finite(POS)]
  raw[,RA_INPUT_ORDER:=seq_len(.N)]
  filtered<-pm_clip_cache_select(raw,windows)
  # Normalize to the exact fields used by the D1 clip module, plus internal
  # sortable keys. Avoid keeping unrelated input columns in each fork.
  filtered<-filtered[,c(needed,'CHROM','POS','RA_INPUT_ORDER'),with=FALSE]
  rm(raw)
  # Preserve the stable tie order from the original file for target selection.
  target<-copy(filtered);setorder(target,CHROM,POS,RA_INPUT_ORDER)
  source<-copy(filtered);setorder(source,CHROM,POS,RA_INPUT_ORDER)
  source_by_chr<-split(source,source$CHROM)
  target_by_chr<-split(target,target$CHROM)
  ans<-list(source=source_by_chr,target=target_by_chr,
    read_rows=read_rows,retained_rows=nrow(filtered))
  if(isTRUE(cfg$source_clip_cache_save_rds)) {
    # Number-based filenames avoid exposing sample identifiers in paths.
    cache_file<-file.path(cache_dir,sprintf('softclip_cache_%03d.rds',as.integer(row$RA_CACHE_INDEX[[1L]])))
    saveRDS(ans,cache_file,compress=FALSE)
  }
  bytes<-as.numeric(object.size(ans))
  message('[POP][CLIP_CACHE_READY] sample=',sid,' | read_rows=',read_rows,
    ' | retained=',nrow(filtered),' | cache_MB=',round(bytes/1024^2,2),
    ' | seconds=',round(as.numeric(difftime(Sys.time(),t0,units='secs')),2),
    ' | source_target_cached=TRUE')
  ans
}
pm_clip_source_slice <- function(cache,chr,pos,window) {
  # The D1 source reader returned rows in ORIGINAL gzip order, whereas
  # target preloading returned coordinate-sorted rows. Restore this order.
  z<-pm_slice(cache,chr,pos,window)
  if(!nrow(z)) return(data.frame(CLIP_ID=character(),CHR=character(),
    POS0=numeric(),POS1=numeric(),SOFTCLIP_LEN=numeric(),SEQUENCE=character(),
    stringsAsFactors=FALSE))
  setorder(z,RA_INPUT_ORDER)
  result<-as.data.frame(z[,.(CLIP_ID,CHR=CHROM,POS0,POS1,SOFTCLIP_LEN,SEQUENCE)])
  result$POS0<-suppressWarnings(as.numeric(result$POS0))
  result$POS1<-suppressWarnings(as.numeric(result$POS1))
  result$SOFTCLIP_LEN<-suppressWarnings(as.numeric(result$SOFTCLIP_LEN))
  result$SEQUENCE<-toupper(as.character(result$SEQUENCE))
  result
}
pm_clip_cache_preflight <- function(cfg,e) {
  dt<-data.table(CLIP_ID=c('A','B','C','D','E','F'),
    CHR=rep('chr1',6),POS0=c(999,400,1499,2000,998,1498),
    POS1=c(1000,401,1500,2001,999,1499),
    SOFTCLIP_LEN=rep(60,6),SEQUENCE=rep(strrep('A',60),6))
  dt[,`:=`(CHROM='chr1',POS=as.numeric(POS1),RA_INPUT_ORDER=seq_len(.N))]
  w<-data.table(CHROM='chr1',start=500,end=1500)
  z<-pm_clip_cache_select(dt,w)
  stopifnot(identical(z$CLIP_ID,c('A','C','E','F')))
  src<-copy(z);setorder(src,CHROM,POS,RA_INPUT_ORDER)
  cached<-pm_clip_source_slice(split(src,src$CHROM),'chr1',1000,500)
  reference<-as.data.frame(dt[is.finite(POS1) & abs(POS1-1000)<=500,
    .(CLIP_ID,CHR,POS0,POS1,SOFTCLIP_LEN,SEQUENCE)])
  stopifnot(identical(as.character(cached$CLIP_ID),as.character(reference$CLIP_ID)))
  message('[POP][CLIP_CACHE_PREFLIGHT] Inclusive boundaries and source-order parity OK')
  invisible(TRUE)
}
pm_empty_state <- function(lid,sid) data.table(locus_id=lid,sample_id=sid,
  state='NOT_CALLABLE',reason='NOT_EVALUATED',rescue_type=NA_character_,
  primary_uid=NA_character_,source_uid=NA_character_,source_sample=NA_character_,
  GT=NA_character_,dosage=NA_integer_,genotype_qc=NA_character_,
  ploidy=NA_integer_,sex_region=NA_character_,
  callable=FALSE,callability_status=NA_character_,frac_ge20=NA_real_,mean_depth=NA_real_,
  quality_gate_skipped=FALSE,
  vcf_status=NA_character_,vcf_callers=NA_character_,vcf_evidence_ids='',
  clip_status=NA_character_,clip_good_ids='',clip_n_good=NA_integer_,
  clip_error_reason=NA_character_,
  fragment_status=NA_character_,fragment_ids='',
  linked_loci='',linked_primary_uids='',evidence_level='NONE',
  assumption_one_alt=FALSE,possible_duplicate=FALSE)
pm_read_sample_vcf <- function(row,e,cfg) {
  out<-lapply(c('cutesv','severus','sniffles2'),function(caller) {
    e$ra_read_vcf(row[[paste0('vcf_',caller)]][[1L]],row$sample_id[[1L]],caller,cfg)
  })
  v<-rbindlist(lapply(out,`[[`,'table'),fill=TRUE)
  if(nrow(v)) setorder(v,CHROM,POS)
  list(table=v,audit=rbindlist(lapply(out,`[[`,'audit'),fill=TRUE))
}
pm_seq_competitors <- function(target_primary,source,cfg) {
  if(!nrow(target_primary)) return(NULL)
  z<-target_primary[CHROM==source$CHROM & is.finite(POS) &
    abs(POS-source$POS)<=cfg$rescue_window & seq_valid==TRUE & !is.na(ALT_SEQUENCE),ALT_SEQUENCE]
  z<-setdiff(unique(as.character(z)),source$ALT_SEQUENCE)
  if(length(z)>6L) z<-z[seq_len(6L)]
  if(length(z)) z else NULL
}
pm_positive <- function(row) row$state %in% c('PRIMARY','RESCUE_GT','RESCUE_ASSUMED')
# Sample x locus state. A positive evidence call cannot bypass a plausible
# already-existing PRIMARY representation in a related locus.
pm_one_state <- function(locus,source,source_sample,target_row,target_primary,
                         edges,quality,vcf,frags,clips_target,source_clip_path,e,cfg,
                         missing_vcf,source_clip_cache=NULL) {
  lid<-locus$locus_id[[1L]];sid<-target_row$sample_id[[1L]]
  ret<-pm_empty_state(lid,sid)
  ret[,`:=`(source_uid=source$uid,source_sample=source$sample_id,
    ploidy=as.integer(e$sa_ploidy(e$sa_region(locus$CHROM[[1L]],locus$POS[[1L]],
      cfg$reference_build),target_row$sex_karyotype[[1L]])),
    sex_region=as.character(e$sa_region(locus$CHROM[[1L]],locus$POS[[1L]],cfg$reference_build)),
    callable=isTRUE(quality$callable[[1L]]),
    callability_status=as.character(quality$quality_status[[1L]]),
    frac_ge20=as.numeric(quality$frac_ge20[[1L]]),
    mean_depth=as.numeric(quality$mean_depth[[1L]]))]
  if(is.na(ret$ploidy[[1L]]) || ret$ploidy[[1L]]==0L ||
     (!cfg$include_sex_chromosomes && ret$sex_region[[1L]]!='AUTOSOME')) {
    ret[,`:=`(state='SEX_EXCLUDED',reason='SEX_PLOIDY_NOT_EVALUABLE')]
    return(ret)
  }
  # PRE-RESCUE REGIONAL QUALITY GATE. The same ra_quality() quality BED
  # decision used for negative callability is applied before ALL positive
  # secondary evidence calls. Fail closed on missing/unpaired BED intervals.
  # Direct PRIMARY rows are handled separately and are not altered here.
  if(isTRUE(cfg$rescue_require_callable) && !isTRUE(ret$callable[[1L]])) {
    qc_status<-ret$callability_status[[1L]]
    if(is.na(qc_status) || !nzchar(qc_status)) qc_status<-'QUALITY_STATUS_UNAVAILABLE'
    ret[,`:=`(state='NOT_CALLABLE',
      reason=paste0('PRE_RESCUE_QUALITY_GATE:',qc_status),
      evidence_level='REGIONAL_QUALITY_GATE',
      quality_gate_skipped=TRUE,
      clip_status='NOT_EVALUATED_QC_GATE')]
    return(ret)
  }
  # Already has a PRIMARY at this locus: the caller should not invoke this function.
  own<-target_primary[CHROM==locus$CHROM[[1L]] &
    abs(POS-locus$POS[[1L]])<=cfg$reconcile_window]
  if(nrow(edges)) {
    z<-edges[plausible==TRUE | unresolved_nearby==TRUE]
    related<-unique(ifelse(z$locus_a==lid,z$locus_b,z$locus_a))
    if(nrow(own)) own<-own[locus_id %in% related] # safe if no edge
  } else own<-own[0L]
  if(nrow(own)) {
    ret[,`:=`(state='REVIEW',reason='NEARBY_REPRESENTED_PRIMARY_OR_ALLELE',
      possible_duplicate=TRUE,evidence_level='PRIMARY_NEARBY',
      linked_loci=paste(sort(unique(own$locus_id)),collapse=';'),
      linked_primary_uids=paste(sort(unique(own$uid)),collapse=';'))]
    return(ret)
  }
  # Ambiguous nearby locus whose identity has not been tested: do not make a
  # confident new rescue without a unique locus assignment.
  local_links<-edges[plausible==TRUE | unresolved_nearby==TRUE]
  unresolved_links<-local_links[unresolved_nearby==TRUE]
  # Source is a directly reconstructed ALT from an ACTUAL sample, never a consensus.
  zvcf<-pm_slice(vcf,source$CHROM,source$POS,cfg$vcf_window)
  v<-e$ra_vcf_evaluate(source,zvcf,cfg,ret$ploidy[[1L]])
  f<-e$ra_fragment_evaluate(source,pm_slice(frags,source$CHROM,source$POS,cfg$fragment_window),cfg)
  ret[,`:=`(vcf_status=v$status,vcf_callers=v$alt_callers,
    vcf_evidence_ids=v$selected_record_ids,fragment_status=f$status,
    fragment_ids=f$ids)]
  clip<-list(resolution='NOT_EVALUATED_BY_POLICY',n_good_alt_clips=NA_integer_,details=NULL)
  do_clip<-cfg$clip_enabled &&
    (cfg$clip_evaluation=='all' || v$status!='VCF_ALT_GT_RESCUE')
  if(do_clip) {
    competitor_alts<-pm_seq_competitors(target_primary,source,cfg)
    # IMPORTANT FIX: the softclip module uses base data.frame column subsetting.
    # pm_slice() returns data.table; passing it through unchanged caused failures
    # for nonempty windows and made technical errors look like biological REVIEW.
    nearby_clips<-if(is.list(clips_target) && !is.data.frame(clips_target))
      pm_slice(clips_target,source$CHROM,source$POS,cfg$clip_window) else clips_target
    if(!is.null(nearby_clips)) {
      if(!is.data.frame(nearby_clips))
        stop('[POP][SOFTCLIP_INPUT] Unexpected clip-table type for sample=',sid,
          ' locus=',lid,' class=',paste(class(nearby_clips),collapse='/'))
      if(nrow(nearby_clips)) {
        needed<-c('CLIP_ID','CHR','POS0','POS1','SOFTCLIP_LEN','SEQUENCE')
        if(!all(needed %in% names(nearby_clips)))
          stop('[POP][SOFTCLIP_INPUT] Missing columns for sample=',sid,
            ' locus=',lid,': ',paste(setdiff(needed,names(nearby_clips)),collapse=','))
      }
      # A small per-locus slice, not the full cached TSV: cheap and safe.
      nearby_clips<-as.data.frame(nearby_clips)
    }
    if(!is.null(nearby_clips) && nrow(nearby_clips)==0L) {
      clip<-list(resolution='NO_CLIP_SUPPORT',reason='EMPTY_TARGET_WINDOW',
        n_good_alt_clips=0L,details=NULL)
    } else {
      source_tbl<-if(!is.null(source_clip_cache))
        pm_clip_source_slice(source_clip_cache,source$CHROM,source$POS,cfg$clip_window) else NULL
      clip<-e$ra_softclip_evaluate(source,source_clip_path,target_row$softclips[[1L]],
        source_tbl,nearby_clips,cfg,competitors=competitor_alts)
    }
  }
  clip_ids<-''
  if(!is.null(clip$details) && !is.null(clip$details$clips) && nrow(clip$details$clips)) {
    cd<-as.data.table(clip$details$clips)
    if(all(c('CLIP_ID','alt_match') %in% names(cd)))
      clip_ids<-paste(unique(as.character(cd$CLIP_ID[cd$alt_match %in% TRUE])),collapse=';')
  }
  clip_error<-if(identical(clip$resolution,'CLIP_ERROR')) {
    z<-as.character(clip$reason)
    if(!length(z) || all(is.na(z))) z<-'ERROR_MESSAGE_NOT_AVAILABLE'
    gsub('[\\r\\n\\t]+',' ',paste(z,collapse='; '))
  } else NA_character_
  ret[,`:=`(clip_status=clip$resolution,clip_good_ids=clip_ids,
    clip_n_good=as.integer(clip$n_good_alt_clips),clip_error_reason=clip_error)]
  # Fail before any REVIEW or AF assignment: a technical error is not biological evidence.
  if(identical(clip$resolution,'CLIP_ERROR')) {
    context<-paste0('sample=',sid,' locus=',lid,' source=',source$uid,
      ' chr=',source$CHROM,' pos=',source$POS,' reason=',clip_error)
    if(identical(cfg$clip_error_policy,'abort'))
      stop('[POP][SOFTCLIP_FATAL] ',context,call.=FALSE)
    ret[,`:=`(state='NOT_CALLABLE',reason='SOFTCLIP_TECHNICAL_ERROR',
      evidence_level='TECHNICAL_ERROR')]
    return(ret)
  }
  # A nearby unresolved locus is not proof of identity, but blocks assignment
  # of weak secondary evidence without resolving the competing structure.
  has_local_evidence <- v$status!='NO_VCF_SUPPORT' ||
    clip$resolution %in% c('RESCUE_CANDIDATE','REVIEW') ||
    f$status %in% c('FRAGMENT_RESCUE_CANDIDATE','FRAGMENT_AMBIGUOUS')
  if(nrow(unresolved_links) && has_local_evidence) {
    ret[,`:=`(state='REVIEW',reason='COMPETING_LINKED_LOCUS_UNRESOLVED',
      linked_loci=paste(unique(ifelse(unresolved_links$locus_a==lid,
        unresolved_links$locus_b,unresolved_links$locus_a)),collapse=';'),
      evidence_level='AMBIGUOUS_LOCUS_LINK')]
  } else if(v$status=='VCF_ALT_GT_RESCUE') {
    ret[,`:=`(state='RESCUE_GT',reason=v$reason,rescue_type='VCF',
      evidence_level=if(v$n_alt_callers>=2L) 'VCF_MULTICALLER' else 'VCF_SINGLE_CALLER',
      GT=v$gt,dosage=as.integer(v$dosage),genotype_qc='VCF_OBSERVED')]
  } else if(clip$resolution=='RESCUE_CANDIDATE') {
    ret[,`:=`(state='RESCUE_ASSUMED',reason=paste0('SOFTCLIP_',clip$reason),
      rescue_type='SOFTCLIP',assumption_one_alt=TRUE,evidence_level='SOFTCLIP')]
  } else if(f$status=='FRAGMENT_RESCUE_CANDIDATE') {
    ret[,`:=`(state='RESCUE_ASSUMED',reason='FRAGMENT_SUM_WITH_RESOLVED_ANCHOR',
      rescue_type='FRAGMENT',assumption_one_alt=TRUE,evidence_level='FRAGMENT')]
  } else if(v$status %in% c('VCF_0_0_WITH_ALT_SUPPORT','VCF_UNTYPED_ALT_SUPPORT',
      'VCF_SINGLE_CALLER_REVIEW','VCF_ALT_CONTRADICTS_PRIMARY',
      'VCF_NEAR_SIZE_DISCORDANT','VCF_NEAR_NO_ALT') ||
      clip$resolution %in% c('REVIEW') ||
      f$status=='FRAGMENT_AMBIGUOUS') {
    ret[,`:=`(state='REVIEW',reason='CONFLICTING_OR_INCOMPLETE_LOCAL_EVIDENCE',
      evidence_level='SECONDARY_AMBIGUOUS')]
  } else if(!isTRUE(quality$callable[[1L]])) {
    ret[,`:=`(state='NOT_CALLABLE',reason=as.character(quality$quality_status[[1L]]))]
  } else if(cfg$require_all_vcf_for_differential && missing_vcf) {
    ret[,`:=`(state='NOT_CALLABLE',reason='REQUIRED_VCF_SOURCE_NOT_AVAILABLE')]
  } else if(!cfg$clip_enabled) {
    ret[,`:=`(state='NOT_CALLABLE',reason='SOFTCLIP_DISABLED_NO_NEGATIVE_EVIDENCE') ]
  } else if(cfg$clip_enabled && clip$resolution %in% c('CLIP_ERROR','NO_CLIP_EVALUATION')) {
    ret[,`:=`(state='NOT_CALLABLE',reason='SOFTCLIP_EVALUATION_UNAVAILABLE')]
  } else {
    ret[,`:=`(state='REF_CALLABLE',reason='CALLABLE_AND_NO_SUPPORTING_ALT_EVIDENCE')]
  }
  ret
}
pm_resolve_reuse <- function(states,links) {
  d<-copy(states)
  # A VCF record/softclip read is allowed to support at most ONE locus in a sample.
  claims<-list();k<-0L
  idx<-which(d$state %in% c('RESCUE_GT','RESCUE_ASSUMED'))
  for(i in idx) {
    kind<-d$rescue_type[[i]]
    raw<-if(kind=='VCF') d$vcf_evidence_ids[[i]] else if(kind=='SOFTCLIP')
      d$clip_good_ids[[i]] else d$fragment_ids[[i]]
    if(is.na(raw) || !nzchar(raw)) next
    ids<-strsplit(raw,';',fixed=TRUE)[[1L]]
    ids<-unique(ids[nzchar(ids)])
    if(!length(ids)) next
    k<-k+1L
    claims[[k]]<-data.table(row_id=i,locus_id=d$locus_id[[i]],
      sample_id=d$sample_id[[i]],kind=kind,evidence_id=ids)
  }
  ledger<-if(length(claims)) rbindlist(claims) else data.table()
  reused<-data.table()
  if(nrow(ledger)) {
    reused<-ledger[,.(n_loci=uniqueN(locus_id),
      loci=paste(sort(unique(locus_id)),collapse=';')),
      by=.(sample_id,kind,evidence_id)][n_loci>1L]
    if(nrow(reused)) {
      hits<-unique(ledger[reused,on=.(sample_id,kind,evidence_id)]$row_id)
      d[hits,`:=`(state='REVIEW',reason='EVIDENCE_REUSED_BETWEEN_POPULATION_LOCI',
        evidence_level='AMBIGUOUS_SHARED_EVIDENCE',GT=NA_character_,
        dosage=NA_integer_,assumption_one_alt=FALSE)]
    }
  }
  # Two rescues cannot both increase the count in linked, plausibly same loci.
  # No automatic promotion of one candidate; competing evidence remains REVIEW.
  blocked<-list()
  if(nrow(links)) {
    z<-links[plausible==TRUE]
    for(i in seq_len(nrow(z))) {
      a<-z$locus_a[[i]];b<-z$locus_b[[i]]
      both<-d[locus_id %in% c(a,b) & state %in% c('RESCUE_GT','RESCUE_ASSUMED')]
      if(!nrow(both)) next
      groups<-both[,.(N=uniqueN(locus_id)),by=sample_id][N>1L]
      if(nrow(groups)) {
        ii<-which(d$sample_id %in% groups$sample_id &
            d$locus_id %in% c(a,b) & d$state %in% c('RESCUE_GT','RESCUE_ASSUMED'))
        if(length(ii)) blocked[[length(blocked)+1L]]<-ii
      }
    }
  }
  lost<-if(length(blocked)) sort(unique(unlist(blocked,use.names=FALSE))) else integer()
  if(length(lost)) d[lost,`:=`(state='REVIEW',
    reason='RESCUE_COMPETITION_ACROSS_PLAUSIBLY_LINKED_LOCI',
    evidence_level='AMBIGUOUS_LINKED_LOCI',GT=NA_character_,
    dosage=NA_integer_,assumption_one_alt=FALSE)]
  list(states=d,reused_evidence=reused,
    competing_rescues=data.table(row_id=lost,locus_id=d$locus_id[lost],
      sample_id=d$sample_id[lost]))
}
pm_frequency <- function(states,loci,metadata,e,cfg) {
  d<-copy(states)
  d[,`:=`(AC_PRIMARY=0L,AN_PRIMARY=0L,AC_GT=0L,AN_GT=0L,
    AC_ADJUSTED=0L,AN_ADJUSTED=0L)]
  eligible<-which(!is.na(d$ploidy) & d$ploidy>0L)
  for(i in eligible) {
    state<-d$state[[i]];pl<-as.integer(d$ploidy[[i]])
    if(state=='PRIMARY') {
      if(is.finite(d$dosage[[i]])) {
        dose<-as.integer(d$dosage[[i]])
        d[i,`:=`(AC_PRIMARY=dose,AN_PRIMARY=pl,AC_GT=dose,AN_GT=pl,
          AC_ADJUSTED=dose,AN_ADJUSTED=pl)]
      }
    } else if(state=='RESCUE_GT' && is.finite(d$dosage[[i]])) {
      dose<-as.integer(d$dosage[[i]])
      d[i,`:=`(AC_GT=dose,AN_GT=pl,AC_ADJUSTED=dose,AN_ADJUSTED=pl)]
    } else if(state=='RESCUE_ASSUMED') {
      d[i,`:=`(AC_ADJUSTED=1L,AN_ADJUSTED=pl)]
    } else if(state=='REF_CALLABLE') {
      d[i,`:=`(AN_PRIMARY=pl,AN_GT=pl,AN_ADJUSTED=pl)]
    }
  }
  if(any(d$AC_ADJUSTED>d$AN_ADJUSTED | d$AC_GT>d$AN_GT | d$AC_PRIMARY>d$AN_PRIMARY))
    stop('Allele-count invariant violated (AC > AN)')
  sm<-d[,.(N_PRIMARY=sum(state=='PRIMARY'),N_RESCUE_GT=sum(state=='RESCUE_GT'),
    N_RESCUE_ASSUMED=sum(state=='RESCUE_ASSUMED'),N_REVIEW=sum(state=='REVIEW'),
    N_NOT_CALLABLE=sum(state=='NOT_CALLABLE'),N_REF_CALLABLE=sum(state=='REF_CALLABLE'),
    N_SEX_EXCLUDED=sum(state=='SEX_EXCLUDED'),
    N_PRIMARY_HOM_ALT=sum(state=='PRIMARY' & ploidy==2L & dosage==2L,na.rm=TRUE),
    N_PRIMARY_HET=sum(state=='PRIMARY' & ploidy==2L & dosage==1L,na.rm=TRUE),
    AC_PRIMARY=sum(AC_PRIMARY),AN_PRIMARY=sum(AN_PRIMARY),
    AC_GT=sum(AC_GT),AN_GT=sum(AN_GT),
    AC_ADJUSTED=sum(AC_ADJUSTED),AN_ADJUSTED=sum(AN_ADJUSTED)),by=locus_id]
  sm[,`:=`(AF_PRIMARY=ifelse(AN_PRIMARY>0,AC_PRIMARY/AN_PRIMARY,NA_real_),
    AF_GT=ifelse(AN_GT>0,AC_GT/AN_GT,NA_real_),
    AF_RESCUE_ADJUSTED=ifelse(AN_ADJUSTED>0,AC_ADJUSTED/AN_ADJUSTED,NA_real_))]
  # N_PRIMARY is recomputed from exclusive locus x sample states; avoid .x/.y.
  base_loci<-copy(loci)
  base_loci[,N_PRIMARY:=NULL]
  q<-merge(base_loci,sm,by='locus_id',all.x=TRUE,sort=FALSE)
  q[,N_SAMPLES:=as.integer(nrow(metadata))]
  # A singleton has a confirmed *absence in every other eligible sample*
  # only when no rescue/review/noncallability remains. In the public catalogue
  # it is always just another locus with observed allele frequency.
  q[,locus_interpretation:=fifelse(N_PRIMARY>=2L,'SHARED_PRIMARY',
    fifelse(N_PRIMARY==1L & N_RESCUE_GT==0L & N_RESCUE_ASSUMED==0L &
      N_REVIEW==0L & N_NOT_CALLABLE==0L,
      'CONFIRMED_SINGLETON_AMONG_EVALUABLE','SINGLETON_OR_UNCERTAIN'))]
  # Population group distributions are computed only from manifest metadata.
  m<-metadata[,.(sample_id,sex,sex_karyotype,origin)]
  z<-merge(d,m,by='sample_id',all.x=TRUE,sort=FALSE)
  if(anyNA(z$sex)) stop('Missing metadata in frequency input')
  private_strata<-list()
  for(group in c('sex','origin')) {
    v<-z[,.(N_PARTICIPANTS=.N,
      N_PRIMARY=sum(state=='PRIMARY'),N_RESCUE=sum(state %in% c('RESCUE_GT','RESCUE_ASSUMED')),
      AC_GT=sum(AC_GT),AN_GT=sum(AN_GT),
      AC_ADJUSTED=sum(AC_ADJUSTED),AN_ADJUSTED=sum(AN_ADJUSTED)),
      by=c('locus_id',group)]
    setnames(v,group,'stratum')
    v[,stratum_type:=group]
    v[,`:=`(AF_GT=ifelse(AN_GT>0,AC_GT/AN_GT,NA_real_),
      AF_ADJUSTED=ifelse(AN_ADJUSTED>0,AC_ADJUSTED/AN_ADJUSTED,NA_real_))]
    private_strata[[length(private_strata)+1L]]<-v
  }
  strata<-rbindlist(private_strata,fill=TRUE)
  # Rare stratum counts are not exported by default even without sample IDs.
  public_strata<-strata[N_PARTICIPANTS>=cfg$min_public_stratum_n]
  list(summary=q,private_states=d,private_strata=strata,public_strata=public_strata)
}
pm_catalog <- function(summary,primary,metadata,public_strata,cfg) {
  # PRIMARY-only observations: rescues do not modify POS/SVLEN/ALT vectors.
  p<-merge(copy(primary),metadata[,.(sample_id,sex,origin)],
    by='sample_id',all.x=TRUE,sort=FALSE)
  set.seed(as.integer(cfg$seed))
  p[,random_key:=runif(.N)]
  setorder(p,locus_id,random_key)
  # All individual measurements within an observation are shuffled TOGETHER.
  # Independently shuffling POS, SVLEN, ALT would destroy scientific pairing.
  counts<-function(x) {
    z<-as.character(x);z<-z[!is.na(z) & nzchar(z) & z!='NA']
    if(!length(z)) return('')
    t<-sort(table(z),decreasing=TRUE)
    paste(paste(names(t),as.integer(t),sep=':'),collapse=';')
  }
  agg<-p[,.(
    N_PRIMARY_OBSERVATIONS=.N,
    POS_MEDIAN=stats::median(POS),POS_MIN=if(.N) min(POS) else NA_real_,POS_MAX=if(.N) max(POS) else NA_real_,
    SVLEN_MEDIAN=stats::median(SVLEN),SVLEN_MIN=if(.N) min(SVLEN) else NA_real_,SVLEN_MAX=if(.N) max(SVLEN) else NA_real_,
    CLASS_COUNTS=counts(event_class),EVENT_COUNTS=counts(event),
    COMPLETENESS_COUNTS=counts(completeness),ELEMENT_TYPE_COUNTS=counts(element_type),
    POLYA_COUNTS=counts(polyA),POLYA_LENGTH_COUNTS=counts(polyA_length),
    TSD_COUNTS=counts(tsd_detected),TSD_LENGTH_COUNTS=counts(tsd_length),
    N_PRIMARY_DIPLOID_HOM_ALT=sum(ploidy==2L & GT %in% c('1/1','1|1'),na.rm=TRUE),
    POS_VALUES=list(as.numeric(POS)),SVLEN_VALUES=list(as.numeric(SVLEN)),
    ALT_SEQUENCES=list(as.character(ALT_SEQUENCE)),
    CLASS_VALUES=list(as.character(event_class)),EVENT_VALUES=list(as.character(event)),
    COMPLETENESS_VALUES=list(as.character(completeness)),
    ELEMENT_TYPE_VALUES=list(as.character(element_type)),
    POLYA_VALUES=list(as.character(polyA)),POLYA_LENGTH_VALUES=list(as.numeric(polyA_length)),
    TSD_VALUES=list(as.character(tsd_detected)),TSD_LENGTH_VALUES=list(as.numeric(tsd_length)),
    SEX_VALUES=list(as.character(sex)),ORIGIN_VALUES=list(as.character(origin))
  ),by=locus_id]
  selected<-c('locus_id','CHROM','POS','SVLEN','N_PRIMARY','N_SAMPLES','Q_MIN',
    'locus_interpretation','N_RESCUE_GT','N_RESCUE_ASSUMED',
    'N_REVIEW','N_NOT_CALLABLE','N_REF_CALLABLE','N_SEX_EXCLUDED',
    'N_PRIMARY_HOM_ALT','N_PRIMARY_HET','AC_PRIMARY','AN_PRIMARY','AF_PRIMARY',
    'AC_GT','AN_GT','AF_GT','AC_ADJUSTED','AN_ADJUSTED','AF_RESCUE_ADJUSTED')
  s<-summary[,..selected]
  out<-merge(s,agg,by='locus_id',all.x=TRUE,sort=FALSE)
  if(cfg$public_detail=='summary') {
    # Neither ALT sequences nor raw per-individual coordinate vectors are public.
    remove<-c('POS_VALUES','SVLEN_VALUES','ALT_SEQUENCES','CLASS_VALUES','EVENT_VALUES',
      'COMPLETENESS_VALUES','ELEMENT_TYPE_VALUES','POLYA_VALUES',
      'POLYA_LENGTH_VALUES','TSD_VALUES','TSD_LENGTH_VALUES',
      'SEX_VALUES','ORIGIN_VALUES')
    out[,(remove):=NULL]
  }
  setorder(out,CHROM,POS,locus_id)
  list(format_version='RetroAtlas_Population_Catalog_v1',
    metadata=list(catalog_id=cfg$catalog_id,reference_build=cfg$reference_build,
      n_samples=nrow(metadata),n_loci=nrow(out),seed=cfg$seed,
      detail=cfg$public_detail,
      privacy='DEIDENTIFIED_NOT_ANONYMOUS',
      warning='Rare loci, precise coordinates and especially ALT sequences may allow reidentification; do not treat this catalogue as anonymous'),
    clusters=out,strata=public_strata)
}
pm_tsv <- function(x,path) {
  z<-copy(as.data.table(x))
  # fwrite(NULL) neither creates nor truncates a file. Empty audits are outputs too.
  if(!ncol(z)) {
    con<-file(path,'wt');close(con)
    return(invisible(NULL))
  }
  for(n in names(z)) if(is.list(z[[n]])) z[,(n):=vapply(get(n),function(q)
    paste(as.character(q),collapse=';'),character(1))]
  fwrite(z,path,sep='\t',na='NA')
}
# Prove that failing the regional gate short-circuits ALL secondary
# evidence APIs. The mock environment deliberately provides no VCF,
# fragment or softclip methods: reaching any of them must fail the test.
pm_regional_gate_preflight <- function() {
  cfg<-pm_cfg_defaults()
  mock<-new.env(parent=emptyenv())
  mock$sa_region<-function(chr,pos,build) 'AUTOSOME'
  mock$sa_ploidy<-function(region,karyotype) 2L
  loc<-data.table(locus_id='QC_PROBE_LOCUS',CHROM='chr1',POS=1000)
  src<-list(uid='QC_PROBE_SOURCE',sample_id='QC_SOURCE',CHROM='chr1',POS=1000)
  target<-data.table(sample_id='QC_TARGET',sex_karyotype='XX')
  q<-data.table(callable=FALSE,quality_status='QUALITY_NO_OVERLAP',
    frac_ge20=NA_real_,mean_depth=NA_real_)
  out<-pm_one_state(loc,src,NULL,target,data.table(),data.table(),q,
    NULL,NULL,NULL,NULL,mock,cfg,TRUE)
  if(nrow(out)!=1L || out$state[[1L]]!='NOT_CALLABLE' ||
     !isTRUE(out$quality_gate_skipped[[1L]]) ||
     out$clip_status[[1L]]!='NOT_EVALUATED_QC_GATE' ||
     !is.na(out$vcf_status[[1L]]) || !is.na(out$fragment_status[[1L]]))
    stop('[POP][QC_GATE_SELFTEST_FAILED] Bad-quality locus was evaluated/rescued')
  message('[POP][QC_GATE_SELFTEST] Noncallable locus bypasses VCF, fragments and softclip')
  invisible(TRUE)
}

pm_audit_selfcheck <- function() {
  # Check the structure of the state model; no biological observations required.
  d<-pm_empty_state('LocusAtlas_0000001','TEST')
  stopifnot(nrow(d)==1L, d$state[[1L]]=='NOT_CALLABLE',
    identical(pm_cfg_defaults()$primary_window,500),
    pm_cfg_defaults()$public_detail=='summary')
  invisible(TRUE)
}
# Dynamic chunks balance expensive softclip loci against fast PRIMARY/QC paths.
# Forks inherit the read-only caches; no cluster export or nested thread pools.
# Results always retain input order. A failed worker aborts instead of losing loci.
pm_evidence_parallel <- function(indices,fun,cores=1L,scheduler='dynamic',chunk_loci=4L) {
  if(!length(indices)) return(list())
  cores<-max(1L,min(as.integer(cores),length(indices)))
  if(.Platform$OS.type!='unix') cores<-1L
  if(scheduler=='static') {
    chunks<-split(indices,rep(seq_len(cores),length.out=length(indices)))
  } else {
    size<-max(1L,min(as.integer(chunk_loci),ceiling(length(indices)/(2L*cores))))
    chunks<-split(indices,(seq_along(indices)-1L)%/%size)
  }
  task<-function(ix) {
    data.table::setDTthreads(1L)
    start<-proc.time()[['elapsed']]
    list(values=lapply(ix,fun),seconds=proc.time()[['elapsed']]-start)
  }
  start<-proc.time()[['elapsed']]
  parts<-if(cores==1L) lapply(chunks,task) else parallel::mclapply(chunks,task,
    mc.cores=cores,mc.preschedule=scheduler=='static',mc.set.seed=FALSE)
  bad<-vapply(seq_along(parts),function(i) is.null(parts[[i]]) ||
    inherits(parts[[i]],'try-error') || length(parts[[i]]$values)!=length(chunks[[i]]),logical(1))
  if(any(bad)) stop('[POP][WORKER_FAILURE] Evidence worker failed in chunk ',
    paste(which(bad),collapse=','),': ',paste(as.character(parts[[which(bad)[1L]]]),collapse=' '))
  out<-vector('list',length(indices))
  for(i in seq_along(chunks)) out[match(chunks[[i]],indices)]<-parts[[i]]$values
  seconds<-vapply(parts,`[[`,numeric(1),'seconds')
  attr(out,'performance')<-data.table(scheduler=scheduler,workers=cores,chunks=length(chunks),
    loci=length(indices),wall_seconds=proc.time()[['elapsed']]-start,
    worker_seconds=sum(seconds),chunk_p50_seconds=unname(quantile(seconds,.5)),
    chunk_max_seconds=max(seconds))
  out
}

pm_run <- function(cfg) {
  started<-Sys.time()
  message('[POP][BUILD] ',PM_BUILD)
  pm_audit_selfcheck()
  if(isTRUE(cfg$rescue_require_callable)) pm_regional_gate_preflight()
  cfg$manifest<-pm_path(cfg$manifest)
  cfg$outdir<-pm_path(cfg$outdir,must=FALSE)
  dir.create(cfg$outdir,recursive=TRUE,showWarnings=FALSE)
  final_path<-file.path(cfg$outdir,'Atlas_merge.rds')
  if(file.exists(final_path) && !isTRUE(cfg$overwrite))
    stop('Refusing to overwrite existing Atlas_merge.rds: ',final_path,
      '. Use a NEW --outdir or explicitly --overwrite true.')
  checkpoint_path<-file.path(cfg$outdir,'00_stage_checkpoints.tsv')
  if(file.exists(checkpoint_path)) unlink(checkpoint_path)
  performance_path<-file.path(cfg$outdir,'00_evidence_performance_PRIVATE.tsv')
  if(file.exists(performance_path)) unlink(performance_path)
  fwrite(data.table(sample_id=character(),batch=integer(),scheduler=character(),workers=integer(),
    chunks=integer(),loci=integer(),wall_seconds=double(),worker_seconds=double(),
    chunk_p50_seconds=double(),chunk_max_seconds=double()),performance_path,sep='\t')
  pm_checkpoint(cfg$outdir,'START',paste0('catalog=',cfg$catalog_id),start=started)
  data.table::setDTthreads(1L)
  if(nzchar(cfg$softclip_config)) cfg$softclip_options<-jsonlite::fromJSON(cfg$softclip_config)
  e<-pm_modules(cfg)
  if(isTRUE(cfg$clip_enabled) && isTRUE(cfg$source_clip_cache_enabled))
    pm_clip_cache_preflight(cfg,e)
  m<-pm_manifest(cfg$manifest,e,cfg)
  # Only one sample contributes PRIMARY per locus; sex=U excludes X/Y from AF.
  saveRDS(m,file.path(cfg$outdir,'00_manifest_private.rds'))
  pm_checkpoint(cfg$outdir,'MANIFEST',paste0('samples=',nrow(m)),start=started)
  loaded<-vector('list',nrow(m))
  for(i in seq_len(nrow(m))) {
    loaded[[i]]<-pm_load_sample(m[i],e,cfg)
    message('[POP][LOAD_RDS] ',m$sample_id[[i]],' rows=',nrow(loaded[[i]]$loci))
  }
  all_loci<-rbindlist(lapply(loaded,`[[`,'loci'),fill=TRUE)
  all_frags<-rbindlist(lapply(loaded,`[[`,'fragments'),fill=TRUE)
  p0<-all_loci[primary_eligible==TRUE]
  sex_excluded<-all_loci[is.na(ploidy)|ploidy==0L,
    .(sample_id,uid,CHROM,POS,SVLEN,sex_region,ploidy,GT)]
  rm(loaded)
  message('[POP][PRIMARY] eligible=',nrow(p0),' | excluded sex/unmodeled=',nrow(sex_excluded))
  if(!nrow(p0)) return(pm_empty_run(m,all_loci,all_frags,sex_excluded,e,cfg,started))
  primary_edges<-pm_primary_edges(p0,e,cfg)
  clusters<-pm_cluster(p0,primary_edges,cfg)
  loci<-clusters$loci;p<-clusters$primary
  if(anyDuplicated(p$uid) || anyDuplicated(p[,.(locus_id,sample_id)]))
    stop('PRIMARY identity invariant violated')
  pm_checkpoint(cfg$outdir,'PRIMARY_CLUSTER',sprintf('%d loci, %d observations, %d accepted edges',
      nrow(loci),nrow(p),sum(clusters$edge_decisions$accepted)),start=started)
  graph<-pm_links(loci,p,e,cfg)
  pm_checkpoint(cfg$outdir,'NEARBY_LOCUS_GRAPH',paste0('candidate_links=',nrow(graph),
    '; plausible=',sum(graph$plausible)),start=started)
  # Fast lookup: all primary members per sample, one representative per locus.
  by_sample<-split(p,p$sample_id)
  reps<-p[, .SD[order(-fifelse(is.finite(RQ_SCORE),RQ_SCORE,0),-fifelse(is.finite(alt_support),alt_support,0),abs(POS-median(POS)),uid)][1L],by=locus_id]
  rep_idx<-setNames(seq_len(nrow(reps)),reps$locus_id)
  by_locus<-split(p,p$locus_id)
  # Build a single cache per TSV.gz BEFORE forking any evidence workers.
  # All samples are included regardless of regional callability because every
  # sample can contribute source clips for a different target sample.
  clip_caches<-list()
  if(cfg$clip_enabled && isTRUE(cfg$source_clip_cache_enabled)) {
    clip_windows<-pm_clip_cache_windows(reps,cfg)
    message('[POP][CLIP_CACHE_PLAN] samples=',nrow(m),
      ' | representatives=',nrow(reps),' | merged_windows=',nrow(clip_windows),
      ' | window_bp=',cfg$clip_window,
      ' | evidence_independent_of_callability=TRUE')
    cache_dir<-file.path(cfg$outdir,'00_softclip_cache_PRIVATE')
    if(isTRUE(cfg$source_clip_cache_save_rds))
      dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
    m[,RA_CACHE_INDEX:=seq_len(.N)]
    for(k in seq_len(nrow(m))) {
      cache<-pm_clip_cache_prepare(m[k],clip_windows,e,cfg,cache_dir)
      clip_caches[[m$sample_id[[k]]]]<-cache
    }
    m[,RA_CACHE_INDEX:=NULL]
    message('[POP][CLIP_CACHE_GLOBAL_READY] samples=',length(clip_caches),
      ' | no_repeated_source_gzip_scans=TRUE')
  }
  graph_by_locus<-lapply(loci$locus_id,function(lid) pm_neighbors(lid,graph))
  names(graph_by_locus)<-loci$locus_id
  fragments_by_sample<-if(nrow(all_frags)) split(all_frags,all_frags$sample_id) else list()
  fragments_by_sample<-lapply(fragments_by_sample,function(x) {
    setorder(x,CHROM,POS);split(x,x$CHROM)
  })
  # Each target sample is loaded once for QC/VCF/softclip, then freed.
  rows<-vector('list',nrow(m));vcf_audits<-vector('list',nrow(m))
  for(t in seq_len(nrow(m))) {
    target<-m[t];sid<-target$sample_id[[1L]]
    started_sample<-Sys.time()
    prim<-by_sample[[sid]]
    if(is.null(prim)) prim<-p[0L]
    prim_i<-setNames(seq_len(nrow(prim)),prim$locus_id)
    v<-pm_read_sample_vcf(target,e,cfg)
    vcf_audits[[t]]<-v$audit
    missing_vcf<-any(!v$audit$status %in% c('OK','EMPTY'))
    vcf_split<-if(nrow(v$table)) split(v$table,v$table$CHROM) else list()
    frags<-fragments_by_sample[[sid]]
    if(is.null(frags)) frags<-list()
    q<-e$ra_quality(target$quality_regions[[1L]],cfg$chromosome_map,
      loci[,.(uid=locus_id,CHROM,POS)],cfg,target$sex[[1L]])
    q[,frac_ge20:=NA_real_]  # Regional mean depth; no artificial >=20x fraction.
    q_i<-setNames(seq_len(nrow(q)),q$uid)
    # Quality is computed once from matched BED intervals and reused by all
    # batches. PRIMARY observations remain governed by the original D1 RDS QC;
    # log BED disagreement without silently redefining PRIMARY membership.
    q_nonprimary<-q[!uid %in% prim$locus_id]
    n_q_pass<-sum(q_nonprimary$callable %in% TRUE)
    n_q_fail<-nrow(q_nonprimary)-n_q_pass
    n_prim_bed_fail<-sum(q$uid %in% prim$locus_id & !(q$callable %in% TRUE))
    q_fail_status<-q_nonprimary[!(callable %in% TRUE),.(n=.N),by=quality_status][order(-n)]
    message('[POP][REGIONAL_QC_PREFLIGHT] sample=',sid,
      ' | rescue_require_callable=',cfg$rescue_require_callable,
      ' | min_depth=',cfg$callable_min_depth,
      ' | max_depth_median_ratio=',cfg$callable_max_depth_ratio,
      ' | nonprimary_callable=',n_q_pass,
      ' | nonprimary_not_callable=',n_q_fail,
      ' | PRIMARY_with_noncallable_BED=',n_prim_bed_fail,
      ' | noncallable_reasons=',paste(head(paste0(q_fail_status$quality_status,':',
        q_fail_status$n),6L),collapse=';'))
    rm(q_nonprimary,q_fail_status)
    clips<-NULL
    if(cfg$clip_enabled && isTRUE(cfg$source_clip_cache_enabled) &&
       !is.null(clip_caches[[sid]])) {
      clips<-clip_caches[[sid]]$target
      message('[POP][CLIP_CACHE_REUSE] target=',sid,
        ' | source_and_target=shared | no_second_gzip_scan=TRUE')
    } else if(cfg$clip_enabled && cfg$clip_input_mode=='preload') {
      path<-target$softclips[[1L]]
      if(!grepl('\\.(fa|fasta)(\\.gz)?$',path,ignore.case=TRUE)) {
        tab<-as.data.table(e$load_softclip_table(path))
        if(!all(c('CHR','POS1') %in% names(tab))) stop('Malformed softclip table for ',sid)
        tab[,`:=`(CHROM=e$ra_chr(CHR),POS=suppressWarnings(as.numeric(POS1)))]
        tab<-tab[is.finite(POS)]
        setorder(tab,CHROM,POS)
        clips<-split(tab,tab$CHROM)
        rm(tab)
      }
    }
    workers<-min(as.integer(cfg$clip_cores),as.integer(cfg$cores),nrow(loci))
    n_clips_preloaded<-if(!is.null(clips) && is.list(clips))
      sum(vapply(clips,nrow,integer(1))) else NA_integer_
    message('[POP][EVIDENCE_BEGIN] sample=',sid,' | loci=',nrow(loci),
      ' | workers=',workers,' | clip_input=',cfg$clip_input_mode,
      ' | clip_rows=',if(is.na(n_clips_preloaded)) 'REGIONAL' else n_clips_preloaded,
      ' | clip_error_policy=',cfg$clip_error_policy,
      ' | progress_batch_loci=',cfg$progress_batch_loci,
      ' | scheduler=',cfg$evidence_scheduler,' | worker_chunk_loci=',cfg$worker_chunk_loci)
    work_one <- function(j) {
      if(workers>1L) data.table::setDTthreads(1L)
      loc<-loci[j];lid<-loc$locus_id[[1L]]
      reg<-as.character(e$sa_region(loc$CHROM[[1L]],loc$POS[[1L]],cfg$reference_build))
      pl<-as.integer(e$sa_ploidy(reg,target$sex_karyotype[[1L]]))
      if(!cfg$include_sex_chromosomes && reg!='AUTOSOME') pl<-NA_integer_
      if(lid %in% names(prim_i)) {
        src<-prim[prim_i[[lid]]]
        gt<-e$ra_gt_info(src$GT[[1L]],pl)
        dose<-gt$dosage
        if(is.finite(dose) && dose==0L) {
          dose<-NA_integer_;gt$gt<-NA_character_
          gt$qc<-'PRIMARY_WITH_REFERENCE_GT_REQUIRES_REVIEW'
        }
        x<-pm_empty_state(lid,sid)
        x[,`:=`(state='PRIMARY',reason='DIRECT_PRIMARY_MEMBER',
          primary_uid=src$uid[[1L]],source_uid=src$uid[[1L]],source_sample=sid,
          GT=gt$gt,dosage=as.integer(dose),genotype_qc=gt$qc,
          ploidy=pl,sex_region=reg,callable=TRUE,callability_status='PRIMARY_PRESENT',
          evidence_level='PRIMARY')]
        return(x)
      } else {
        k<-q_i[[lid]]
        qr<-if(length(k) && is.finite(k)) q[k] else
          data.table(callable=FALSE,quality_status='QUALITY_LOOKUP_MISSING',
            frac_ge20=NA_real_,mean_depth=NA_real_)
        source<-as.list(reps[rep_idx[[lid]]])
        source_row<-m[sample_id==source$sample_id]
        return(pm_one_state(loc,source,source_row,target,prim,
          graph_by_locus[[lid]],qr,vcf_split,frags,clips,
          source_row$softclips[[1L]],e,cfg,missing_vcf,
          source_clip_cache=if(cfg$source_clip_cache_enabled &&
            !is.null(clip_caches[[source$sample_id]]))
              clip_caches[[source$sample_id]]$source else NULL))
      }
    }
    # Parallel batches give the parent a genuine progress heartbeat (workers are
    # forked within each batch). Never claim that an in-flight worker is done.
    one<-vector('list',nrow(loci))
    processed<-0L
    clip_checked<-0L
    clip_bad<-0L
    qc_gate_skipped<-0L
    batches<-split(seq_len(nrow(loci)),
      (seq_len(nrow(loci))-1L) %/% as.integer(cfg$progress_batch_loci))
    for(batch_i in seq_along(batches)) {
      indices<-batches[[batch_i]]
      batch_start<-Sys.time()
      result<-pm_evidence_parallel(indices,work_one,cores=workers,
        scheduler=cfg$evidence_scheduler,chunk_loci=cfg$worker_chunk_loci)
      perf<-attr(result,'performance');perf[,`:=`(sample_id=sid,batch=batch_i)]
      setcolorder(perf,c('sample_id','batch',setdiff(names(perf),c('sample_id','batch'))))
      fwrite(perf,performance_path,sep='\t',append=TRUE)
      one[indices]<-result
      zz<-rbindlist(result,fill=TRUE)
      processed<-processed+nrow(zz)
      clip_checked<-clip_checked+sum(!is.na(zz$clip_status) &
        !zz$clip_status %in% c('NOT_EVALUATED_BY_POLICY','NO_CLIP_EVALUATION',
          'NOT_EVALUATED_QC_GATE'))
      clip_bad<-clip_bad+sum(zz$clip_status %in% 'CLIP_ERROR',na.rm=TRUE)
      qc_gate_batch<-sum(zz$quality_gate_skipped %in% TRUE)
      qc_gate_skipped<-qc_gate_skipped+qc_gate_batch
      summary_state<-zz[,.(n=.N),by=state][order(-n)]
      state_text<-paste(paste0(summary_state$state,':',summary_state$n),collapse=';')
      message('[POP][EVIDENCE_PROGRESS] sample=',sid,
        ' | batch=',batch_i,'/',length(batches),
        ' | completed=',processed,'/',nrow(loci),
        ' | batch_sec=',round(as.numeric(difftime(Sys.time(),batch_start,units='secs')),1),
        ' | clip_evaluated_cumulative=',clip_checked,
        ' | clip_errors_cumulative=',clip_bad,
        ' | qc_gate_skipped_batch=',qc_gate_batch,
        ' | qc_gate_skipped_cumulative=',qc_gate_skipped,
        ' | batch_states=',state_text)
      if(clip_bad>0L) {
        examples<-unique(zz[clip_status=='CLIP_ERROR' & !is.na(clip_error_reason),
          clip_error_reason])
        message('[POP][SOFTCLIP_ERROR_SUMMARY] sample=',sid,
          ' | errors=',clip_bad,' | example_reasons=',
          paste(head(examples,3L),collapse=' | '))
        if(identical(cfg$clip_error_policy,'abort'))
          stop('Technical softclip errors in ',sid,'; see diagnostic messages above')
      }
      rm(result,zz)
    }
    rows[[t]]<-rbindlist(one,fill=TRUE)
    if(nrow(rows[[t]])!=nrow(loci)) stop('Missing locus states for ',sid)
    sample_clip_summary<-rows[[t]][!is.na(clip_status),
      .(n=.N),by=clip_status][order(-n)]
    message('[POP][SOFTCLIP_SAMPLE_SUMMARY] sample=',sid,' | ',
      paste(paste0(sample_clip_summary$clip_status,':',sample_clip_summary$n),collapse=';'))
    message('[POP][REGIONAL_QC_GATE_SUMMARY] sample=',sid,
      ' | skipped_before_VCF_fragment_softclip=',
      sum(rows[[t]]$quality_gate_skipped %in% TRUE),
      ' | rescue_in_noncallable=',sum(rows[[t]]$state %in% c('RESCUE_GT','RESCUE_ASSUMED') &
        !(rows[[t]]$callable %in% TRUE)))
    if(isTRUE(cfg$rescue_require_callable) &&
       any(rows[[t]]$state %in% c('RESCUE_GT','RESCUE_ASSUMED') &
         !(rows[[t]]$callable %in% TRUE)))
      stop('[POP][REGIONAL_QC_GATE_INVARIANT] Rescue bypassed regional quality gate: ',sid)
    rm(v,vcf_split,clips,one,q)
    gc(verbose=FALSE)
    pm_checkpoint(cfg$outdir,paste0('EVIDENCE_',sid),sprintf('loci=%d | seconds=%.1f | clip_errors=%d',
      nrow(loci),as.numeric(difftime(Sys.time(),started_sample,units='secs')),clip_bad),start=started)
  }
  states<-rbindlist(rows,fill=TRUE)
  n_clip_errors<-sum(states$clip_status %in% 'CLIP_ERROR',na.rm=TRUE)
  message('[POP][SOFTCLIP_GLOBAL_SUMMARY] errors=',n_clip_errors,
    ' | rescue_softclip=',sum(states$rescue_type %in% 'SOFTCLIP',na.rm=TRUE),
    ' | evaluated=',sum(!is.na(states$clip_status) &
      !states$clip_status %in% c('NOT_EVALUATED_BY_POLICY','NO_CLIP_EVALUATION',
        'NOT_EVALUATED_QC_GATE')),
    ' | skipped_by_regional_QC=',sum(states$quality_gate_skipped %in% TRUE))
  if(n_clip_errors>0L && identical(cfg$clip_error_policy,'abort'))
    stop('[POP][SOFTCLIP_GLOBAL_FATAL] Technical clip errors persisted; refusing AF output')
  if(nrow(states)!=nrow(m)*nrow(loci) || anyDuplicated(states[,.(locus_id,sample_id)]))
    stop('Each population locus must have exactly one state for every sample')
  rescue_qc<-pm_resolve_reuse(states,graph)
  states<-rescue_qc$states
  if(isTRUE(cfg$rescue_require_callable) &&
     any(states$state %in% c('RESCUE_GT','RESCUE_ASSUMED') &
       !(states$callable %in% TRUE)))
    stop('[POP][REGIONAL_QC_GATE_INVARIANT] Noncallable rescue before AF')
  pm_checkpoint(cfg$outdir,'RESCUE_COMPETITION',paste0('reused=',nrow(rescue_qc$reused_evidence),
    '; linked_competing=',nrow(rescue_qc$competing_rescues)),start=started)
  freq<-pm_frequency(states,loci,m,e,cfg)
  sm<-freq$summary
  if(any(sm$N_PRIMARY+sm$N_RESCUE_GT+sm$N_RESCUE_ASSUMED+
         sm$N_REVIEW+sm$N_NOT_CALLABLE+sm$N_REF_CALLABLE+sm$N_SEX_EXCLUDED!=nrow(m)))
    stop('Population states do not sum to N samples for every locus')
  # Public table is SUMMARY-ONLY by default; optional vectors are NOT anonymous.
  catalog<-pm_catalog(sm,p,m,freq$public_strata,cfg)
  # Preserve de-identified links (no Atlas source IDs, read IDs, or sample IDs).
  catalog$plausible_locus_links<-graph[,.(locus_a,locus_b,CHROM,distance_bp,
    size_sim,same_class,same_event,alt_qc,link_basis,plausible,unresolved_nearby)]
  internal<-list(format_version='RetroAtlas_Population_Internal_v1',
    metadata=list(catalog_id=cfg$catalog_id,build=PM_BUILD,
      reference_build=cfg$reference_build,created_at=as.character(Sys.time()),
      n_samples=nrow(m),n_loci=nrow(loci),softclip_error_count=n_clip_errors,
      softclip_validation=if(n_clip_errors==0L) 'NO_TECHNICAL_CLIP_ERRORS' else 'FAILED_TECHNICAL_CLIP_ERRORS',
      caveat='Rescue evaluations use a best directly reconstructed PRIMARY representative per locus; comparative ALT haplotypes can require targeted review'),
    parameters=cfg,manifest=m,all_intra_loci=all_loci,
    primary_members=p,primary_candidate_edges=primary_edges,
    primary_edge_assignments=clusters$edge_decisions,
    loci=sm,locus_sample_states=states,
    neighboring_loci=graph,fragment_evidence=all_frags,
    vcf_load_audit=rbindlist(vcf_audits,fill=TRUE),
    reused_evidence=rescue_qc$reused_evidence,
    competing_rescues=rescue_qc$competing_rescues,
    sex_excluded=sex_excluded,private_frequency_strata=freq$private_strata)
  pm_export(internal,catalog,cfg,started)
}

# Export the same schema for populated and empty cohorts.
pm_export <- function(internal,catalog,cfg,started) {
  m<-internal$manifest;sm<-internal$loci;p<-internal$primary_members
  states<-internal$locus_sample_states;graph<-internal$neighboring_loci
  clusters<-list(edge_decisions=internal$primary_edge_assignments)
  vcf_audits<-list(internal$vcf_load_audit)
  rescue_qc<-list(reused_evidence=internal$reused_evidence,competing_rescues=internal$competing_rescues)
  freq<-list(public_strata=catalog$strata)
  loci<-sm
  internal_file<-file.path(cfg$outdir,paste0('population_',cfg$catalog_id,'_internal.rds'))
  catalog_file<-file.path(cfg$outdir,'Atlas_merge.rds')
  saveRDS(internal,internal_file,compress='gzip')
  saveRDS(catalog,catalog_file,compress='gzip')
  if(cfg$write_tsv) {
    pm_tsv(sm,file.path(cfg$outdir,'01_population_loci.tsv'))
    pm_tsv(states,file.path(cfg$outdir,'02_locus_sample_states_PRIVATE.tsv'))
    pm_tsv(p[,.(locus_id,sample_id,uid,CHROM,POS,SVLEN,GT,alt_support,
       event_class,event,completeness,element_type,source_ids,haplotype_set)],
       file.path(cfg$outdir,'03_primary_members_PRIVATE.tsv'))
    pm_tsv(graph,file.path(cfg$outdir,'04_nearby_locus_relationships.tsv'))
    pm_tsv(clusters$edge_decisions,file.path(cfg$outdir,'05_primary_edge_audit.tsv'))
    pm_tsv(rbindlist(vcf_audits,fill=TRUE),file.path(cfg$outdir,'06_vcf_load_audit.tsv'))
    pm_tsv(rescue_qc$reused_evidence,file.path(cfg$outdir,'07_evidence_reuse_audit.tsv'))
    pm_tsv(rescue_qc$competing_rescues,file.path(cfg$outdir,'08_competing_rescue_audit.tsv'))
    pm_tsv(catalog$clusters,file.path(cfg$outdir,'Atlas_merge_summary.tsv'))
    pm_tsv(freq$public_strata,file.path(cfg$outdir,'09_public_strata_suppressed.tsv'))
    pm_tsv(data.table(parameter=names(cfg),value=vapply(cfg,toString,character(1))),
      file.path(cfg$outdir,'00_parameters_PRIVATE.tsv'))
  }
  pm_checkpoint(cfg$outdir,'DONE',sprintf('samples=%d loci=%d primary=%d rescue_GT=%d rescue_assumed=%d review=%d',
    nrow(m),nrow(loci),sum(sm$N_PRIMARY),sum(sm$N_RESCUE_GT),
    sum(sm$N_RESCUE_ASSUMED),sum(sm$N_REVIEW)),start=started)
  message('[POP][INTERNAL] ',internal_file)
  message('[POP][CATALOG] ',catalog_file)
  invisible(internal)
}

pm_empty_run <- function(m,all_loci,all_frags,sex_excluded,e,cfg,started) {
  p<-copy(all_loci[0L]);p[,locus_id:=character()]
  loci<-data.table(locus_id=character(),CHROM=character(),POS=double(),POS_MIN=double(),
    POS_MAX=double(),SVLEN=double(),N_PRIMARY=integer(),N_PRIMARY_SAMPLES=integer(),
    CLASS_MODE=character(),EVENT_MODE=character(),Q_MIN=integer())
  states<-pm_empty_state('','')[0L]
  freq<-pm_frequency(states,loci,m,e,cfg)
  catalog<-pm_catalog(freq$summary,p,m,freq$public_strata,cfg)
  graph<-pm_links(loci,p,e,cfg)
  catalog$plausible_locus_links<-graph[,.(locus_a,locus_b,CHROM,distance_bp,
    size_sim,same_class,same_event,alt_qc,link_basis,plausible,unresolved_nearby)]
  edges<-pm_primary_edges(p,e,cfg)
  assignments<-copy(edges);assignments[,`:=`(accepted=logical(),decision=character())]
  internal<-list(format_version='RetroAtlas_Population_Internal_v1',
    metadata=list(catalog_id=cfg$catalog_id,build=PM_BUILD,reference_build=cfg$reference_build,
      created_at=as.character(Sys.time()),n_samples=nrow(m),n_loci=0L,softclip_error_count=0L,
      softclip_validation='NOT_EVALUATED_NO_PRIMARY',status='NO_ELIGIBLE_PRIMARY'),
    parameters=cfg,manifest=m,all_intra_loci=all_loci,primary_members=p,
    primary_candidate_edges=edges,primary_edge_assignments=assignments,loci=freq$summary,
    locus_sample_states=states,neighboring_loci=graph,fragment_evidence=all_frags,
    vcf_load_audit=data.table(sample_id=character(),caller=character(),path=character(),status=character(),n_records=integer()),
    reused_evidence=data.table(),competing_rescues=data.table(),sex_excluded=sex_excluded,
    private_frequency_strata=freq$private_strata)
  pm_checkpoint(cfg$outdir,'NO_ELIGIBLE_PRIMARY','Empty catalogue; no evidence evaluations',start=started)
  pm_export(internal,catalog,cfg,started)
}
if(sys.nframe()==0L) pm_run(pm_parse())
