#!/usr/bin/env Rscript
# RetroAtlas Compare — standalone, directional 1:1, Q10..Q1, 2026-10-12.
# No pre-existing population/catalog required; NEVER source the population merge.
# Requires data.table, Biostrings, pwalign, IRanges and accompanying modules.
# BUILD: RETROCOMPARE_Q10_FINAL_PLAUSIBLE_DUP_GUARD_RESCUE_FIRST_20261019
suppressPackageStartupMessages(library(data.table))

options(datatable.print.nrows=10L)
RA_COMPARE_BUILD <- "RETROCOMPARE_Q10_FINAL_PLAUSIBLE_DUP_GUARD_RESCUE_FIRST_20261019"
# Resolve companions from THIS script's directory, not the SLURM working directory.
# Running: Rscript /absolute/path/RetroAtlas_compare.R ... works from anywhere.
compare_script_dir <- local({
  args <- commandArgs(trailingOnly=FALSE)
  f <- sub('^--file=', '', grep('^--file=', args, value=TRUE))
  if(length(f)) dirname(normalizePath(f[[1L]],mustWork=TRUE)) else getwd()
})

ra_script_path <- function(path, must_exist=TRUE, description='Input') {
  if(length(path)!=1L || is.na(path) || !nzchar(path))
    stop(description, ' path must be a nonempty string')
  # Relative inputs and module names refer to the script directory by default.
  # Absolute paths remain untouched (e.g. user's /pool_sas inputs).
  is_absolute <- startsWith(path, '/') || startsWith(path, '~')
  candidate <- if(is_absolute) path.expand(path) else file.path(compare_script_dir, path)
  if(must_exist && !file.exists(candidate))
    stop(description, ' not found: ', candidate, '\n',
         'Expected companion scripts and relative manifest next to RetroAtlas_compare.R (',
         compare_script_dir, ').')
  if(must_exist) normalizePath(candidate, mustWork=TRUE) else candidate
}

cfg_defaults <- function() list(
  manifest='', outdir='RetroAtlas_COMPARE', compare_id='COMPARE_A_B',
  a='', b='', reference_build='GRCh38', cores=2L, io_cores=2L, clip_cores=2L,
  primary_min_support=3, min_svlen=50,
  require_pass_cover=TRUE, require_pass_region=TRUE,
  require_resolved_class=TRUE, require_resolved_event=TRUE,
  require_resolved_completeness=TRUE,
  primary_window=500, primary_min_size_sim=0.80,
  primary_seq_mode='conflicts', primary_seq_identity=95,
  primary_seq_coverage=90, primary_seq_veto=TRUE,
  primary_seq_max_alignments=5000L, primary_seq_max_cells=80000000,
  q10_pos=5, q10_size=0.98, q9_pos=10, q9_size=0.95,
  q8_pos=25, q8_size=0.90, q7_pos=50, q7_size=0.90,
  q6_pos=100, q6_size=0.90, q5_pos=250, q5_size=0.90,
  q4_pos=500, q4_size=0.85, q3_pos=500, q3_size=0.80,
  rescue_window=1000, rescue_size_sim=0.80,
  # Independent protection against over-counting allelic/duplicate PRIMARY representations.
  # The 1 kb window finds competitors; proximity alone NEVER establishes identity.
  duplicate_guard_enabled=TRUE, duplicate_window=1000,
  duplicate_geometry_window=100, duplicate_geometry_size_sim=0.85,
  duplicate_align_enabled=TRUE, duplicate_align_min_size_sim=0.80,
  duplicate_seq_identity=95, duplicate_seq_coverage=90,
  duplicate_max_alignments_per_source=4L,
  undermerge_seq_identity=95, undermerge_seq_coverage=90,
  vcf_window=500, vcf_size_sim=0.80, vcf_min_support=1,
  vcf_require_pass=FALSE, vcf_allow_single_alt_caller=TRUE,
  clip_enabled=TRUE, clip_input_mode='preload', clip_window=5000,
  clip_core_window=1000, clip_min_identity=90,
  clip_min_coverage=90, clip_min_clips=3,
  clip_min_cluster=3, clip_min_union_bp=200,
  fragment_enabled=TRUE, fragment_window=1000,
  fragment_max_components=4, fragment_max_span=1000,
  fragment_size_sim=0.80,
  callable_min_depth=0, callable_max_depth_ratio=3,
  chromosome_map='', reference_helpers='../reference_helpers.R',
  require_all_vcf_for_differential=TRUE,
  include_sex_chromosomes=TRUE, write_tsv=TRUE,
  softclip_config='', softclip_module='retroatlas_softclip_rescue_v2.R',
  sequence_module='RetroAtlas_primary_ALT_pairwise.R',
  sex_module='RetroAtlas_sexaware_frequency.R'
)
ra_bool <- function(x) {
  z <- tolower(as.character(x))
  if(z %in% c('true','t','1','yes','y')) return(TRUE)
  if(z %in% c('false','f','0','no','n')) return(FALSE)
  stop('Invalid logical flag: ',x)
}
ra_parse_args <- function(defaults=cfg_defaults(),argv=commandArgs(trailingOnly=TRUE)) {
  if('--help' %in% argv || '-h' %in% argv) {
    cat('RetroAtlas_compare.R --manifest two_samples.tsv --outdir COMPARE_DIR [--compare-id ID] [--a ID --b ID] [--option value ...]\n')
    cat('Options (all configurable):\n')
    for(n in names(defaults)) cat(sprintf('  --%-40s %s\n',gsub('_','-',n),as.character(defaults[[n]])))
    quit(status=0)
  }
  if(length(argv) %% 2L) stop('CLI expects --flag value pairs; a value is missing')
  x <- defaults
  if(!length(argv)) stop('--manifest is required; use --help')
  for(i in seq(1L,length(argv),by=2L)) {
    key <- gsub('-','_',sub('^--','',argv[[i]]))
    if(!startsWith(argv[[i]],'--') || !key %in% names(x)) stop('Unknown option: ',argv[[i]])
    old <- defaults[[key]]; value <- argv[[i+1L]]
    x[[key]] <- if(is.logical(old)) ra_bool(value) else if(is.numeric(old)) {
      y <- suppressWarnings(as.numeric(value)); if(!is.finite(y)) stop('Invalid number for ',key)
      y
    } else as.character(value)
  }
  x$cores <- as.integer(x$cores); x$io_cores <- as.integer(x$io_cores); x$clip_cores <- as.integer(x$clip_cores)
  if(x$cores<1 || x$io_cores<1 || x$clip_cores<1) stop('cores must be positive')
  if(x$primary_window<=0 || x$primary_window>5000 || x$rescue_window<x$primary_window ||
     x$primary_min_size_sim<=0 || x$primary_min_size_sim>1 ||
     x$rescue_size_sim<=0 || x$rescue_size_sim>1) stop('Invalid search/size limits')
  if(!x$primary_seq_mode %in% c('off','conflicts','all')) stop('primary-seq-mode: off/conflicts/all')
  if(!x$clip_input_mode %in% c('preload','region')) stop('clip-input-mode: preload/region')
  if(x$vcf_min_support<0) stop('vcf-min-support must be >=0')
  if(x$duplicate_window<0 || x$duplicate_window>100000 ||
     x$duplicate_geometry_window<0 || x$duplicate_geometry_window>x$duplicate_window ||
     x$duplicate_geometry_size_sim<=0 || x$duplicate_geometry_size_sim>1 ||
     x$duplicate_align_min_size_sim<=0 || x$duplicate_align_min_size_sim>1 ||
     x$duplicate_seq_identity<0 || x$duplicate_seq_identity>100 ||
     x$duplicate_seq_coverage<0 || x$duplicate_seq_coverage>100 ||
     x$duplicate_max_alignments_per_source<0 || x$duplicate_max_alignments_per_source>1000)
    stop('Invalid duplicate guard thresholds')
  if(x$primary_seq_identity<0 || x$primary_seq_identity>100 ||
     x$primary_seq_coverage<0 || x$primary_seq_coverage>100) stop('Invalid sequence thresholds')
  if(!nzchar(x$manifest)) stop('--manifest is required')
  x
}

ra_exists <- function(path) !is.na(path) && nzchar(path) && toupper(path)!='NA' && file.exists(path)
ra_chr <- function(z) {
  a <- toupper(sub('^CHR','',as.character(z),ignore.case=TRUE))
  a[a=='23'] <- 'X'; a[a=='24'] <- 'Y'; a[a %in% c('25','MT')] <- 'M'
  paste0('chr',a)
}
ra_first <- function(dt,aliases,required=FALSE,default=NA) {
  nm <- aliases[aliases %in% names(dt)]
  if(length(nm)) return(dt[[nm[[1L]]]])
  if(required) stop('Missing RDS field: ',paste(aliases,collapse='/'))
  rep(default,nrow(dt))
}
ra_number <- function(x) suppressWarnings(as.numeric(x))
ra_size_sim <- function(a,b) {
  if(!is.finite(a) || !is.finite(b) || a<=0 || b<=0) return(NA_real_)
  min(a,b)/max(a,b)
}
ra_gt_info <- function(gt,ploidy=2L) {
  out <- list(gt=NA_character_,dosage=NA_integer_,qc='GT_MISSING')
  if(length(gt)!=1L || is.na(gt) || !nzchar(gt) || gt %in% c('.','./.','.|.')) return(out)
  z <- strsplit(sub(':.*$','',gsub('|','/',as.character(gt),fixed=TRUE)),'/',fixed=TRUE)[[1L]]
  if(any(z=='.') || any(!grepl('^[0-9]+$',z))) {out$qc<-'GT_UNRESOLVED';return(out)}
  alleles <- as.integer(z)
  if(!is.finite(ploidy) || ploidy<=0L) {out$qc<-'PLOIDY_UNMODELED';return(out)}
  if(ploidy==1L && length(alleles)==2L && alleles[[1L]]==alleles[[2L]]) {
    out$gt<-as.character(as.integer(alleles[[1L]]>0L))
    out$dosage<-as.integer(alleles[[1L]]>0L)
    out$qc<-'DIPLOID_HOMO_NORMALIZED_HAPLOID';return(out)
  }
  if(length(alleles)!=ploidy) {out$qc<-'GT_PLOIDY_CONFLICT';return(out)}
  out$gt<-paste(as.integer(alleles>0L),collapse='/')
  out$dosage<-sum(alleles>0L)
  out$qc<-'OK'
  out
}
ra_resolved <- function(x) {
  z <- tolower(trimws(as.character(x)))
  !is.na(z) & !z %in% c('','.','na','nan','none','null','unknown',
                         'unresolved','not_resolved','not resolved')
}
ra_iupac <- function(s) length(s)==1L && !is.na(s) && nzchar(s) && grepl('^[ACGTRYSWKMBDHVN]+$',toupper(s))
ra_parallel <- function(X,FUN,cores=1L) {
  if(length(X)==0L) return(list())
  cores<-min(as.integer(cores),length(X))
  if(cores<=1L || .Platform$OS.type=='windows') return(lapply(X,FUN))
  ans <- parallel::mclapply(X,FUN,mc.cores=cores,mc.preschedule=TRUE)
  bad <- vapply(ans,inherits,logical(1),'try-error')
  if(any(bad)) {
    idx <- which(bad)
    detail <- vapply(idx,function(j) paste0('task ',j,': ',as.character(ans[[j]])[[1L]]),character(1))
    stop('Worker failure: ',paste(detail,collapse=' | '))
  }
  ans
}
# Partition any number of candidates; cut(breaks=1) is invalid in R.
ra_candidate_chunks <- function(n, cores) {
  n <- as.integer(n); cores <- as.integer(cores)
  if(length(n)!=1L || length(cores)!=1L || is.na(n) || is.na(cores) || n<0L || cores<1L)
    stop('Invalid candidate chunk dimensions')
  if(n==0L) return(list())
  nworker <- min(n, cores)
  groups <- ((seq_len(n)-1L) %% nworker)+1L
  unname(split(seq_len(n), groups))
}

# Minimal regression preflight run by the primary script itself.
ra_compare_preflight <- function() {
  if(length(ra_candidate_chunks(0L,16L))!=0L ||
     length(ra_candidate_chunks(1L,16L))!=1L ||
     !identical(ra_candidate_chunks(1L,16L)[[1L]],1L) ||
     !identical(sort(unlist(ra_candidate_chunks(7L,3L),use.names=FALSE)),seq_len(7L)))
    stop('INTERNAL_SELF_CHECK_FAILED: primary candidate partition')
  # Fail before loading multi-gigabyte inputs if the PRIMARY edge/Q schema is
  # internally inconsistent (regression for the Q-vs-q production failure).
  cfg <- cfg_defaults()
  a <- list(CHROM='chr1',POS=1000,SVLEN=100,event_class='Alu',
            event='single',completeness='near_full',seq_valid=TRUE)
  b <- a
  edge <- ra_pair_q(a,b,cfg)
  if(is.null(edge) || !identical(edge$Q,10L))
    stop('INTERNAL_SELF_CHECK_FAILED: expected Q10 primary pair')
  dtA <- data.table(uid='MOCK_A',GT='0/1',ploidy=2L)
  dtB <- data.table(uid='MOCK_B',GT='0/1',ploidy=2L)
  e <- as.data.table(c(list(uid_A='MOCK_A',uid_B='MOCK_B'),edge))
  chk <- ra_make_primary(e,dtA,dtB,NULL)
  if(nrow(chk)!=1L || chk$Q[[1L]]!=10L ||
     !identical(chk$genotype_identity[[1L]],TRUE) ||
     chk$Q_LABEL[[1L]]!='Q10')
    stop('INTERNAL_SELF_CHECK_FAILED: primary Q and genotype export')
  invisible(TRUE)
}

ra_normalize_intra <- function(path,sid,cfg) {
  if(!ra_exists(path)) stop('Missing final_rds for ',sid,': ',path)
  x <- readRDS(path)
  bundle <- is.list(x) && 'loci' %in% names(x)
  d <- as.data.table(if(bundle) x$loci else x)
  d[,RESCUE_ONLY:=FALSE]
  if(bundle && is.data.frame(x$rescue_candidates) && nrow(x$rescue_candidates)) {
    candidates <- as.data.table(copy(x$rescue_candidates));candidates[,RESCUE_ONLY:=TRUE]
    d <- rbindlist(list(d,candidates),fill=TRUE,use.names=TRUE)
  }
  n <- nrow(d)
  g <- function(a,def=NA) ra_first(d,a,FALSE,def)
  dt <- data.table(sample_id=sid,rds_row=seq_len(n),
    rescue_only=d$RESCUE_ONLY,
    CHROM=ra_chr(ra_first(d,c('FINAL_CHR','CHROM','chrom','chr'),TRUE)),
    ORIGINAL_CHROM=as.character(g(c('ORIGINAL_CHROM','CHROM'))),
    POS=ra_number(ra_first(d,c('POS','POS_MEDIAN','pos'),TRUE)),
    SVLEN=abs(ra_number(ra_first(d,c('svlen','SVLEN','SVLEN_MEDIAN'),TRUE))),
    GT=as.character(g(c('GT','genotype'))),
    LR=as.character(g(c('LR','lr'))),
    event_class=as.character(g(c('event_class','EVENT_CLASS','class','mei_class'))),
    event=as.character(g(c('event','EVENT'))),
    completeness=as.character(g(c('completeness','COMPLETENESS'))),
    ALT_SEQUENCE=as.character(g(c('ALT','alt','inserted_sequence'))),
    cover_status=as.character(g(c('cover_status','COVER_STATUS'),'PASS')),
    region_status=as.character(g(c('region_status','REGION_STATUS'),'PASS')),
    fragment_mode=as.character(g(c('FRAGMENT_RECONSTRUCTION_MODE'),'')),
    fragment_caution=as.logical(g(c('FRAGMENT_RECONSTRUCTION_CAUTION'),FALSE)),
    population_primary_eligible=as.logical(g(c('POPULATION_PRIMARY_ELIGIBLE','PRIMARY_ELIGIBLE'),TRUE)),
    RQ_SCORE=ra_number(g(c('RECONSTRUCTION_QUALITY_SCORE','RQ_SCORE'),5)))
  dt[,source_ids:=as.character(g(c('SOURCE_IDS'),NA_character_))]
  # These fields preserve intramuestra event identity and phase: no re-merging.
  dt[,source_event_ids:=as.character(g(c('SOURCE_EVENT_IDS'),NA_character_))]
  dt[,haplotype_set:=as.character(g(c('HAPLOTYPE_SET'),NA_character_))]
  dt[,alt_support:=primary_alt_support(d)]
  dt[is.na(population_primary_eligible),population_primary_eligible:=TRUE]
  dt[,uid:=paste(sample_id,rds_row,sep='::')]
  dt[,is_synthetic:=grepl('HIGHCONF|DROPOUT|COMPOUND|SYNTHETIC',ifelse(is.na(fragment_mode),'',fragment_mode),ignore.case=TRUE)]
  dt[,seq_valid:=vapply(ALT_SEQUENCE,ra_iupac,logical(1)) &
      is.finite(SVLEN) & abs(nchar(ALT_SEQUENCE)-SVLEN)/pmax(SVLEN,1)<0.20]
  dt[,primary_eligible:=!rescue_only & population_primary_eligible %in% TRUE &
    !is_synthetic & !(fragment_caution %in% TRUE) & is.finite(POS) & is.finite(SVLEN) & SVLEN>=cfg$min_svlen &
    is.finite(alt_support) & alt_support>=cfg$primary_min_support &
    (!cfg$require_resolved_class | ra_resolved(event_class)) &
    (!cfg$require_resolved_event | ra_resolved(event)) &
    (!cfg$require_resolved_completeness | ra_resolved(completeness))]
  fragment <- if(bundle && 'fragment_evidence' %in% names(x)) as.data.table(x$fragment_evidence) else data.table()
  if(nrow(fragment)) {
    f <- function(a,def=NA) ra_first(fragment,a,FALSE,def)
    fragment<-data.table(sample_id=sid,
      CHROM=ra_chr(f(c('CHROM','chrom','chr'))),POS=ra_number(f(c('POS','pos'))),
      SVLEN=abs(ra_number(f(c('svlen','SVLEN','size')))),
      event_class=as.character(f(c('event_class','EVENT_CLASS','class'))),
      event=as.character(f(c('event','EVENT'))),
      fragment_group_id=as.character(f(c('FRAGMENT_GROUP_ID','.fragment_group_id'))),
      fragment_haplotype=as.character(f(c('FRAGMENT_HAPLOTYPE','.fragment_haplotype'))),
      structural_pass=as.logical(f(c('structural_pass'),NA)))
    fragment[,fragment_uid:=paste(sample_id,'FRAG',seq_len(.N),sep='::')]
  }
  list(loci=dt,fragments=fragment)
}

ra_read_vcf <- function(path,sid,caller,cfg) {
  audit<-data.table(sample_id=sid,caller=caller,path=path,status='OK',n_records=0L)
  if(!ra_exists(path)) {audit[,status:='MISSING'];return(list(table=data.table(),audit=audit))}
  # fread() supports plain VCF and gz paths when zlib support is available;
  # use a gzip connection if not. No shell pipelines or unquoted paths.
  con <- if(grepl('\\.gz$',path,ignore.case=TRUE)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con))
  lines<-readLines(con,warn=FALSE)
  head<-which(startsWith(lines,'#CHROM'))
  if(!length(head)) {audit[,status:='INVALID_HEADER'];return(list(table=data.table(),audit=audit))}
  h<-head[[length(head)]]
  records<-lines[seq.int(h,length(lines))]
  if(length(records)<=1L) {audit[,status:='EMPTY'];return(list(table=data.table(),audit=audit))}
  d<-tryCatch(fread(text=paste(records,collapse='\n'),sep='\t',header=TRUE,fill=TRUE,
      showProgress=FALSE,colClasses='character'),error=function(e) NULL)
  if(is.null(d) || !all(c('#CHROM','POS','REF','ALT','INFO') %in% names(d))) {
    audit[,status:='MALFORMED'];return(list(table=data.table(),audit=audit))
  }
  inf <- function(s,k) {
    z <- regmatches(s,regexpr(paste0('(^|;)',k,'=[^;]+'),s,perl=TRUE))
    if(!length(z) || is.na(z) || !nzchar(z)) return(NA_character_)
    sub(paste0('^.*',k,'='),'',z)
  }
  n<-nrow(d)
  svt<-vapply(d$INFO,inf,character(1),'SVTYPE')
  sl<-vapply(d$INFO,inf,character(1),'SVLEN')
  svlen<-vapply(strsplit(ifelse(is.na(sl),'',sl),',',fixed=TRUE),function(z) if(length(z)) ra_number(z[[1L]]) else NA_real_,numeric(1))
  explicit <- !grepl('^<',d$ALT) & nchar(d$ALT)>nchar(d$REF)
  svlen[!is.finite(svlen) & explicit] <- (nchar(d$ALT)-nchar(d$REF))[!is.finite(svlen) & explicit]
  samplecols<-setdiff(names(d),c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'))
  gt<-rep(NA_character_,n);support<-rep(NA_real_,n)
  if(length(samplecols) && 'FORMAT' %in% names(d)) {
    if(length(samplecols)!=1L) stop('Ambiguous multi-sample VCF: ',path)
    samp<-as.character(d[[samplecols[[1L]]]])
    fmt<-as.character(d$FORMAT)
    getfmt<-function(k) vapply(seq_len(n),function(i){
      a<-strsplit(fmt[[i]],':',fixed=TRUE)[[1L]]
      v<-strsplit(samp[[i]],':',fixed=TRUE)[[1L]]
      j<-match(k,a);if(is.na(j)||j>length(v)) return(NA_character_)
      v[[j]]
    },character(1))
    gt<-getfmt('GT')
    for(k in c('DV','RE','SUPPORT','AD')) {
      z<-getfmt(k)
      if(k=='AD') z<-vapply(strsplit(ifelse(is.na(z),'',z),',',fixed=TRUE),
        function(a) if(length(a)>1L) a[[2L]] else NA_character_,character(1))
      fill<-!is.finite(support)
      support[fill]<-ra_number(z[fill])
    }
  }
  gtmiss<-is.na(gt) | gt %in% c('.','./.','.|.','')
  if(any(gtmiss)) gt[gtmiss]<-vapply(d$INFO[gtmiss],inf,character(1),'GENOTYPE')
  for(k in c('SUPPORT','RE','DV')) {
    z<-ra_number(vapply(d$INFO,inf,character(1),k))
    support[!is.finite(support)]<-z[!is.finite(support)]
  }
  out<-data.table(sample_id=sid,caller=caller,CHROM=ra_chr(d[['#CHROM']]),
    POS=ra_number(d$POS),SVLEN=abs(svlen),GT=gt,support=support,
    FILTER=as.character(d$FILTER),SVTYPE=svt,
    ALT=as.character(d$ALT),REF=as.character(d$REF),
    VCF_ID=as.character(d$ID))
  out<-out[(SVTYPE=='INS' | (!is.na(SVLEN) & explicit)) &
    is.finite(POS) & is.finite(SVLEN) & SVLEN>0]
  if(cfg$vcf_require_pass) out<-out[FILTER %in% c('PASS','.',NA_character_)]
  out[,vcf_alt_seq:=ifelse(!grepl('^<',ALT) &
    nchar(ALT)>nchar(REF) & startsWith(ALT,REF),substring(ALT,nchar(REF)+1L),NA_character_)]
  audit[,n_records:=nrow(out)]
  list(table=out,audit=audit)
}

# Same regional-depth contract as Single. Never gate positive evidence.
ra_quality <- function(regpath, mapping, loci, cfg, sex='U') {
  cov <- read_coverage(regpath, read_chromosome_map(mapping),
                       cfg$callable_min_depth, cfg$callable_max_depth_ratio, sex)
  z <- annotate_coverage(loci, cov)
  data.table(uid=as.character(z$uid), callable=is.na(z$EXCLUSION_REASON),
      quality_status=ifelse(is.na(z$EXCLUSION_REASON),'CALLABLE',z$EXCLUSION_REASON),
      mean_depth=z$depth)
}

ra_pair_q <- function(a,b,cfg,seq_module_ready=TRUE) {
  d<-abs(a$POS-b$POS)
  s<-ra_size_sim(a$SVLEN,b$SVLEN)
  if(a$CHROM!=b$CHROM || !is.finite(d) || d>cfg$primary_window ||
     !is.finite(s) || s<cfg$primary_min_size_sim) return(NULL)
  cls<-!is.na(a$event_class) && !is.na(b$event_class) && a$event_class==b$event_class
  evt<-!is.na(a$event) && !is.na(b$event) && a$event==b$event
  comp<-!is.na(a$completeness) && !is.na(b$completeness) && a$completeness==b$completeness
  needs_seq<-cfg$primary_seq_mode=='all' || (cfg$primary_seq_mode=='conflicts' && !cls)
  seq_status<-'NOT_EVALUATED'; identity<-NA_real_;coverage<-NA_real_;orientation<-NA_character_
  if(needs_seq) {
    if(isTRUE(a$seq_valid) && isTRUE(b$seq_valid) && seq_module_ready) {
      z<-tryCatch(retroatlas_alt_pairwise(a$ALT_SEQUENCE,b$ALT_SEQUENCE,
         params=list(min_identity=cfg$primary_seq_identity,
          min_reciprocal_coverage=cfg$primary_seq_coverage,
          max_align_cells=cfg$primary_seq_max_cells)),error=function(e) NULL)
      if(!is.null(z)) {
        seq_status<-as.character(z$seq_qc[[1L]])
        identity<-as.numeric(z$identity_pct[[1L]])
        coverage<-as.numeric(z$reciprocal_coverage_min_pct[[1L]])
        orientation<-as.character(z$orientation[[1L]])
      } else seq_status<-'ALIGNMENT_ERROR'
    } else seq_status<-'ALT_UNAVAILABLE_OR_SYNTHETIC'
  }
  goodseq<-identical(seq_status,'SEQ_CONCORDANT')
  if(!cls && !goodseq) return(NULL)
  if(cls && isTRUE(cfg$primary_seq_veto) && seq_status=='SEQ_DISCORDANT') return(NULL)
  q <- if(!cls) {
    if(d<=cfg$q6_pos && s>=cfg$q6_size) 2L else 1L
  } else if(d<=cfg$q10_pos && s>=cfg$q10_size && evt && comp) 10L
  else if(d<=cfg$q9_pos && s>=cfg$q9_size && evt && comp) 9L
  else if(d<=cfg$q8_pos && s>=cfg$q8_size && evt) 8L
  else if(d<=cfg$q7_pos && s>=cfg$q7_size) 7L
  else if(d<=cfg$q6_pos && s>=cfg$q6_size) 6L
  else if(d<=cfg$q5_pos && s>=cfg$q5_size) 5L
  else if(d<=cfg$q4_pos && s>=cfg$q4_size) 4L
  else if(d<=cfg$q3_pos && s>=cfg$q3_size) 3L
  else return(NULL)
  if(!cls && (d>cfg$primary_window || s<cfg$primary_min_size_sim)) return(NULL)
  list(Q=as.integer(q),reason=if(!cls)'ALT_CONFIRMS_ANNOTATION_CONFLICT'
    else if(evt && comp)'GEOMETRY_ANNOTATION_CONCORDANT'
    else 'GEOMETRY_WITH_EVENT_OR_COMPLETENESS_DISCORDANCE',
    dpos=d,size_sim=s,class_match=cls,event_match=evt,completeness_match=comp,
    seq_status=seq_status,seq_identity=identity,seq_coverage=coverage,
    seq_orientation=orientation, score=as.numeric(q)*1e6+
      1000*s-0.1*d+5*as.numeric(evt)+2*as.numeric(comp))
}

# Exact bipartite assignment via min-cost augmenting Hungarian algorithm.
# Matched edges have positive weights. Invalid pairs and unmatched dummies cost 0.
# Operates only on connected components of the sparse candidate graph.
ra_hungarian <- function(w) {
  n<-nrow(w);m<-ncol(w)
  if(!n || !m) return(integer())
  if(n>m) stop('Hungarian requires nrow<=ncol')
  cost<-max(w)-w
  u<-numeric(n+1L);v<-numeric(m+1L);p<-integer(m+1L);way<-integer(m+1L)
  # dummy index 1 represents 0 in the conventional algorithm.
  for(i in seq_len(n)) {
    p[[1L]]<-i;j0<-1L; minv<-rep(Inf,m+1L);used<-rep(FALSE,m+1L)
    repeat {
      used[[j0]]<-TRUE;i0<-p[[j0]]; delta<-Inf;j1<-1L
      for(j in 2:(m+1L)) if(!used[[j]]) {
        cur<-cost[i0,j-1L]-u[i0+1L]-v[j]
        if(cur<minv[[j]]) {minv[[j]]<-cur;way[[j]]<-j0}
        if(minv[[j]]<delta) {delta<-minv[[j]];j1<-j}
      }
      for(j in seq_len(m+1L)) {
        if(used[[j]]) {u[p[[j]]+1L]<-u[p[[j]]+1L]+delta;v[[j]]<-v[[j]]-delta}
        else minv[[j]]<-minv[[j]]-delta
      }
      j0<-j1
      if(p[[j0]]==0L) break
    }
    repeat {
      j1<-way[[j0]];p[[j0]]<-p[[j1]];j0<-j1
      if(j0==1L) break
    }
  }
  ans<-integer(n)
  for(j in 2:(m+1L)) if(p[[j]]>0L) ans[[p[[j]]]]<-j-1L
  ans
}

# Build maximum-weight bipartite components without transitive identity merge.
ra_assign <- function(edges) {
  if(!nrow(edges)) return(edges[0L])
  a<-sort(unique(edges$uid_A));b<-sort(unique(edges$uid_B))
  aid<-match(edges$uid_A,a);bid<-match(edges$uid_B,b)
  # Union-find over both sides; components can never contain >1 of either side
  # in a final matching but can contain alternative candidate edges.
  nn<-length(a)+length(b);parent<-seq_len(nn)
  root<-function(i) {while(parent[[i]]!=i)i<-parent[[i]];i}
  for(i in seq_len(nrow(edges))) {
    ra<-root(aid[[i]]);rb<-root(length(a)+bid[[i]])
    if(ra!=rb) parent[[rb]]<-ra
  }
  groups<-vapply(aid,root,integer(1))
  res<-list()
  for(ind in split(seq_len(nrow(edges)),groups)) {
    e<-edges[ind]
    ax<-sort(unique(e$uid_A));bx<-sort(unique(e$uid_B))
    k<-length(ax); nb<-length(bx)
    w<-matrix(0,nrow=k,ncol=nb+k)
    # A missing edge must never appear as an eligible match.
    for(ii in seq_len(nrow(e))) {
      w[match(e$uid_A[[ii]],ax),match(e$uid_B[[ii]],bx)]<-e$score[[ii]]
    }
    chosen<-ra_hungarian(w)
    for(i in seq_len(k)) {
      j<-chosen[[i]]
      if(j>0L && j<=nb) {
        hit<-which(e$uid_A==ax[[i]] & e$uid_B==bx[[j]])
        if(length(hit)) res[[length(res)+1L]]<-e[hit[[1L]]]
      }
    }
  }
  if(length(res)) rbindlist(res,fill=TRUE) else edges[0L]
}

ra_gt_equal <- function(gt_a,gt_b,pl_a,pl_b) {
  ia<-ra_gt_info(gt_a,pl_a);ib<-ra_gt_info(gt_b,pl_b)
  if(!is.finite(pl_a) || !is.finite(pl_b) || pl_a!=pl_b ||
     !is.finite(ia$dosage) || !is.finite(ib$dosage)) return(NA)
  identical(as.integer(ia$dosage),as.integer(ib$dosage))
}

# Guard against regression in BED parsing and same-window callability.
# Tiny empty-clip integration test: loads the real softclip module's public API
# and ensures an empty region is not considered a rescue or a positive genotype.
ra_softclip_regression_preflight <- function() {
  empty <- data.frame(CLIP_ID=character(),CHR=character(),POS0=numeric(),
     POS1=numeric(),SOFTCLIP_LEN=numeric(),SEQUENCE=character(),
     stringsAsFactors=FALSE)
  x <- retroatlas_clip_classify(chr='chr1',pos=1000,
       alt_seq=strrep('ACGT',50),query_table=empty,source_table=empty,
       detailed=TRUE)
  if(!identical(as.character(x$decision$resolution[[1L]]),'NO_CLIP_SUPPORT') ||
     !identical(as.integer(x$decision$n_good_alt_clips[[1L]]),0L))
    stop('[INTERNAL_CLIP_SELF_CHECK] empty softclip input must have no support')
  invisible(TRUE)
}

# Test the *full* VCF rescue evaluator (not just primary Q schema) before
# opening large input files. Regression for data.table z1[alt_positive]:
# the explicit logical predicate is required for unambiguous filtering.
ra_vcf_regression_preflight <- function() {
  cfg <- cfg_defaults()
  src <- list(CHROM='chr1', POS=1000, SVLEN=100, seq_valid=FALSE)
  mk <- function(callers, gt, support) data.table::data.table(
    CHROM=rep('chr1',length(callers)), POS=rep(1000,length(callers)),
    SVLEN=rep(100,length(callers)), caller=callers, GT=gt,
    support=as.numeric(support), vcf_alt_seq=rep(NA_character_,length(callers)))
  one <- ra_vcf_evaluate(src,mk('cutesv','0/1',4),cfg,2L)
  if(!identical(one$status,'VCF_ALT_GT_RESCUE') ||
     !identical(one$dosage,1L) || !identical(one$n_alt_callers,1L))
    stop('[INTERNAL_VCF_SELF_CHECK] 0/1 single-caller rescue incorrect')
  two <- ra_vcf_evaluate(src,mk(c('cutesv','severus','sniffles2'),
         c('0/1','1/1','1/1'),c(4,4,5)),cfg,2L)
  if(!identical(two$status,'VCF_ALT_GT_RESCUE') ||
     !identical(two$dosage,2L) || !isTRUE(two$gt_conflict))
    stop('[INTERNAL_VCF_SELF_CHECK] multicaller majority and GT conflict incorrect')
  ref <- ra_vcf_evaluate(src,mk('cutesv','0/0',5),cfg,2L)
  if(!identical(ref$status,'VCF_0_0_WITH_ALT_SUPPORT'))
    stop('[INTERNAL_VCF_SELF_CHECK] 0/0 + ALT reads must remain REVIEW')
  single_policy <- cfg
  single_policy$vcf_allow_single_alt_caller <- FALSE
  deny <- ra_vcf_evaluate(src,mk('severus','0/1',5),single_policy,2L)
  if(!identical(deny$status,'VCF_SINGLE_CALLER_REVIEW'))
    stop('[INTERNAL_VCF_SELF_CHECK] single-caller review policy ignored')
  hap <- ra_vcf_evaluate(src,mk('sniffles2','1',5),cfg,1L)
  if(!identical(hap$dosage,1L) || !identical(hap$gt,'1'))
    stop('[INTERNAL_VCF_SELF_CHECK] haploid ALT normalization incorrect')
  invisible(TRUE)
}

# Assign one category to each unmatched PRIMARY; evidence inspection is always
# directional: ALT / source clips from source, VCF / target clips from target.
ra_vcf_evaluate <- function(source,target_vcf,cfg,ploidy) {
  empty <- list(status='NO_VCF_SUPPORT',reason='NO_GEOMETRIC_INS_RECORD',
    callers='',alt_callers='',gt=NA_character_,dosage=NA_integer_,
    gt_conflict=FALSE,n_alt_callers=0L,n_records=0L,
    best_size_sim=NA_real_,best_dpos=NA_real_,alt_seq_qc='NOT_EVALUATED',
    selected_record_ids='',
    supporting_records=data.table())
  if(!nrow(target_vcf)) return(empty)
  z<-target_vcf[CHROM==source$CHROM &
      abs(POS-source$POS)<=cfg$vcf_window &
      is.finite(SVLEN) & SVLEN>0]
  if(!nrow(z)) return(empty)
  z[,dpos:=abs(POS-source$POS)]
  z[,size_sim:=pmin(SVLEN,source$SVLEN)/pmax(SVLEN,source$SVLEN)]
  z<-z[size_sim>=cfg$vcf_size_sim]
  if(!nrow(z)) {empty$status<-'VCF_NEAR_SIZE_DISCORDANT';empty$reason<-'VCF_INS_NEAR_BUT_SIZE_INCOMPATIBLE';return(empty)}
  # Complete explicit VCF ALT (if present) must not contradict the primary ALT.
  z[,alt_seq_qc:='NOT_EVALUABLE']
  if(isTRUE(source$seq_valid)) {
    idx<-which(!is.na(z$vcf_alt_seq) & nchar(z$vcf_alt_seq)>=cfg$min_svlen)
    for(i in idx) {
      s<-z$vcf_alt_seq[[i]]
      if(!ra_iupac(s)) next
      q<-tryCatch(retroatlas_alt_pairwise(source$ALT_SEQUENCE,s,
        params=list(min_identity=cfg$primary_seq_identity,
          min_reciprocal_coverage=cfg$primary_seq_coverage,
          max_align_cells=cfg$primary_seq_max_cells)),error=function(e)NULL)
      z$alt_seq_qc[[i]]<-if(is.null(q)) 'ALIGN_ERROR' else as.character(q$seq_qc[[1L]])
    }
  }
  z<-z[!alt_seq_qc %in% c('SEQ_DISCORDANT','SEQ_PARTIAL')]
  if(!nrow(z)) {
    empty$status<-'VCF_ALT_CONTRADICTS_PRIMARY';empty$reason<-'SEQ_DISCORDANT_OR_PARTIAL';return(empty)
  }
  z[,gt_norm:=vapply(GT,function(g) ra_gt_info(g,ploidy)$gt,character(1))]
  z[,gt_dose:=vapply(GT,function(g) ra_gt_info(g,ploidy)$dosage,integer(1))]
  z[,alt_positive:=is.finite(gt_dose) & gt_dose>0L &
      is.finite(support) & support>=cfg$vcf_min_support]
  z[,ref_but_alt_support:=!is.na(gt_dose) & gt_dose==0L &
      is.finite(support) & support>=cfg$vcf_min_support]
  z[,caller_priority:=match(caller,c('sniffles2','severus','cutesv'))]
  z[is.na(caller_priority),caller_priority:=99L]
  setorder(z,caller_priority,dpos,-size_sim)
  # One best compatible record per caller: no double-counting the same caller.
  z1<-unique(z,by='caller')
  pos<-z1[alt_positive == TRUE]
  empty$n_records<-nrow(z1)
  empty$callers<-paste(sort(unique(z1$caller)),collapse=';')
  empty$best_size_sim<-max(z1$size_sim,na.rm=TRUE)
  empty$best_dpos<-min(z1$dpos,na.rm=TRUE)
  empty$alt_seq_qc<-paste(unique(z1$alt_seq_qc),collapse=';')
  empty$supporting_records<-z1
  if(nrow(pos)) {
    doses<-pos$gt_dose
    count<-table(doses)
    top<-as.integer(names(count)[count==max(count)])
    if(length(top)>1L) {
      pos_by_rank<-pos[order(caller_priority)]
      chosen<-pos_by_rank$gt_dose[match(TRUE,pos_by_rank$gt_dose %in% top)]
    } else chosen<-top[[1L]]
    empty$status<-'VCF_ALT_GT_RESCUE'
    empty$reason<-if(nrow(pos)>1L) 'ALT_GT_MULTICALLER' else 'ALT_GT_SINGLE_CALLER'
    empty$alt_callers<-paste(sort(unique(pos$caller)),collapse=';')
    empty$n_alt_callers<-nrow(pos)
    empty$gt<-if(ploidy==1L) '1' else if(chosen==2L) '1/1' else '0/1'
    empty$dosage<-as.integer(chosen)
    empty$selected_record_ids<-paste(paste(pos$caller,pos$CHROM,pos$POS,pos$SVLEN,sep=':'),collapse=';')
    empty$gt_conflict<-length(unique(z1$gt_dose[is.finite(z1$gt_dose)]))>1L
    if(!cfg$vcf_allow_single_alt_caller && nrow(pos)<2L) {
      empty$status<-'VCF_SINGLE_CALLER_REVIEW';empty$reason<-'SINGLE_CALLER_POLICY'
    }
    return(empty)
  }
  if(any(z1$ref_but_alt_support)) {
    empty$status<-'VCF_0_0_WITH_ALT_SUPPORT';empty$reason<-'REFERENCE_GT_POSITIVE_READ_EVIDENCE'
  } else if(any(is.finite(z1$support) & z1$support>=cfg$vcf_min_support)) {
    empty$status<-'VCF_UNTYPED_ALT_SUPPORT';empty$reason<-'ALT_READS_WITHOUT_GENOTYPE'
  } else {empty$status<-'VCF_NEAR_NO_ALT';empty$reason<-'VCF_GEOMETRY_WITHOUT_ALT_GT'}
  empty
}

ra_fragment_evaluate <- function(source,frags,cfg) {
  if(!cfg$fragment_enabled || !nrow(frags)) return(list(status='NO_FRAGMENT',ids='',sum_svlen=NA_real_))
  z<-frags[CHROM==source$CHROM & is.finite(POS) &
       abs(POS-source$POS)<=cfg$fragment_window &
       is.finite(SVLEN) & SVLEN>0]
  if(!nrow(z)) return(list(status='NO_FRAGMENT',ids='',sum_svlen=NA_real_))
  if(!'fragment_group_id' %in% names(z)) return(list(status='FRAGMENT_AMBIGUOUS',ids='',sum_svlen=NA_real_))
  groups<-z[!is.na(fragment_group_id) & nzchar(fragment_group_id),
     .(n=.N,sum_svlen=sum(SVLEN),span=max(POS)-min(POS),
       resolved_anchor=any(!is.na(event_class) & event_class==source$event_class),
       ids=paste(fragment_uid,collapse=';')),
     by=fragment_group_id]
  valid<-groups[n>=2L & n<=cfg$fragment_max_components &
       span<=cfg$fragment_max_span & resolved_anchor &
       pmin(sum_svlen,source$SVLEN)/pmax(sum_svlen,source$SVLEN)>=cfg$fragment_size_sim]
  if(nrow(valid)==1L) return(list(status='FRAGMENT_RESCUE_CANDIDATE',
       ids=valid$ids[[1L]],sum_svlen=valid$sum_svlen[[1L]]))
  if(nrow(valid)>1L || nrow(z)) return(list(status='FRAGMENT_AMBIGUOUS',
       ids=paste(z$fragment_uid,collapse=';'),sum_svlen=NA_real_))
  list(status='NO_FRAGMENT',ids='',sum_svlen=NA_real_)
}

ra_softclip_evaluate <- function(source,src_clips,tgt_clips,src_tbl,tgt_tbl,cfg,
                                 competitors=NULL) {
  empty<-list(resolution='NO_CLIP_EVALUATION',reason='',n_good_alt_clips=NA_integer_,
       details=NULL)
  if(!cfg$clip_enabled || !isTRUE(source$seq_valid)) {
    empty$reason<-'ALT_NOT_DIRECTLY_RECONSTRUCTED_OR_CLIP_DISABLED'
    return(empty)
  }
  tryCatch({
    params<-list(window=cfg$clip_window,core_window=cfg$clip_core_window,
       rescue_min_clips=cfg$clip_min_clips,
       rescue_min_cluster=cfg$clip_min_cluster,
       min_identity=cfg$clip_min_identity,
       min_clip_coverage=cfg$clip_min_coverage,
       min_alt_union_bp=cfg$clip_min_union_bp)
    params <- utils::modifyList(if(is.null(cfg$softclip_options)) list() else cfg$softclip_options, params)
    # Explicitly supply source_clips and query_clips in the correct direction.
    full<-retroatlas_clip_classify(chr=source$CHROM,pos=round(source$POS),
       alt_seq=source$ALT_SEQUENCE,query_clips=if(is.null(tgt_tbl)) tgt_clips else NULL,
       source_clips=if(is.null(src_tbl)) src_clips else NULL,
       query_table=tgt_tbl,source_table=src_tbl,params=params,
       competitor_alts=competitors,detailed=TRUE)
    list(resolution=as.character(full$decision$resolution[[1L]]),
       reason=as.character(full$decision$reason[[1L]]),
       n_good_alt_clips=as.integer(full$decision$n_good_alt_clips[[1L]]),
       details=full)
  },error=function(e) list(resolution='CLIP_ERROR',reason=conditionMessage(e),
         n_good_alt_clips=NA_integer_,details=NULL))
}

# =============================================================================
# POST-PRIMARY 1 kb LINK GUARD. Intramuestra haplotype representations remain
# distinct. The link is a caution/provenance graph, NEVER an automatic merge.
# =============================================================================
ra_source_tokens <- function(source_ids,source_event_ids=NA_character_) {
  z <- paste(c(as.character(source_ids),as.character(source_event_ids)),collapse=';')
  if(is.na(z) || !nzchar(z)) return(character())
  # Event labels are meaningful ONLY within the same sample.
  hit <- regmatches(z,gregexpr('Atlas_[0-9]+',z,perl=TRUE))[[1L]]
  unique(hit[nzchar(hit)])
}

ra_dup_link_audit <- function(unmatched,all_source_primary,cfg) {
  empty <- data.table(source_uid=character(),source_sample=character(),
    other_primary_uid=character(),other_primary_pair_id=character(),
    CHROM=character(),source_pos=numeric(),neighbor_pos=numeric(),dpos=numeric(),
    size_sim=numeric(),same_class=logical(),same_event=logical(),
    same_source_event=logical(),shared_source_events=character(),
    source_haplotypes=character(),neighbor_haplotypes=character(),
    alt_exact=logical(),alt_seq_qc=character(),link_basis=character(),
    plausible_duplicate=logical(),neighbor_already_matched=logical())
  if(!isTRUE(cfg$duplicate_guard_enabled) ||
     !nrow(unmatched) || !nrow(all_source_primary)) return(empty)
  candidates <- list()
  for(i in seq_len(nrow(unmatched))) {
    src <- as.list(unmatched[i])
    nearby <- all_source_primary[CHROM==src$CHROM & is.finite(POS) &
      abs(POS-src$POS)<=cfg$duplicate_window & uid!=src$uid]
    if(!nrow(nearby)) next
    setorder(nearby,POS,uid)
    sid1 <- ra_source_tokens(src$source_ids,src$source_event_ids)
    align_count <- 0L
    for(j in seq_len(nrow(nearby))) {
      tgt <- as.list(nearby[j])
      sim <- ra_size_sim(src$SVLEN,tgt$SVLEN)
      dpos <- abs(src$POS-tgt$POS)
      class_ok <- !is.na(src$event_class) && !is.na(tgt$event_class) &&
        tolower(src$event_class)==tolower(tgt$event_class)
      event_ok <- !is.na(src$event) && !is.na(tgt$event) &&
        tolower(src$event)==tolower(tgt$event)
      common <- intersect(sid1,ra_source_tokens(tgt$source_ids,tgt$source_event_ids))
      same_source <- length(common)>0L
      seq_exact <- isTRUE(src$seq_valid) && isTRUE(tgt$seq_valid) &&
        identical(toupper(src$ALT_SEQUENCE),toupper(tgt$ALT_SEQUENCE))
      seq_qc <- if(seq_exact) 'EXACT_ALT' else 'NOT_EVALUATED'
      geometric <- isTRUE(class_ok) && isTRUE(event_ok) && is.finite(sim) &&
        sim>=cfg$duplicate_geometry_size_sim && dpos<=cfg$duplicate_geometry_window
      # Longer-distance associations need a shared original event or direct ALT
      # concordance. Similar classes/sizes within 1 kb are NOT sufficient.
      if(!same_source && !seq_exact && !geometric &&
         isTRUE(cfg$duplicate_align_enabled) && isTRUE(class_ok) &&
         is.finite(sim) && sim>=cfg$duplicate_align_min_size_sim &&
         isTRUE(src$seq_valid) && isTRUE(tgt$seq_valid) &&
         align_count<cfg$duplicate_max_alignments_per_source) {
        align_count <- align_count+1L
        pw <- tryCatch(retroatlas_alt_pairwise(src$ALT_SEQUENCE,tgt$ALT_SEQUENCE,
           params=list(min_identity=cfg$duplicate_seq_identity,
             min_reciprocal_coverage=cfg$duplicate_seq_coverage,
             max_align_cells=cfg$primary_seq_max_cells)),
           error=function(e)NULL)
        seq_qc <- if(is.null(pw)) 'ALIGNMENT_UNAVAILABLE' else
          as.character(pw$seq_qc[[1L]])
      }
      seq_confirmed <- seq_qc %in% c('EXACT_ALT','SEQ_CONCORDANT')
      flagged <- same_source || seq_confirmed || geometric
      basis <- if(same_source) 'SHARED_INTRASAMPLE_SOURCE_EVENT'
        else if(seq_confirmed) 'ALT_SEQUENCE_CONCORDANT'
        else if(geometric) 'NEAR_GEOMETRY_CLASS_EVENT'
        else 'PROXIMITY_ONLY_NOT_PROVEN'
      pair <- if('matched_pair_id' %in% names(nearby)) tgt$matched_pair_id
        else NA_character_
      candidates[[length(candidates)+1L]] <- data.table(
        source_uid=as.character(src$uid),source_sample=as.character(src$sample_id),
        other_primary_uid=as.character(tgt$uid),other_primary_pair_id=as.character(pair),
        CHROM=as.character(src$CHROM),source_pos=as.numeric(src$POS),
        neighbor_pos=as.numeric(tgt$POS),dpos=as.numeric(dpos),
        size_sim=as.numeric(sim),same_class=isTRUE(class_ok),same_event=isTRUE(event_ok),
        same_source_event=isTRUE(same_source),shared_source_events=paste(common,collapse=';'),
        source_haplotypes=as.character(src$haplotype_set),
        neighbor_haplotypes=as.character(tgt$haplotype_set),
        alt_exact=isTRUE(seq_exact),alt_seq_qc=seq_qc,
        link_basis=basis,plausible_duplicate=flagged,
        neighbor_already_matched=!is.na(pair) && nzchar(pair))
    }
  }
  if(!length(candidates)) return(empty)
  rbindlist(candidates,fill=TRUE)
}

# Pure, auditable conservative veto. Does not turn a REVIEW into PRIMARY_MATCH.
ra_duplicate_guard <- function(links) {
  if(!nrow(links)) return(list(block=FALSE,reason='NO_INTERNAL_NEIGHBORS',
    subtype=NA_character_,ids='',pair_ids='',source_events='',bases=''))
  flagged <- links[plausible_duplicate==TRUE]
  if(!nrow(flagged)) return(list(block=FALSE,reason='NO_PROVEN_DUPLICATE_LINK',
    subtype=NA_character_,ids='',pair_ids='',source_events='',bases=''))
  matched <- flagged[neighbor_already_matched==TRUE]
  subtype <- if(nrow(matched)) 'POSSIBLE_MULTIALLELIC_LOCUS_ALREADY_MATCHED'
    else 'POSSIBLE_MULTIALLELIC_OR_REDUNDANT_SOURCE_PRIMARY'
  list(block=TRUE,reason='PLAUSIBLE_INTRASAMPLE_REPRESENTATION_LINK',
    subtype=subtype,
    ids=paste(sort(unique(flagged$other_primary_uid)),collapse=';'),
    pair_ids=paste(sort(unique(flagged$other_primary_pair_id[
      !is.na(flagged$other_primary_pair_id) & nzchar(flagged$other_primary_pair_id)])),collapse=';'),
    source_events=paste(sort(unique(unlist(strsplit(paste(flagged$shared_source_events,collapse=';'),
      ';',fixed=TRUE),use.names=FALSE))),collapse=';'),
    bases=paste(sort(unique(flagged$link_basis)),collapse=';'))
}

ra_duplicate_regression_preflight <- function() {
  cfg <- cfg_defaults()
  cfg$duplicate_align_enabled <- FALSE
  x <- data.table(uid='A::1',sample_id='A',CHROM='chr5',POS=1000,SVLEN=127,
    event_class='SVA',event='single',source_ids='Atlas_13481',
    source_event_ids='Atlas_13481',haplotype_set='1,1',
    ALT_SEQUENCE=strrep('A',127),seq_valid=TRUE)
  y <- copy(x)
  y[,`:=`(uid='A::2',POS=1084,SVLEN=175,source_ids='Atlas_13481',
    haplotype_set='1,0',ALT_SEQUENCE=strrep('C',175),
    matched_pair_id='PAIR::A::2::B::1')]
  r <- ra_dup_link_audit(x,y,cfg)
  if(nrow(r)!=1L || !isTRUE(r$plausible_duplicate[[1L]]) ||
     !isTRUE(ra_duplicate_guard(r)$block) ||
     ra_duplicate_guard(r)$subtype!='POSSIBLE_MULTIALLELIC_LOCUS_ALREADY_MATCHED')
    stop('[DUPLICATE_GUARD_CHECK] split Atlas event not identified')
  y[,`:=`(source_ids='Atlas_999',source_event_ids='Atlas_999',
    POS=1700,SVLEN=130,ALT_SEQUENCE=strrep('C',175))]
  r <- ra_dup_link_audit(x,y,cfg)
  if(nrow(r)!=1L || isTRUE(ra_duplicate_guard(r)$block))
    stop('[DUPLICATE_GUARD_CHECK] proximity must not imply identity')
  y[,`:=`(POS=1080,SVLEN=128,ALT_SEQUENCE=strrep('G',128))]
  r <- ra_dup_link_audit(x,y,cfg)
  if(!isTRUE(ra_duplicate_guard(r)$block) || r$link_basis[[1L]]!='NEAR_GEOMETRY_CLASS_EVENT')
    stop('[DUPLICATE_GUARD_CHECK] strong close geometry not identified')
  invisible(TRUE)
}

ra_compare_region <- function(chr,pos,sa, sb,cfg) {
  if(!isTRUE(cfg$include_sex_chromosomes) &&
       !sub('^chr','',chr) %in% as.character(1:22)) return('EXCLUDED_OPTION')
  reg<-sa_region(chr,pos,cfg$reference_build)
  if(reg=='UNSUPPORTED' || reg=='MT_UNMODELED') return('EXCLUDED_UNSUPPORTED')
  if(reg %in% c('Y_NONPAR','Y_PAR') &&
     (sa=='XX' || sb=='XX' || sa=='UNKNOWN' || sb=='UNKNOWN' || sa!=sb))
     return('EXCLUDED_Y_SEX_COMPARISON')
  if(reg=='Y_PAR') return('EXCLUDED_Y_PAR_DUPLICATE')
  if(reg %in% c('X_NONPAR','X_PAR') &&
      (sa=='UNKNOWN' || sb=='UNKNOWN' || !sa %in% c('XX','XY') || !sb %in% c('XX','XY')))
      return('EXCLUDED_X_UNKNOWN_KARYOTYPE')
  if(reg=='Y_NONPAR' && (sa!='XY' || sb!='XY')) return('EXCLUDED_Y_UNMODELED_KARYOTYPE')
  'INCLUDED'
}

ra_make_primary <- function(e,a_dt,b_dt,meta) {
  expected <- c('uid_A','uid_B','Q','score','reason','dpos','size_sim')
  missing <- setdiff(expected, names(e))
  if(length(missing)) stop('[PRIMARY_SCHEMA] Missing matched edge columns: ',
                           paste(missing,collapse=', '),
                           ' (available: ',paste(names(e),collapse=', '),')')
  if(anyNA(e$Q) || any(e$Q < 1L | e$Q > 10L))
    stop('[PRIMARY_SCHEMA] Q must be an integer from 1 to 10')
  z<-merge(e,a_dt,by.x='uid_A',by.y='uid',all.x=TRUE,sort=FALSE)
  # Preserve *all* traceability columns separately for both PRIMARY observations.
  setnames(z,setdiff(names(a_dt),'uid'),paste0('A_',setdiff(names(a_dt),'uid')))
  z<-merge(z,b_dt,by.x='uid_B',by.y='uid',all.x=TRUE,sort=FALSE)
  setnames(z,setdiff(names(b_dt),'uid'),paste0('B_',setdiff(names(b_dt),'uid')))
  z[,`:=`(category='PRIMARY_MATCH',Q_LABEL=paste0('Q',Q),
     genotype_identity=mapply(ra_gt_equal,A_GT,B_GT,A_ploidy,B_ploidy),
     genotype_identity_scope=ifelse(A_ploidy==B_ploidy,'COMPARABLE_PLOIDY','DIFFERENT_PLOIDY'))]
  ga<-Map(ra_gt_info,z$A_GT,z$A_ploidy)
  gb<-Map(ra_gt_info,z$B_GT,z$B_ploidy)
  z[,`:=`(GT_normalized_A=vapply(ga,`[[`,character(1),'gt'),
     GT_normalized_B=vapply(gb,`[[`,character(1),'gt'),
     GT_status_A=vapply(ga,`[[`,character(1),'qc'),
     GT_status_B=vapply(gb,`[[`,character(1),'qc'),
     GT_dosage_A=vapply(ga,`[[`,integer(1),'dosage'),
     GT_dosage_B=vapply(gb,`[[`,integer(1),'dosage'))]
  z[,primary_GT_ref_conflict:=(!is.na(GT_dosage_A) & GT_dosage_A==0L) |
        (!is.na(GT_dosage_B) & GT_dosage_B==0L)]
  z[primary_GT_ref_conflict==TRUE,genotype_identity:=NA]
  z[]
}

# A reconstructed secondary locus can confirm a callable source, never found one.
ra_reconstructed_evaluate <- function(src,partial,cfg) {
  good <- character()
  eligible <- partial[is.finite(alt_support) & alt_support>=cfg$vcf_min_support &
    seq_valid==TRUE & !is_synthetic & !(fragment_caution %in% TRUE) &
    ra_resolved(event_class) & ra_resolved(event) &
    pmin(SVLEN,src$SVLEN)/pmax(SVLEN,src$SVLEN)>=cfg$rescue_size_sim]
  eligible[,candidate_dosage:=vapply(GT,function(g) ra_gt_info(g,src$ploidy_target)$dosage,integer(1))]
  eligible <- eligible[is.na(candidate_dosage) | candidate_dosage>0L]
  if(isTRUE(src$seq_valid)) for(j in seq_len(nrow(eligible))) {
    q <- tryCatch(retroatlas_alt_pairwise(src$ALT_SEQUENCE,eligible$ALT_SEQUENCE[[j]],
      params=list(min_identity=cfg$primary_seq_identity,min_reciprocal_coverage=cfg$primary_seq_coverage,
                  max_align_cells=cfg$primary_seq_max_cells)),error=function(e)NULL)
    if(!is.null(q) && as.character(q$seq_qc[[1L]])=='SEQ_CONCORDANT') good <- c(good,eligible$uid[[j]])
  }
  unique(good)
}

ra_unmatched_one <- function(i,from,all_tgt,used_target_uids,vcf_tgt,
                             fragments_tgt,clip_src,clip_tgt,table_src,table_tgt,
                             q_target,cfg,vcf_missing_tgt,internal_links=data.table()) {
  src<-as.list(from[i])
  # PRIMARIES already assigned to other primary pairs are still competitors.
  near<-all_tgt[CHROM==src$CHROM &
      abs(POS-src$POS)<=cfg$rescue_window & is.finite(POS)]
  near<-near[order(abs(near$POS-src$POS))]
  near_primary<-near[primary_eligible==TRUE]
  partial<-near[primary_eligible!=TRUE]
  reconstructed_ids <- ra_reconstructed_evaluate(src,partial,cfg)
  links<-if(nrow(internal_links)) internal_links[source_uid==src$uid] else internal_links
  duplicate_guard<-ra_duplicate_guard(links)
  near_target_primary_ids<-paste(near_primary$uid,collapse=';')
  near_target_matched_pairs<-if('matched_pair_id' %in% names(near_primary))
    paste(unique(na.omit(near_primary$matched_pair_id[
      nzchar(near_primary$matched_pair_id)])),collapse=';') else ''
  review_subtype<-NA_character_
  cat<-'NOT_CALLABLE';reason<-'UNDETERMINED'
  linked_uid<-NA_character_;linked_dpos<-NA_real_
  rescue_type<-NA_character_;rescue_gt<-NA_character_;rescue_dosage<-NA_integer_
  vcf<-ra_vcf_evaluate(src,vcf_tgt,cfg,src$ploidy_target)
  fragment<-ra_fragment_evaluate(src,fragments_tgt,cfg)
  competitor_alts<-unique(near_primary[!is.na(ALT_SEQUENCE) & seq_valid==TRUE, ALT_SEQUENCE])
  competitor_alts<-setdiff(competitor_alts,src$ALT_SEQUENCE)
  if(length(competitor_alts)>6L) competitor_alts<-competitor_alts[seq_len(6L)]
  clip<-ra_softclip_evaluate(src,clip_src,clip_tgt,table_src,table_tgt,cfg,
                             competitors=if(length(competitor_alts)) competitor_alts else NULL)
  # Structural competition dominates rescue. An unmatched PRIMARY in the
  # strict window must never be silently re-labelled as a VCF/clip rescue.
  if(isTRUE(duplicate_guard$block)) {
    cat<-'REVIEW';reason<-duplicate_guard$reason
    review_subtype<-duplicate_guard$subtype
  } else if(nrow(near_primary[abs(POS-src$POS)<=cfg$primary_window])) {
    cat<-'REVIEW';reason<-'UNMATCHED_PRIMARY_WITHIN_PRIMARY_WINDOW'
    review_subtype<-'OPPOSITE_PRIMARY_WITHIN_500BP'
  } else if(nrow(near_primary) && any(near_primary$uid %in% used_target_uids)) {
    cat<-'REVIEW';reason<-'NEARBY_PRIMARY_ALREADY_ASSIGNED_IN_MATCH'
    review_subtype<-'OPPOSITE_PRIMARY_ALREADY_MATCHED_WITHIN_1KB'
  } else if(nrow(near_primary)) {
    # >500 and <=rescue_window: never PRIMARY. Sequence check is mandatory.
    compatible<-near_primary[SVLEN>0 &
        pmin(SVLEN,src$SVLEN)/pmax(SVLEN,src$SVLEN)>=cfg$rescue_size_sim]
    good<-list()
    if(isTRUE(src$seq_valid)) for(j in seq_len(nrow(compatible))) {
      t<-as.list(compatible[j])
      if(!isTRUE(t$seq_valid)) next
      q<-tryCatch(retroatlas_alt_pairwise(src$ALT_SEQUENCE,t$ALT_SEQUENCE,
        params=list(min_identity=cfg$undermerge_seq_identity,
           min_reciprocal_coverage=cfg$undermerge_seq_coverage,
           max_align_cells=cfg$primary_seq_max_cells)),error=function(e)NULL)
      if(!is.null(q) && as.character(q$seq_qc[[1L]])=='SEQ_CONCORDANT') {
        good[[length(good)+1L]]<-list(uid=t$uid,dist=abs(t$POS-src$POS))
      }
    }
    if(length(good)==1L && nrow(near_primary)==1L) {
      cat<-'RESCUE';reason<-'UNDERMERGE_SEQ_CONFIRMED_BUT_GT_UNRESOLVED'
      linked_uid<-good[[1L]]$uid;linked_dpos<-good[[1L]]$dist
      rescue_type<-'UNDERMERGE_PRIMARY_SEQUENCE_SUPPORTED'
      rescue_gt<-ra_gt_info(near_primary$GT[[1L]],src$ploidy_target)$gt
      rescue_dosage<-ra_gt_info(near_primary$GT[[1L]],src$ploidy_target)$dosage
    } else {
      cat<-'REVIEW';reason<-'NEARBY_PRIMARY_IDENTITY_OR_COMPETITION_UNRESOLVED'
      review_subtype<-'OPPOSITE_PRIMARY_501_1000BP_IDENTITY_UNRESOLVED'
    }
  } else if(length(reconstructed_ids)==1L) {
    cat<-'RESCUE';reason<-'RECONSTRUCTED_ALT_SEQUENCE_CONFIRMED'
    rescue_type<-'RECONSTRUCTED_CANDIDATE'
    candidate <- partial[uid==reconstructed_ids[[1L]]]
    rescue_gt<-ra_gt_info(candidate$GT[[1L]],src$ploidy_target)$gt
    rescue_dosage<-ra_gt_info(candidate$GT[[1L]],src$ploidy_target)$dosage
  } else if(vcf$status=='VCF_ALT_GT_RESCUE') {
    cat<-'RESCUE';reason<-vcf$reason;rescue_type<-'VCF_ALT_GT';
    rescue_gt<-vcf$gt;rescue_dosage<-vcf$dosage
  } else if(clip$resolution=='RESCUE_CANDIDATE') {
    cat<-'RESCUE';reason<-paste0('SOFTCLIP_',clip$reason)
    rescue_type<-'SOFTCLIP_CANDIDATE_GT_MISSING'
  } else if(fragment$status=='FRAGMENT_RESCUE_CANDIDATE') {
    cat<-'RESCUE';reason<-'FRAGMENT_SUM_GEOMETRY_ANCHORED'
    rescue_type<-'FRAGMENT_CANDIDATE_GT_MISSING'
  } else if(nrow(partial) || vcf$status %in% c('VCF_0_0_WITH_ALT_SUPPORT',
       'VCF_UNTYPED_ALT_SUPPORT','VCF_SINGLE_CALLER_REVIEW','VCF_ALT_CONTRADICTS_PRIMARY',
       'VCF_NEAR_SIZE_DISCORDANT','VCF_NEAR_NO_ALT') ||
       clip$resolution %in% c('REVIEW','CLIP_ERROR') ||
       fragment$status=='FRAGMENT_AMBIGUOUS') {
    cat<-'REVIEW';reason<-'LOCAL_STRUCTURAL_EVIDENCE_AMBIGUOUS'
    review_subtype<-'SECONDARY_EVIDENCE_AMBIGUOUS'
  } else if(!isTRUE(q_target$callable[[1L]])) {
    cat<-'NOT_CALLABLE';reason<-as.character(q_target$quality_status[[1L]])
  } else if(cfg$require_all_vcf_for_differential && isTRUE(vcf_missing_tgt)) {
    cat<-'NOT_CALLABLE';reason<-'VCF_CALLER_MISSING_DIFFERENTIAL_UNSAFE'
  } else if(clip$resolution=='NO_CLIP_EVALUATION' ||
            clip$resolution=='CLIP_ERROR') {
    cat<-'NOT_CALLABLE';reason<-'SOFTCLIP_NOT_EVALUABLE'
  } else {
    cat<-'DIFFERENTIAL';reason<-'TARGET_CALLABLE_NO_ALT_NO_AMBIGUOUS_EVIDENCE'
  }
  # Positive evidence survives low/unknown depth; structural ambiguity still wins.
  rescue_without_callability <- identical(cat,'RESCUE') && !isTRUE(q_target$callable[[1L]])
  suppressed_rescue <- if(isTRUE(duplicate_guard$block)) {
    if(length(reconstructed_ids)==1L) 'RECONSTRUCTED_CANDIDATE'
    else if(vcf$status=='VCF_ALT_GT_RESCUE') 'VCF_ALT_GT'
    else if(clip$resolution=='RESCUE_CANDIDATE') 'SOFTCLIP'
    else if(fragment$status=='FRAGMENT_RESCUE_CANDIDATE') 'FRAGMENT'
    else 'NONE'
  } else 'NOT_APPLICABLE'
  list(row=data.table(source_uid=src$uid,source_sample=src$sample_id,
       target_sample=if(src$sample_id==cfg$A_id) cfg$B_id else cfg$A_id,
       category=cat,reason=reason,CHROM=src$CHROM,POS=src$POS,
       SVLEN=src$SVLEN,ORIGINAL_CHROM=src$ORIGINAL_CHROM,event_class=src$event_class,event=src$event,
       completeness=src$completeness,GT_source=src$GT,
       source_ALT_SEQUENCE=src$ALT_SEQUENCE,source_support=src$alt_support,
       source_ploidy=src$ploidy, target_ploidy=src$ploidy_target,
       rescue_type=rescue_type,rescue_gt=rescue_gt,rescue_dosage=rescue_dosage,
       rescue_without_callability=rescue_without_callability,
       reconstructed_candidate_ids=paste(reconstructed_ids,collapse=';'),
       linked_primary_uid=linked_uid,linked_primary_dpos=linked_dpos,
       near_primary_count=nrow(near_primary),near_low_count=nrow(partial),
       review_subtype=review_subtype,
       same_sample_primary_near_count=nrow(links),
       plausible_internal_link_count=sum(links$plausible_duplicate %in% TRUE),
       linked_internal_primary_uids=duplicate_guard$ids,
       linked_existing_primary_pair_ids=duplicate_guard$pair_ids,
       shared_intrasample_source_events=duplicate_guard$source_events,
       internal_link_evidence=duplicate_guard$bases,
       rescue_blocked_by_internal_link=isTRUE(duplicate_guard$block),
       suppressed_rescue_channel=suppressed_rescue,
       nearby_opposite_primary_uids=near_target_primary_ids,
       nearby_opposite_matched_pair_ids=near_target_matched_pairs,
       duplicate_search_window_bp=cfg$duplicate_window,
       vcf_status=vcf$status,vcf_reason=vcf$reason,vcf_callers=vcf$callers,
       vcf_alt_callers=vcf$alt_callers,vcf_n_alt_callers=vcf$n_alt_callers,
       vcf_gt_conflict=vcf$gt_conflict,vcf_selected_record_ids=vcf$selected_record_ids,
       vcf_best_size_sim=vcf$best_size_sim,
       vcf_best_dpos=vcf$best_dpos,
       clip_resolution=clip$resolution,clip_reason=clip$reason,
       clip_n_good_alt=clip$n_good_alt_clips,
       clip_good_ids=if(!is.null(clip$details) && nrow(clip$details$clips))
         paste(unique(as.character(clip$details$clips$CLIP_ID[clip$details$clips$alt_match %in% TRUE])),collapse=';') else '',
       fragment_status=fragment$status,fragment_ids=fragment$ids,
       fragment_sum_svlen=fragment$sum_svlen,
       target_callable=isTRUE(q_target$callable[[1L]]),
       target_quality_status=as.character(q_target$quality_status[[1L]]),
       target_mean_depth=as.numeric(q_target$mean_depth[[1L]])),
       vcf_evidence=vcf$supporting_records,clip_audit=clip$details)
}

ra_output_tsv <- function(t,outdir) {
  if(!is.data.frame(t)) return(invisible(NULL))
  x<-copy(as.data.table(t))
  for(n in names(x)) if(is.list(x[[n]])) x[,(n):=vapply(get(n),function(v)
    paste(as.character(v),collapse=';'),character(1))]
  fwrite(x,file.path(outdir,paste0(toupper(deparse(substitute(t))),'.tsv')),sep='\t',na='NA')
}

ra_run <- function(cfg) {
  started<-Sys.time()
  message('[COMPARE][BUILD] ',RA_COMPARE_BUILD)
  ra_compare_preflight()
  ra_vcf_regression_preflight()
  message('[COMPARE][SELF_CHECK] candidate partition + PRIMARY Q/schema/GT + VCF 0/1,1/1,0/0 OK')
  # Both the manifest and output folder are resolved against this R script,
  # regardless of sbatch submission path or current shell working directory.
  cfg$manifest <- ra_script_path(cfg$manifest, description='Manifest')
  cfg$outdir <- ra_script_path(cfg$outdir, must_exist=FALSE, description='Output directory')
  dir.create(cfg$outdir,showWarnings=FALSE,recursive=TRUE)
  message('[COMPARE][SCRIPT_DIR] ',compare_script_dir)
  # Import the companion modules automatically. No module flags in SLURM are
  # needed; by default all three .R files live alongside RetroAtlas_compare.R.
  # Source into the functions' defining environment (not ra_run's frame).
  core_env<-environment(ra_pair_q)
  seq_path <- ra_script_path(cfg$sequence_module, description='ALT pairwise R module')
  sex_path <- ra_script_path(cfg$sex_module, description='Sex-aware R module')
  message('[COMPARE][AUTO_MODULE] ALT=',basename(seq_path))
  sys.source(seq_path,envir=core_env)
  message('[COMPARE][AUTO_MODULE] SEX=',basename(sex_path))
  sys.source(sex_path,envir=core_env)
  if(nzchar(cfg$softclip_config)) cfg$softclip_options <- jsonlite::fromJSON(cfg$softclip_config)
  if(cfg$clip_enabled) {
    # Softclip v2 defines its own ra_chr() (scalar): isolate its namespace to
    # prevent overwriting compare.R's vectorized ra_chr() function.
    clip_path <- ra_script_path(cfg$softclip_module, description='Softclip v2 R module')
    message('[COMPARE][AUTO_MODULE] SOFTCLIP=',basename(clip_path))
    clip_env<-new.env(parent=core_env)
    sys.source(clip_path,envir=clip_env)
    for(fun in c('retroatlas_clip_classify','load_softclip_table'))
      if(!exists(fun,envir=clip_env,inherits=FALSE,mode='function'))
        stop('Missing expected function ',fun,' in ',clip_path)
    assign('retroatlas_clip_classify',clip_env$retroatlas_clip_classify,envir=core_env)
    assign('load_softclip_table',clip_env$load_softclip_table,envir=core_env)
  }
  for(fun in c('retroatlas_alt_pairwise','sa_manifest','sa_region','sa_ploidy'))
    if(!exists(fun,envir=core_env,inherits=TRUE,mode='function'))
      stop('Missing expected module function: ',fun)
  sys.source(ra_script_path(cfg$reference_helpers), envir=core_env)
  data.table::setDTthreads(1L)  # Parallelism is controlled by ra_parallel, not nested worker pools.
  ra_duplicate_regression_preflight()
  message('[COMPARE][SELF_CHECK] regional depth contract + duplicate link guard OK')
  if(cfg$clip_enabled) {
    ra_softclip_regression_preflight()
    message('[COMPARE][SELF_CHECK] softclip module empty-region API OK')
  }
  m<-fread(cfg$manifest,na.strings=c('NA','N/A',''),colClasses=c(sample_id='character'))
  required<-c('sample_id','final_rds','vcf_cutesv','vcf_severus','vcf_sniffles2',
     'quality_regions','softclips','assembly','procedence')
  miss<-setdiff(required,names(m))
  if(length(miss)) stop('Missing manifest columns: ',paste(miss,collapse=', '))
  if(nrow(m)!=2L || anyDuplicated(m$sample_id)) stop('Compare requires exactly 2 distinct sample_id rows')
  if(!'sex' %in% names(m)) m[,sex:='U']
  if(!'sex_karyotype' %in% names(m)) m[,sex_karyotype:='UNKNOWN']
  if(!nzchar(cfg$a)) cfg$a<-m$sample_id[[1L]]
  if(!nzchar(cfg$b)) cfg$b<-m$sample_id[[2L]]
  if(cfg$a==cfg$b || !setequal(c(cfg$a,cfg$b),m$sample_id))
    stop('Invalid --a/--b: must be distinct IDs from manifest')
  m[,order_index:=match(sample_id,c(cfg$a,cfg$b))]
  setorder(m,order_index)
  m[,order_index:=NULL]
  m<-sa_manifest(m)
  kA<-m$sex_karyotype[[1L]];kB<-m$sex_karyotype[[2L]]
  if(!cfg$reference_build %in% c('GRCh38','GRCh38.p14','hg38','T2T-CHM13v2.0'))
      stop('Unsupported reference build')
  if(any(m$assembly != cfg$reference_build)) stop('Mixed/reference assembly mismatch')
  cfg$A_id<-cfg$a;cfg$B_id<-cfg$b
  for(i in 1:2) {
    for(col in c('final_rds','quality_regions'))
      if(!ra_exists(m[[col]][[i]])) stop('Missing required ',col,': ',m[[col]][[i]])
    if(cfg$clip_enabled && !ra_exists(m$softclips[[i]]))
      stop('Missing softclips: ',m$softclips[[i]])
  }
  message('[COMPARE][LOAD_RDS] ',cfg$a,' / ',cfg$b)
  raw<-ra_parallel(1:2,function(i) ra_normalize_intra(m$final_rds[[i]],m$sample_id[[i]],cfg),
     cores=min(cfg$io_cores,cfg$cores))
  a<-raw[[1L]]$loci;b<-raw[[2L]]$loci
  refresh_primary <- function(d,i) {
    quality <- ra_quality(m$quality_regions[[i]],cfg$chromosome_map,d,cfg,m$sex[[i]])
    d[,source_callable:=quality$callable[match(uid,quality$uid)]]
    d[,source_quality_status:=quality$quality_status[match(uid,quality$uid)]]
    d[,primary_eligible:=primary_eligible & source_callable %in% TRUE]
    d
  }
  a <- refresh_primary(a,1L); b <- refresh_primary(b,2L)
  fa<-raw[[1L]]$fragments;fb<-raw[[2L]]$fragments
  # Sex-filter all evidence before match/rescue, but preserve a QC audit.
  # Explicitly return both annotated tables. Do not modify transient
  # copies in a loop and assume the originals contain sex_filter.
  ra_tag_sex <- function(dt) {
    dt <- data.table::copy(dt)
    dt[, sex_region := sa_region(CHROM, POS, cfg$reference_build)]
    dt[, sex_filter := vapply(seq_len(.N), function(i)
      ra_compare_region(CHROM[[i]], POS[[i]], kA, kB, cfg), character(1))]
    dt[, ploidy := sa_ploidy(sex_region, ifelse(sample_id == cfg$a, kA, kB))]
    dt[, ploidy_target := sa_ploidy(sex_region, ifelse(sample_id == cfg$a, kB, kA))]
    dt
  }
  a <- ra_tag_sex(a)
  b <- ra_tag_sex(b)
  message('[COMPARE][SEX_FILTER] populated for A=',nrow(a),' and B=',nrow(b),
          ' | A columns=',paste(intersect(c('sex_filter','sex_region','ploidy'),names(a)),collapse=','),
          ' | B columns=',paste(intersect(c('sex_filter','sex_region','ploidy'),names(b)),collapse=','))
  if(!all(c('sex_filter','sex_region','ploidy','ploidy_target') %in% names(a)) ||
     !all(c('sex_filter','sex_region','ploidy','ploidy_target') %in% names(b)))
    stop('Sex-aware annotation was not attached to both sample tables')
  sex_excluded<-rbindlist(list(a[sex_filter!='INCLUDED'],b[sex_filter!='INCLUDED']),fill=TRUE)
  a<-a[sex_filter=='INCLUDED'];b<-b[sex_filter=='INCLUDED']
  # Input QC: event reconstructions tagged HIGHCONF/DROPOUT/COMPOUND/SYNTHETIC
  # cannot found a PRIMARY match. Keep for local fragment/review rescue.
  primary_a<-a[primary_eligible==TRUE]
  primary_b<-b[primary_eligible==TRUE]
  message('[COMPARE][PRIMARY] A=',nrow(primary_a),' B=',nrow(primary_b),
     ' | excluded sex=',nrow(sex_excluded))
  # Pair candidates by chromosome and coordinate intervals, avoiding A x B Cartesian product.
  candidates<-data.table(uid_A=character(),uid_B=character())
  if(nrow(primary_a) && nrow(primary_b)) {
    aidx<-primary_a[,.(uid,CHROM,start=as.numeric(POS),end=as.numeric(POS))]
    bidx<-primary_b[,.(uid,CHROM,start=as.numeric(POS-cfg$primary_window),
                         end=as.numeric(POS+cfg$primary_window))]
    setkey(bidx,CHROM,start,end)
    z<-foverlaps(aidx,bidx,which=TRUE,nomatch=NA)
    z<-z[!is.na(yid)]
    if(nrow(z)) candidates<-data.table(uid_A=aidx$uid[z$xid],uid_B=bidx$uid[z$yid])
  }
  setkey(primary_a,uid);setkey(primary_b,uid)
  if(nrow(candidates) && cfg$primary_seq_mode!='off') {
    aa<-primary_a[match(candidates$uid_A,primary_a$uid)]$event_class
    bb<-primary_b[match(candidates$uid_B,primary_b$uid)]$event_class
    needed<-if(cfg$primary_seq_mode=='all') nrow(candidates) else sum(is.na(aa) | is.na(bb) | aa!=bb)
    if(needed>cfg$primary_seq_max_alignments) stop('Primary sequence alignment guard exceeded: ',needed,
       ' > ',cfg$primary_seq_max_alignments,'. Increase --primary-seq-max-alignments or use conflicts mode')
  }
  # Evaluate quality by pairs. The optional alignment guard is global and
  # deterministic: candidates are ordered by coordinates and IDs first.
  if(nrow(candidates)) {
    candidates[,rA:=match(uid_A,primary_a$uid)]
    candidates[,rB:=match(uid_B,primary_b$uid)]
    setorder(candidates,rA,rB)
    # Safe for one candidate, one worker or many candidates.
    cand_chunks <- ra_candidate_chunks(nrow(candidates),cfg$cores)
    message('[COMPARE][CANDIDATES] ',nrow(candidates),
            ' | partitioned into ',length(cand_chunks),' chunks')
    # Pair-scoring is independent. Each pwalign call is bounded by max DP cells.
    evaluated<-ra_parallel(cand_chunks,function(ix) lapply(ix,function(i) {
      aa<-as.list(primary_a[candidates$rA[[i]]]);bb<-as.list(primary_b[candidates$rB[[i]]])
      pr<-ra_pair_q(aa,bb,cfg)
      if(is.null(pr)) return(NULL)
      as.data.table(c(list(uid_A=aa$uid,uid_B=bb$uid),pr))
    }),cores=min(cfg$cores,length(cand_chunks)))
    edgelist<-unlist(evaluated,recursive=FALSE)
    edgelist<-Filter(Negate(is.null),edgelist)
    edges<-if(length(edgelist)) rbindlist(edgelist,fill=TRUE) else data.table()
  } else {
    message('[COMPARE][CANDIDATES] ',nrow(candidates),' | no pairwise scoring needed')
    edges<-data.table()
  }
  if(nrow(edges)) {
    # Limit large connected components to prevent pathological O(n^3) memory/time.
    if(anyDuplicated(edges[,.(uid_A,uid_B)])) stop('Duplicate primary candidate edge')
    if(!'Q' %in% names(edges)) stop('[PRIMARY_SCHEMA] Pair scoring did not output Q; available fields: ',paste(names(edges),collapse=','))
    if(!is.numeric(edges$Q) || anyNA(edges$Q)) stop('[PRIMARY_SCHEMA] Invalid Q values')
    chosen<-ra_assign(edges)
  } else chosen<-data.table()
  if(nrow(chosen)) {
    primary_matches<-ra_make_primary(chosen,primary_a,primary_b,m)
    primary_matches[,`:=`(source='RDS_PRIMARY_BOTH',
         genotype_identity_reason=ifelse(is.na(genotype_identity),
           'GT_UNRESOLVED_DIFFERENT_PLOIDY_OR_PRIMARY_REF_CONFLICT',
           ifelse(genotype_identity,'DOSAGE_CONCORDANT','DOSAGE_DISCORDANT')))]
  } else primary_matches<-data.table()
  if(nrow(primary_matches) &&
     (anyDuplicated(primary_matches$uid_A) || anyDuplicated(primary_matches$uid_B)))
       stop('[PRIMARY_INVARIANT] More than one match per sample-side PRIMARY')
  matchedA<-if(nrow(chosen)) chosen$uid_A else character()
  matchedB<-if(nrow(chosen)) chosen$uid_B else character()
  # Link each original PRIMARY to the completed pair, without reusing or merging UIDs.
  primary_a[,matched_pair_id:=NA_character_]
  primary_b[,matched_pair_id:=NA_character_]
  a[,matched_pair_id:=NA_character_]
  b[,matched_pair_id:=NA_character_]
  if(nrow(primary_matches)) {
    primary_matches[,pair_id:=paste0(cfg$compare_id,'::',uid_A,'::',uid_B)]
    primary_a[primary_matches,on=.(uid=uid_A),matched_pair_id:=i.pair_id]
    primary_b[primary_matches,on=.(uid=uid_B),matched_pair_id:=i.pair_id]
    # Enrich the full target evidence tables as well: they are passed to
    # ra_unmatched_one(), while primary_a/b are used only for link construction.
    a[primary_matches,on=.(uid=uid_A),matched_pair_id:=i.pair_id]
    b[primary_matches,on=.(uid=uid_B),matched_pair_id:=i.pair_id]
  }
  unmatchedA<-primary_a[!uid %in% matchedA]
  unmatchedB<-primary_b[!uid %in% matchedB]
  if(nrow(primary_matches)) message('[COMPARE][Q_DISTRIBUTION] ',paste(names(table(primary_matches$Q)),as.integer(table(primary_matches$Q)),sep=':',collapse=' '))
  message('[COMPARE][MATCHED] ',nrow(primary_matches),
          ' | unmatched A=',nrow(unmatchedA),' B=',nrow(unmatchedB))
  # BEFORE expensive rescue, examine internal relationships of every unmatched
  # PRIMARY against all PRIMARY in the same sample (including matched ones).
  internal_a <- ra_dup_link_audit(unmatchedA,primary_a,cfg)
  internal_b <- ra_dup_link_audit(unmatchedB,primary_b,cfg)
  internal_links <- rbindlist(list(internal_a,internal_b),fill=TRUE)
  if(nrow(internal_links)) {
    nprotected <- uniqueN(internal_links[plausible_duplicate==TRUE]$source_uid)
  } else nprotected <- 0L
  message('[COMPARE][DUP_GUARD] window=',cfg$duplicate_window,
    ' | nearby pairs=',nrow(internal_links),
    ' | unmatched PRIMARY with plausible internal link=',nprotected)
  # No molecule is fused here: links are audit-only and can veto a new rescue.
  # VCF load: every caller for both samples; missing caller blocks clean DIFFERENTIAL.
  vcf_jobs<-ra_parallel(1:2,function(i) {
    vv<-lapply(c('cutesv','severus','sniffles2'),function(caller)
      ra_read_vcf(m[[paste0('vcf_',caller)]][[i]],m$sample_id[[i]],caller,cfg))
    list(table=rbindlist(lapply(vv,`[[`,'table'),fill=TRUE),
         audit=rbindlist(lapply(vv,`[[`,'audit'),fill=TRUE))
  },cores=min(cfg$io_cores,cfg$cores))
  va<-vcf_jobs[[1L]]$table;vb<-vcf_jobs[[2L]]$table
  vcf_audit<-rbindlist(lapply(vcf_jobs,`[[`,'audit'),fill=TRUE)
  vcf_missingA<-any(vcf_audit[sample_id==cfg$a]$status %in% c('MISSING','INVALID_HEADER','MALFORMED'))
  vcf_missingB<-any(vcf_audit[sample_id==cfg$b]$status %in% c('MISSING','INVALID_HEADER','MALFORMED'))
  # Quality at the other individual's PRIMARY positions, only for unmatched loci.
  qjobs<-ra_parallel(list(list(m=m[2L],targets=unmatchedA),
                          list(m=m[1L],targets=unmatchedB)),function(t) {
    ra_quality(t$m$quality_regions[[1L]],cfg$chromosome_map,t$targets,cfg,t$m$sex[[1L]])
  },cores=min(cfg$io_cores,cfg$cores))
  qa<-qjobs[[1L]];qb<-qjobs[[2L]]
  # Softclip cache in 1 pass per sample when preload; memory vs I/O selectable.
  ctA<-ctB<-NULL
  if(cfg$clip_enabled && cfg$clip_input_mode=='preload' &&
        (nrow(unmatchedA)||nrow(unmatchedB))) {
    message('[COMPARE][SOFTCLIP_PRELOAD] Reads may be large; use --clip-input-mode region for low RAM')
    caches<-ra_parallel(1:2,function(i) load_softclip_table(m$softclips[[i]]),
        cores=min(cfg$io_cores,cfg$cores))
    ctA<-caches[[1L]];ctB<-caches[[2L]]
  }
  # Run each unmatched primary in correct direction. Source clips are from the
  # individual's own file, query clips are from the opposite individual.
  inspect_direction<-function(src,from,to,target_vcf,fragment_target,q,
                               clips_source,clips_target,tab_source,tab_target,
                               missing_target,used_target_ids,source_internal_links) {
    if(!nrow(src)) return(list(rows=data.table(),audits=list(),vcfs=list()))
    job<-function(i) {
      qr<-q[uid==src$uid[[i]]]
      if(!nrow(qr)) qr<-data.table(callable=FALSE,quality_status='MISSING_QUALITY_LOOKUP',
                          mean_depth=NA_real_)
      ra_unmatched_one(i,src,to,used_target_ids,target_vcf,fragment_target,
           clips_source,clips_target,tab_source,tab_target,qr,cfg,missing_target,
           source_internal_links)
    }
    message('[COMPARE][RESCUE_DIRECTION] ',from$sample_id[[1L]],' -> ',to$sample_id[[1L]],
            ' | source PRIMARY=',nrow(src),' | workers=',min(cfg$clip_cores,cfg$cores,nrow(src)))
    res<-ra_parallel(seq_len(nrow(src)),job,cores=min(cfg$clip_cores,cfg$cores))
    list(rows=rbindlist(lapply(res,`[[`,'row'),fill=TRUE),
         audits=lapply(res,`[[`,'clip_audit'),
         vcfs=lapply(res,`[[`,'vcf_evidence'))
  }
  ar<-inspect_direction(unmatchedA,a,b,vb,fb,qa,m$softclips[[1L]],m$softclips[[2L]],
           ctA,ctB,vcf_missingB,matchedB,internal_a)
  br<-inspect_direction(unmatchedB,b,a,va,fa,qb,m$softclips[[2L]],m$softclips[[1L]],
           ctB,ctA,vcf_missingA,matchedA,internal_b)
  classified<-rbindlist(list(ar$rows,br$rows),fill=TRUE)
  if(!nrow(classified) && !'category' %in% names(classified))
    classified<-data.table(category=character(),rescue_type=character(),
       source_uid=character(),target_sample=character(),
       linked_primary_uid=character(),source_sample=character(),pair_id=character())
  # Reconcile reciprocal undermerges; A->B and B->A are ONE relationship,
  # not two independent rescues. Unreciprocated link is an identity REVIEW.
  if(nrow(classified)) {
    under<-which(classified$rescue_type=='UNDERMERGE_PRIMARY_SEQUENCE_SUPPORTED' &
          !is.na(classified$linked_primary_uid))
    remove<-integer()
    for(i in under) {
      if(i %in% remove) next
      j<-which(classified$source_uid==classified$linked_primary_uid[[i]] &
         classified$linked_primary_uid==classified$source_uid[[i]] &
         classified$rescue_type=='UNDERMERGE_PRIMARY_SEQUENCE_SUPPORTED')
      j<-setdiff(j,i)
      if(length(j)!=1L) {
        classified[i,`:=`(category='REVIEW',reason='UNDERMERGE_NONRECIPROCAL_OR_COMPETING',
                          review_subtype='OPPOSITE_UNDERMERGE_LINK_NONRECIPROCAL',
                          rescue_type='UNDERMERGE_AMBIGUOUS')]
      } else if(classified$source_sample[[i]]==cfg$a) {
        remove<-c(remove,j[[1L]])
      }
    }
    if(length(remove)) classified<-classified[-unique(remove)]
    classified[,pair_id:=paste0(cfg$compare_id,'::',source_uid)]
  }
  # A VCF record / clip ID may not rescue two different source loci in
  # the same query sample. Quarantine competing rescues instead of boosting AF.
  reuse_audit<-data.table()
  if(nrow(classified)) {
    supporting<-classified[category=='RESCUE' &
       rescue_type %in% c('VCF_ALT_GT','SOFTCLIP_CANDIDATE_GT_MISSING','RECONSTRUCTED_CANDIDATE')]
    if(nrow(supporting)) {
      rec<-rbindlist(lapply(seq_len(nrow(supporting)),function(i) {
        z<-supporting[i]
        raw<-if(z$rescue_type=='VCF_ALT_GT') z$vcf_selected_record_ids else if(z$rescue_type=='RECONSTRUCTED_CANDIDATE') z$reconstructed_candidate_ids else z$clip_good_ids
        keys<-strsplit(raw,';',fixed=TRUE)[[1L]]
        keys<-keys[nzchar(keys)]
        if(!length(keys)) return(data.table())
        data.table(source_uid=z$source_uid,target_sample=z$target_sample,
            evidence_id=keys,kind=z$rescue_type)
      }),fill=TRUE)
      if(nrow(rec)) {
        reuse_audit<-rec[,.(n_distinct_loci=uniqueN(source_uid),
            sources=paste(sort(unique(source_uid)),collapse=';')),
             by=.(target_sample,kind,evidence_id)][n_distinct_loci>1L]
        if(nrow(reuse_audit)) {
          reused<-unique(rec[reuse_audit,on=.(target_sample,kind,evidence_id)]$source_uid)
          classified[source_uid %in% reused &
             rescue_type %in% c('VCF_ALT_GT','SOFTCLIP_CANDIDATE_GT_MISSING','RECONSTRUCTED_CANDIDATE'),
             `:=`(category='REVIEW',reason='EVIDENCE_REUSED_ACROSS_SOURCE_LOCI',
                  review_subtype='EVIDENCE_REUSED_BY_MULTIPLE_SOURCE_LOCI')]
        }
      }
    }
  }
  # An internal link must NEVER be silently counted as RESCUE/DIFFERENTIAL.
  if(nrow(classified) && any(classified$rescue_blocked_by_internal_link==TRUE &
      classified$category!='REVIEW'))
    stop('[DUP_GUARD_INVARIANT] Linked source was not quarantined as REVIEW')
  rescue_protection_summary <- data.table(metric=c(
    'UNMATCHED_PRIMARY', 'NEARBY_INTERNAL_LINK_PAIRS',
    'UNMATCHED_WITH_PLAUSIBLE_INTERNAL_LINK',
    'RESCUES_BLOCKED_BY_INTERNAL_LINK',
    'VCF_RESCUES_BLOCKED', 'SOFTCLIP_RESCUES_BLOCKED',
    'FRAGMENT_RESCUES_BLOCKED'),
    value=c(nrow(unmatchedA)+nrow(unmatchedB),nrow(internal_links),nprotected,
      sum(classified$rescue_blocked_by_internal_link %in% TRUE),
      sum(classified$suppressed_rescue_channel=='VCF_ALT_GT',na.rm=TRUE),
      sum(classified$suppressed_rescue_channel=='SOFTCLIP',na.rm=TRUE),
      sum(classified$suppressed_rescue_channel=='FRAGMENT',na.rm=TRUE)))
  message('[COMPARE][DUP_GUARD_RESULT] REVIEW protected=',
    sum(classified$rescue_blocked_by_internal_link %in% TRUE),
    ' | suppressed ALT channels=',
    sum(classified$suppressed_rescue_channel %in% c('VCF_ALT_GT','SOFTCLIP','FRAGMENT'),na.rm=TRUE))
  if(nrow(classified) && any(classified$category=='DIFFERENTIAL' &
                             !(classified$target_callable %in% TRUE)))
    stop('[CALLABILITY_INVARIANT] DIFFERENTIAL needs a callable target')
  rescues<-classified[category=='RESCUE']
  review<-classified[category=='REVIEW']
  not_callable<-classified[category=='NOT_CALLABLE']
  differentials<-classified[category=='DIFFERENTIAL']
  if(nrow(primary_matches)) {
    primary_matches[,pair_id:=paste0(cfg$compare_id,'::',uid_A,'::',uid_B)]
  }
  # A single result table with exclusive categories, after preserving full
  # PRIMARY_MATCH detail in its own block.
  matches_flat<-if(nrow(primary_matches)) primary_matches[,.(pair_id,category,
      source_uid=uid_A,source_sample=cfg$a,target_sample=cfg$b,
      reason,Q,Q_LABEL,CHROM=A_CHROM,POS=A_POS,SVLEN=A_SVLEN)] else data.table()
  unmatched_flat<-if(nrow(classified)) classified[,.(pair_id,category,source_uid,
      source_sample,target_sample,reason,Q=NA_integer_,Q_LABEL=NA_character_,
      CHROM,POS,SVLEN,review_subtype,linked_internal_primary_uids,
      linked_existing_primary_pair_ids,rescue_blocked_by_internal_link,
      suppressed_rescue_channel,nearby_opposite_primary_uids,
      nearby_opposite_matched_pair_ids)] else data.table()
  loci_classified<-rbindlist(list(matches_flat,unmatched_flat),fill=TRUE)
  if(nrow(loci_classified) && (anyDuplicated(loci_classified$pair_id) ||
      any(!loci_classified$category %in% c('PRIMARY_MATCH','RESCUE','REVIEW','NOT_CALLABLE','DIFFERENTIAL'))))
      stop('Output classification invariant violated')
  if(!'category' %in% names(loci_classified)) loci_classified<-data.table(category=character())
  summary<-data.table(metric=c('PRIMARY_MATCH','RESCUE','REVIEW','NOT_CALLABLE','DIFFERENTIAL'),
       value=vapply(c('PRIMARY_MATCH','RESCUE','REVIEW','NOT_CALLABLE','DIFFERENTIAL'),
           function(k) sum(loci_classified$category==k),integer(1)))
  clip_audit_raw<-c(ar$audits,br$audits)
  clip_summaries<-Filter(Negate(is.null),lapply(seq_along(clip_audit_raw),function(i) {
    z<-clip_audit_raw[[i]]
    if(is.null(z)) return(NULL)
    d<-as.data.table(z$audit)
    d[,query_index:=i]
    d
  }))
  clip_decisions<-if(length(clip_summaries)) rbindlist(clip_summaries,fill=TRUE) else data.table()
  # Sample metadata accompanies exported categories, including empty tables.
  for(tbl in c('loci_classified','primary_matches','rescues','review','not_callable','differentials')) {
    d <- get(tbl)
    d[, `:=`(ASSEMBLY=rep(cfg$reference_build,.N), SEX_A=rep(m$sex[[1L]],.N),
      SEX_B=rep(m$sex[[2L]],.N), PROCEDENCE_A=rep(m$procedence[[1L]],.N),
      PROCEDENCE_B=rep(m$procedence[[2L]],.N))]
    assign(tbl,d)
  }
  # Exact inputs, settings and all audit evidence embedded in the comparison RDS.
  result<-list(
    format_version='RetroAtlas_compare_standalone_Q10_2026-10-14_rescue_first',
    metadata=list(compare_id=cfg$compare_id,sample_A=cfg$a,sample_B=cfg$b,
      karyotype_A=kA,karyotype_B=kB,reference_build=cfg$reference_build,
      sex_A=m$sex[[1L]],sex_B=m$sex[[2L]],
      procedence_A=m$procedence[[1L]],procedence_B=m$procedence[[2L]],
      chromosome_lengths=read_chromosome_map(cfg$chromosome_map)[,c('FINAL_CHR','expected_length')],
      rescue_policy='Positive rescue evidence precedes target callability',
      created_at=as.character(Sys.time()),
      genotype_identity_definition='normalized ALT dosage equal only at equal ploidy',
      no_population_input=TRUE),
    category_summary=summary, loci_classified=loci_classified,
    primary_matches=primary_matches,
    rescues=rescues,review=review,not_callable=not_callable,
    differentials=differentials,
    audit=list(primary_all_candidate_pairs=candidates,
       excluded_primary_loci=rbindlist(list(a[primary_eligible!=TRUE],b[primary_eligible!=TRUE]),fill=TRUE),
       primary_accepted_edges=edges,
       primary_chosen_edges=chosen,
       primary_unmatched_A=unmatchedA,
       primary_unmatched_B=unmatchedB,
       input_rds_audit=data.table(sample_id=c(cfg$a,cfg$b),
           n_rows=c(nrow(a),nrow(b)),n_primary=c(nrow(primary_a),nrow(primary_b))),
       vcf_load_audit=vcf_audit,
       vcf_evidence_by_unmatched=c(ar$vcfs,br$vcfs),
       clip_decisions=clip_decisions,
       rescue_evidence_reuse_audit=reuse_audit,
       plausible_internal_primary_links=internal_links,
       duplicate_guard_summary=rescue_protection_summary,
       clip_details=clip_audit_raw,
       sex_excluded_original_rows=sex_excluded,
       quality_at_unmatched=list(A_to_B=qa,B_to_A=qb)),
    input_manifest=m[,setdiff(names(m),c('ancestry_group','ethnicity')),with=FALSE],
    parameters=cfg,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
  out<-file.path(cfg$outdir,paste0('Compare_',gsub('[^a-zA-Z0-9_.-]','_',cfg$compare_id),'.rds'))
  tmp<-paste0(out,'.tmp');saveRDS(result,tmp,compress='xz')
  if(!file.rename(tmp,out)) stop('Failed to move completed RDS into place')
  if(cfg$write_tsv) {
    fwrite(summary,file.path(cfg$outdir,'01_category_summary.tsv'),sep='\t')
    fwrite(loci_classified,file.path(cfg$outdir,'02_loci_classified.tsv'),sep='\t')
    for(item in list(list(name='03_primary_matches',x=primary_matches),
       list(name='04_rescues',x=rescues),list(name='05_review',x=review),
       list(name='06_not_callable',x=not_callable),
       list(name='07_differentials',x=differentials),
       list(name='08_vcf_load_audit',x=vcf_audit),
       list(name='09_primary_edges',x=edges),
       list(name='10_softclip_decisions',x=clip_decisions),
       list(name='11_sex_excluded',x=sex_excluded),
       list(name='12_evidence_reuse_qc',x=reuse_audit),
       list(name='13_plausible_internal_locus_links',x=internal_links),
       list(name='14_duplicate_guard_summary',x=rescue_protection_summary))) {
      if(is.data.frame(item$x)) {
        path <- file.path(cfg$outdir,paste0(item$name,'.tsv'))
        if(!ncol(item$x)) writeLines(character(),path)
        else fwrite(item$x,path,sep='\t',na='NA')
      }
    }
    param<-data.table(parameter=names(cfg),value=vapply(cfg,toString,character(1)))
    fwrite(param,file.path(cfg$outdir,'00_parameters.tsv'),sep='\t')
  }
  message('[COMPARE][DONE] ',out)
  print(summary)
  invisible(result)
}

if(sys.nframe()==0L) {
  conf<-ra_parse_args()
  ra_run(conf)
}
