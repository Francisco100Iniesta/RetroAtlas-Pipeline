#!/usr/bin/env Rscript
# ATLAS · Pairwise Insertion Explorer
# Generates one self-contained HTML report from RetroAtlas Compare D1 RDS.
# Requirements: R >= 4.1 and jsonlite. No external JS/CSS or network needed.
# This is a REPORTER: it does not re-run matching, rescue or genotype inference.

suppressPackageStartupMessages({
  if (!requireNamespace('jsonlite', quietly=TRUE))
    stop('Missing package jsonlite. Install with install.packages("jsonlite").')
})

parse_cli <- function(args) {
  o <- list(rds=NULL,parameters=NULL,glossary=NULL,mane=NULL,
            out='compare_report.html',template=NULL,promoter_bp=2000,
            regions_a=NULL,regions_b=NULL,mapping=NULL,helpers=NULL,catalog=NULL,mane_assembly=NULL)
  i <- 1L
  while (i<=length(args)) {
    arg <- args[[i]]
    if (arg %in% c('-h','--help')) {
      cat('Usage: report_compare.R --rds Compare.rds --regions-a A.regions.bed --regions-b B.regions.bed --mapping chromosomes.tsv --helpers reference_helpers.R --catalog chromosomes.tsv --out report.html [--mane BED --mane-assembly BUILD]\n')
      quit(save='no',status=0)
    }
    if (!startsWith(arg,'--')) stop('Unexpected CLI argument: ',arg)
    if (grepl('=',arg,fixed=TRUE)) {
      a <- strsplit(substring(arg,3L),'=',fixed=TRUE)[[1L]]
      key <- a[[1L]]; val <- paste(a[-1L],collapse='=')
    } else {
      key <- substring(arg,3L)
      if (i==length(args)) stop('Missing value for --',key)
      i <- i+1L; val <- args[[i]]
    }
    key <- gsub('-','_',key,fixed=TRUE)
    if (!key %in% names(o)) stop('Unknown parameter: --',key)
    o[[key]] <- val
    i <- i+1L
  }
  if (is.null(o$rds)) stop('Required argument --rds is missing')
  for (p in c('rds','parameters','glossary','mane','template','regions_a','regions_b','mapping','helpers','catalog')) {
    f <- o[[p]]
    if (!is.null(f) && !file.exists(f)) stop('File not found: ',f)
  }
  if(any(vapply(c('regions_a','regions_b','mapping','helpers','catalog'),function(k)is.null(o[[k]]),logical(1))))
    stop('Require regions-a, regions-b, mapping, helpers, catalog')
  o$promoter_bp <- as.integer(o$promoter_bp)
  if (is.na(o$promoter_bp) || o$promoter_bp<0L) stop('--promoter-bp must be nonnegative')
  o
}
opts <- parse_cli(commandArgs(trailingOnly=TRUE))
dir.create(dirname(opts$out),recursive=TRUE,showWarnings=FALSE)
writeLines('chrom\tEXCLUSION_REASON',paste0(opts$out,'.mane_excluded.tsv'))
cmdline <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
dir_script <- if (length(cmdline)) dirname(normalizePath(sub('^--file=','',cmdline[[1]]),mustWork=FALSE)) else getwd()
template <- if (!is.null(opts$template)) opts$template else file.path(dir_script,'report_template.html')
if (!file.exists(template)) stop('Cannot find template: ',template)
if (is.null(opts$glossary) && file.exists(file.path(dir_script,'COMPARE_D1_GLOSSARY.md')))
  opts$glossary <- file.path(dir_script,'COMPARE_D1_GLOSSARY.md')

cat('[ATLAS] Loading D1 comparison RDS\n')
x <- readRDS(opts$rds)
required <- c('metadata','parameters','primary_matches','rescues','review','not_callable','differentials')
if (!is.list(x) || !all(required %in% names(x)))
  stop('The RDS does not match RetroAtlas Compare D1: missing ',paste(setdiff(required,names(x)),collapse=', '))
source(opts$helpers)
data.table::setDTthreads(1L)
meta <- x$metadata
mapping <- read_chromosome_map(opts$mapping)
if(any(mapping$ASSEMBLY != meta$reference_build)) stop('Report chromosome assembly mismatch')
if(!is.null(opts$mane) && !identical(opts$mane_assembly,meta$reference_build)) stop('MANE assembly mismatch')

if (is.null(meta$sample_A) || is.null(meta$sample_B)) stop('metadata must identify sample_A and sample_B')
sa <- as.character(meta$sample_A[[1]])
sb <- as.character(meta$sample_B[[1]])
if (!nzchar(sa) || !nzchar(sb) || identical(sa,sb)) stop('This reporter requires exactly two different sample IDs')

load_mane <- function(bed,promoter_bp) {
  b <- read_gene_annotation(bed)
  names(b) <- c('chrom','chromStart','chromEnd','name','score','strand',
    'thickStart','thickEnd','rgb','blockCount','blockSizes','chromStarts',
    'transcriptID','cdsStartStat','cdsEndStat','exonFrames','featureType',
    'geneID','geneSymbol','geneType','proteinID','refseqTranscriptID',
    'geneIDNCBI','refseqProteinID','maneStatus')
  b <- b[b$maneStatus %in% c('MANE Select','MANE Plus Clinical','RefSeq Select Curated') &
           !is.na(b$geneSymbol) & nzchar(b$geneSymbol), ,drop=FALSE]
  ck <- function(y) toupper(sub('^chr','',as.character(y),ignore.case=TRUE))
  catalog <- read.delim(opts$catalog, stringsAsFactors=FALSE)
  alias <- setNames(rep(catalog$FINAL_CHR,lengths(strsplit(catalog$aliases,';',fixed=TRUE))),unlist(strsplit(catalog$aliases,';',fixed=TRUE)))
  canonical <- unname(alias[b$chrom])
  if(anyNA(canonical)) {
    # MANE contains alternate loci too; keep these outside the primary report.
    excluded <- b[is.na(canonical),,drop=FALSE]
    excluded$EXCLUSION_REASON <- 'NONPRIMARY_OR_UNKNOWN_CHROMOSOME'
    utils::write.table(excluded,paste0(opts$out,'.mane_excluded.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    warning(nrow(excluded),' nonprimary MANE records excluded; see annotation audit')
  }
  b <- b[!is.na(canonical),,drop=FALSE]
  canonical <- canonical[!is.na(canonical)]
  if(!nrow(b)) stop('No primary MANE records match this assembly')
  limits <- catalog$length[match(canonical,catalog$FINAL_CHR)]
  if(any(!is.finite(as.numeric(b$chromStart)) | as.numeric(b$chromStart)<0 |
         as.numeric(b$chromEnd)>limits | as.numeric(b$chromEnd)<=as.numeric(b$chromStart))) stop('MANE coordinates outside reference')
  b$chrom <- canonical
  b$key <- ck(b$chrom)
  b$priority <- ifelse(b$maneStatus=='MANE Select',0L,1L)
  b <- b[order(b$priority,b$key,b$geneSymbol),,drop=FALSE]
  b <- b[!duplicated(paste(b$key,b$geneSymbol)),,drop=FALSE]
  b$start <- as.numeric(b$chromStart)+1
  b$end <- as.numeric(b$chromEnd)
  b$pstart <- ifelse(b$strand=='+',pmax(1,b$start-promoter_bp),b$end+1)
  b$pend <- ifelse(b$strand=='+',b$start-1,b$end+promoter_bp)
  grouped <- split(seq_len(nrow(b)),b$key)
  parse_blocks <- function(t) {
    if (is.na(t)||!nzchar(t)) return(numeric())
    as.numeric(strsplit(gsub(',$','',t),',',fixed=TRUE)[[1]])
  }
  b$exons <- lapply(seq_len(nrow(b)),function(i) {
    st <- parse_blocks(b$chromStarts[[i]]); ln <- parse_blocks(b$blockSizes[[i]])
    if (!length(st)||length(st)!=length(ln)) return(matrix(numeric(),ncol=2L))
    cbind(as.numeric(b$chromStart[[i]])+st+1,as.numeric(b$chromStart[[i]])+st+ln)
  })
  list(b=b,grouped=grouped)
}
mane_data <- if(is.null(opts$mane)) NULL else load_mane(opts$mane,opts$promoter_bp)

# Optional MANE annotation; MANE Select preferred per gene and Plus Clinical fallback.
annotate_loci <- function(positions, chroms, bed, promoter_bp) {
  n <- length(positions)
  ans <- vector('list',n)
  if (is.null(bed)) {
    for (i in seq_len(n)) ans[[i]] <- list(Gene=NULL,Location=NULL)
    return(ans)
  }
  b <- mane_data$b; grouped <- mane_data$grouped
  ck <- function(y) toupper(sub('^chr','',as.character(y),ignore.case=TRUE))
  for (i in seq_len(n)) {
    pos <- suppressWarnings(as.numeric(positions[[i]]))
    ids <- grouped[[ck(chroms[[i]])]]
    z <- list(Gene=NULL,Location='intergenic')
    if (!length(ids)||!is.finite(pos)) {ans[[i]] <- z;next}
    m <- (pos>=b$start[ids] & pos<=b$end[ids]) |
         (pos>=b$pstart[ids] & pos<=b$pend[ids])
    candidates <- ids[which(m %in% TRUE)]; top <- 0L
    for (j in candidates) {
      inside <- pos>=b$start[[j]] && pos<=b$end[[j]]
      if (inside) {
        e <- b$exons[[j]]
        hit <- if(nrow(e)) which(pos>=e[,1] & pos<=e[,2]) else integer()
        if(length(hit)) {
          en <- if(b$strand[[j]]=='+')hit[[1]] else nrow(e)-hit[[1]]+1L
          loc <- paste0('exon',en);score <- 3L
        } else {loc <- 'intron';score <- 2L}
      } else {loc <- 'promoter';score <- 1L}
      if (score>top) {z <- list(Gene=b$geneSymbol[[j]],Location=loc);top <- score}
    }
    ans[[i]] <- z
  }
  ans
}

as_df <- function(y) {
  if(is.null(y)) return(data.frame())
  as.data.frame(y,stringsAsFactors=FALSE)
}
primary <- as_df(x$primary_matches)
rescues <- as_df(x$rescues)
review <- as_df(x$review)
notcall <- as_df(x$not_callable)
differential <- as_df(x$differentials)
count_vec <- c(PRIMARY_MATCH=nrow(primary),RESCUE=nrow(rescues),REVIEW=nrow(review),
               NOT_CALLABLE=nrow(notcall),DIFFERENTIAL=nrow(differential))
message('[ATLAS] Loci: ',paste(names(count_vec),count_vec,collapse='; '))

# Scalar values for safe, small JSON payload. No full ALT strings are embedded in HTML.
scalar <- function(v) {
  if (is.null(v)||!length(v)) return(NULL)
  if (is.factor(v)) v <- as.character(v)
  if (is.list(v)) return(as.character(jsonlite::toJSON(v,auto_unbox=TRUE)))
  if (!length(v) || all(is.na(v))) return(NULL)
  v <- v[!is.na(v)]
  if (length(v)>1L) return(paste(as.character(v),collapse='; '))
  y <- unname(v[[1L]])
  if (is.numeric(y) && !is.finite(y)) return(NULL)
  y
}
record <- function(d,i,cols) {
  z <- setNames(vector('list',length(cols)),cols)
  for (k in cols) if (k %in% names(d)) z[k] <- list(scalar(d[[k]][i]))
  z
}
primary_fields <- c('pair_id','uid_A','uid_B','Q','Q_LABEL','reason','dpos','size_sim',
  'class_match','event_match','completeness_match','seq_status','seq_identity','seq_coverage','seq_orientation',
  'A_CHROM','A_POS','A_SVLEN','A_GT','A_event_class','A_event','A_completeness',
  'A_alt_support','A_region_status','A_is_synthetic','A_sex_region','A_ploidy',
  'B_CHROM','B_POS','B_SVLEN','B_GT','B_event_class','B_event','B_completeness',
  'B_alt_support','B_region_status','B_is_synthetic','B_sex_region','B_ploidy',
  'genotype_identity','genotype_identity_scope','genotype_identity_reason',
  'GT_normalized_A','GT_normalized_B','primary_GT_ref_conflict')
common_fields <- c('pair_id','source_uid','source_sample','target_sample','category',
  'reason','CHROM','POS','SVLEN','event_class','event','completeness','GT_source',
  'reconstructed_candidate_ids','source_support','source_ploidy','target_ploidy','rescue_type','rescue_gt','rescue_dosage',
  'linked_primary_uid','linked_primary_dpos','near_primary_count','near_low_count',
  'review_subtype','linked_internal_primary_uids','linked_existing_primary_pair_ids',
  'internal_link_evidence','rescue_blocked_by_internal_link','suppressed_rescue_channel',
  'nearby_opposite_primary_uids','nearby_opposite_matched_pair_ids',
  'vcf_status','vcf_reason','vcf_callers','vcf_alt_callers','vcf_n_alt_callers',
  'vcf_gt_conflict','vcf_selected_record_ids','clip_resolution','clip_reason',
  'clip_n_good_alt','clip_good_ids','fragment_status','fragment_ids','fragment_sum_svlen',
  'target_callable','rescue_without_callability','target_mean_depth')

metadata_fields <- c('ASSEMBLY','SEX_A','SEX_B','PROCEDENCE_A','PROCEDENCE_B')
primary_fields <- c(primary_fields,metadata_fields,'A_ORIGINAL_CHROM','B_ORIGINAL_CHROM')
common_fields <- c(common_fields,metadata_fields,'ORIGINAL_CHROM')
mk_primary <- function(d) {
  n <- nrow(d)
  if (!n) return(list())
  ann <- annotate_loci(d$A_POS,d$A_CHROM,opts$mane,opts$promoter_bp)
  lapply(seq_len(n),function(i) {
    z <- record(d,i,primary_fields)
    z$id <- z$pair_id
    z$CHROM <- z$A_CHROM;z$POS <- z$A_POS;z$B_POS <- z$B_POS
    z$SVLEN <- z$A_SVLEN;z$event_class <- z$A_event_class
    z$Gene <- ann[[i]]$Gene;z$Location <- ann[[i]]$Location
    gi <- z$genotype_identity
    z$genotype_identity <- if(is.null(gi)) NULL else as.logical(as.integer(gi))
    z$gt_view <- if(isFALSE(z$genotype_identity))'REVIEW' else if(isTRUE(z$genotype_identity))'MATCH' else 'UNKNOWN'
    z
  })
}
mk_other <- function(d) {
  n <- nrow(d)
  if (!n) return(list())
  ann <- annotate_loci(d$POS,d$CHROM,opts$mane,opts$promoter_bp)
  lapply(seq_len(n),function(i) {
    z <- record(d,i,common_fields)
    z$id <- z$pair_id
    z$Gene <- ann[[i]]$Gene;z$Location <- ann[[i]]$Location
    z
  })
}
message('[ATLAS] Building compact JSON evidence tables')
p <- mk_primary(primary); r <- mk_other(rescues);rv <- mk_other(review)
nc <- mk_other(notcall);di <- mk_other(differential)

read_doc <- function(p) if(!is.null(p))paste(readLines(p,warn=FALSE,encoding='UTF-8'),collapse='\n') else ''
method <- paste(c(
 '# RetroAtlas Pairwise Comparison and Classification',
 'RetroAtlas compares mobile-element insertions between two samples using high-confidence reconstructed insertions, structural variant calls, soft-clipped reads, fragment evidence, and regional coverage.',
 'Insertions independently reconstructed in both samples are classified as **PRIMARY_MATCH**, with matching quality determined by breakpoint proximity, insertion-size similarity, structural annotations, and sequence agreement when evaluated.',
 'Insertions found in only one primary reconstruction set are examined for additional evidence in the other sample. Supported insertions are classified as **RESCUE**, while ambiguous signals, competing nearby loci, and potentially redundant allelic representations are assigned to **REVIEW**.',
 'After searching rescue evidence, regions without sufficient coverage or quality for reliable absence assessment are classified as **NOT_CALLABLE**. Insertions with no supporting evidence in the opposite sample, despite adequate callability and the required negative evidence checks, are classified as **DIFFERENTIAL**.',
 'To reduce false-positive rescues and avoid double-counting, RetroAtlas examines potentially related PRIMARY representations within a configurable 1,000-bp window. Such representations are linked for review rather than automatically merged or counted as independent insertions.',
 '## Genotype identity audit',
 'Genotype identity compares normalized ALT dosage only for positions with comparable expected ploidy. A PRIMARY match with `GT=0/1` in one sample and `GT=1/1` in the other is retained in PRIMARY_MATCH and flagged for secondary genotype review; it is not counted again as a classified REVIEW locus.',
 'PRIMARY matching follows the complete Q10–Q1 tier system. Q10–Q3 require compatible class annotations with tier-specific breakpoint/size limits; Q2/Q1 allow different classes only with concordant ALT sequence alignment.',
 'Matching quality Q does not imply measured sequence identity if sequence alignment was not evaluated. The primary alignment mode and its thresholds are listed in the run parameters.'
),collapse='\n\n')
params <- if (is.list(x$parameters))x$parameters else as.list(x$parameters)
effective <- setNames(lapply(names(params),function(k){
  a <- scalar(params[[k]])
  if(is.null(a))'NA' else as.character(a)
}),names(params))
provided <- list()
if(!is.null(opts$parameters)) {
  q <- utils::read.delim(opts$parameters,stringsAsFactors=FALSE,check.names=FALSE,quote='')
  if(!all(c('parameter','value')%in%names(q)))stop('--parameters must contain parameter and value columns')
  provided <- as.list(stats::setNames(as.character(q$value),as.character(q$parameter)))
}

# Same regional-depth criterion as Single and Compare, over validated chromosomes.
build_quality <- function(opts,x,sa,sb) {
  suppressPackageStartupMessages(library(data.table))
  minimum <- x$parameters$callable_min_depth
  ratio <- x$parameters$callable_max_depth_ratio
  A <- read_coverage(opts$regions_a,mapping,minimum,ratio,if(is.null(meta$sex_A)) 'U' else meta$sex_A)
  B <- read_coverage(opts$regions_b,mapping,minimum,ratio,if(is.null(meta$sex_B)) 'U' else meta$sex_B)
  one <- function(v) v$regions[,.(chr=CHROM,start,end,depth,ratio=depth_ratio,band=depth_band,
                                  call=cover_status=='PASS' & region_status=='PASS')]
  d <- merge(one(A),one(B),by=c('chr','start','end'),suffixes=c('_a','_b'))
  if(!nrow(d)) return(list(present=FALSE,reason='No shared coverage windows'))
  d[,bp:=end-start]
  s <- d[,.(bp=sum(bp),n=.N,callable_A=100*sum(bp[call_a])/sum(bp),
       callable_B=100*sum(bp[call_b])/sum(bp),callable_both=100*sum(bp[call_a & call_b])/sum(bp),
       depth_A=weighted_depth_quantiles(depth_a,bp,.5),depth_B=weighted_depth_quantiles(depth_b,bp,.5),
       max_depth_A=max(depth_a),max_depth_B=max(depth_b),
       lower_A=if(sum(bp[call_a])>0) 100*sum(bp[band_a=='LowerIQR'])/sum(bp[call_a]) else 0,iqr_A=if(sum(bp[call_a])>0) 100*sum(bp[band_a=='IQR'])/sum(bp[call_a]) else 0,
       upper_A=if(sum(bp[call_a])>0) 100*sum(bp[band_a=='UpperIQR'])/sum(bp[call_a]) else 0,
       lower_B=if(sum(bp[call_b])>0) 100*sum(bp[band_b=='LowerIQR'])/sum(bp[call_b]) else 0,iqr_B=if(sum(bp[call_b])>0) 100*sum(bp[band_b=='IQR'])/sum(bp[call_b]) else 0,
       upper_B=if(sum(bp[call_b])>0) 100*sum(bp[band_b=='UpperIQR'])/sum(bp[call_b]) else 0),by=chr]
  s[,delta_pp:=callable_B-callable_A]
  s[,length_bp:=mapping$expected_length[match(chr,mapping$CHROM)]]
  s[,evaluated_pct:=100*bp/length_bp]
  full <- sum(d$bp)
  summary <- list(bp=full,windows=nrow(d),callable_A=100*sum(d$bp[d$call_a])/full,
    callable_B=100*sum(d$bp[d$call_b])/full,callable_both=100*sum(d$bp[d$call_a & d$call_b])/full,
    disagree=100*sum(d$bp[d$call_a!=d$call_b])/full)
  # Deterministic positions and discordances for rendering; summaries use all windows.
  setorder(d,chr,start)
  d[,deviation:=abs(depth_a-depth_b)]
  ix <- sort(unique(c(as.integer(round(seq(1,nrow(d),length.out=min(4000,nrow(d))))),
                      order(d$deviation,decreasing=TRUE)[seq_len(min(1000,nrow(d)))])))
  pts <- lapply(ix,function(i)unname(list(d$chr[i],d$start[i]/1e6,d$depth_a[i],d$depth_b[i],d$ratio_a[i],d$ratio_b[i],as.integer(d$call_a[i]),as.integer(d$call_b[i]))))
  list(present=TRUE,sample_A=sa,sample_B=sb,genome=summary,chr=as.data.frame(s),scatter=pts,
    params=list(min_depth=minimum,max_depth_ratio=ratio,baseline_A=A$baseline,baseline_B=B$baseline),
    inputs=list(A=basename(opts$regions_a),B=basename(opts$regions_b)),
    unpaired_windows_A=nrow(A$regions)-nrow(d),unpaired_windows_B=nrow(B$regions)-nrow(d),
    formula=paste0('Coverage is normalized to the autosomal median. In males, X and Y each use their own chromosome median; in females X uses the autosomal baseline and Y is excluded. Callable regions require ',
                   if(minimum>0) paste0('coverage ≥ ',minimum,'×') else 'strictly positive coverage',
                   ' and coverage ≤ ',ratio,' times the applicable median. Outside these limits, regions are not callable.',
                   ' Normalization medians use all evaluated regions; IQR quartiles and percentages use callable bases only.',
                   ' Statistics use the evaluated regions only; in Adaptive Sampling these are restricted to the target panel.'))
}

counts <- as.list(as.integer(count_vec));names(counts) <- names(count_vec)
quality <- build_quality(opts,x,sa,sb)
output <- list(
  meta=list(sample_A=sa,sample_B=sb,compare_id=scalar(meta$compare_id),
            reference_build=scalar(meta$reference_build),created_at=scalar(meta$created_at),
            sex_A=meta$sex_A,sex_B=meta$sex_B,procedence_A=meta$procedence_A,procedence_B=meta$procedence_B,
            chromosome_lengths=meta$chromosome_lengths,mane_available=!is.null(opts$mane),rescue_policy=meta$rescue_policy),
  counts=counts,primary=p,rescue=r,review=rv,not_callable=nc,differential=di,
  effective_params=effective,provided_params=provided,
  method=method,glossary=read_doc(opts$glossary),quality=quality
)
json <- jsonlite::toJSON(output,auto_unbox=TRUE,null='null',na='null',digits=10,
                         POSIXt='ISO8601',force=TRUE)
enc <- jsonlite::base64_enc(charToRaw(enc2utf8(json)))
# jsonlite::base64_enc() wraps long output; raw line breaks inside JS string literals break Chrome.
enc <- gsub('[[:space:]]', '', enc)
page <- paste(readLines(template,warn=FALSE,encoding='UTF-8'),collapse='\n')
if(!grepl('__BASE64__',page,fixed=TRUE))stop('Template missing __BASE64__ placeholder')
page <- sub('__BASE64__',enc,page,fixed=TRUE)
dir.create(dirname(opts$out),recursive=TRUE,showWarnings=FALSE)
writeLines(enc2utf8(page),opts$out,useBytes=TRUE)
message('[ATLAS] HTML: ',normalizePath(opts$out,mustWork=FALSE))
message('[ATLAS] Genotype audit flags: ',sum(primary$genotype_identity==0,na.rm=TRUE))
