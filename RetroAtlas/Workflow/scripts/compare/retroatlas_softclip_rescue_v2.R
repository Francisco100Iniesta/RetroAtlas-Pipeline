#!/usr/bin/env Rscript
# RetroAtlas V2: three-state soft-clip ALT evidence classifier.
# Output is explicitly CANDIDATE-level: no molecular genotype confirmation.
# Resolution: RESCUE_CANDIDATE / REVIEW / NO_CLIP_SUPPORT.
# Compact interface: retroatlas_clip_classify(...) -> one-row, three-column data.frame.
# For audit: retroatlas_clip_classify(..., detailed=TRUE) -> list(decision, audit, clips, clusters...).
# Two execution modes:
#   source('retroatlas_softclip_rescue_v2.R')
#   salida <- evaluate_clip_rescue(chr='chr15', pos=29561003, alt_seq=alt,
#                                  query_clips='muestra_B.detail.tsv.gz')
#   salida$summary  # data.frame de una fila; salida$clips y salida$clusters
# CLI:
#   Rscript retroatlas_softclip_rescue_v2.R --chr chr15 --pos 29561003 \
#     --alt-rds Twin_0040_intra_final_TSD.rds --alt-id Atlas_5192 \
#     --query-clips muestra_B.softclips_ge50.detail.tsv.gz \
#     --source-clips muestra_A.softclips_ge50.detail.tsv.gz \
#     --window 5000 --output-dir resultados
#
# LIMITATIONS: CLIP_ID does not represent QNAME or independent molecules; no MAPQ,
# CIGAR, genomic junction anchoring or genome-wide ALT specificity is available._CANDIDATE
# RESCUE_CANDIDATE is computational only: validate against BAM/CRAM before altering GT.
# NO_CLIP_SUPPORT does not imply ALT absence or a confirmed differential.

suppressPackageStartupMessages({
  if (!requireNamespace('Biostrings', quietly=TRUE) ||
      !requireNamespace('pwalign', quietly=TRUE) ||
      !requireNamespace('IRanges', quietly=TRUE)) {
    stop('Instala BiocManager::install(c("Biostrings", "pwalign", "IRanges"))', call.=FALSE)
  }
})

ra_defaults <- function() list(
  window=5000L, core_window=1000L, cluster_bp=10L,
  shared_peak_bp=25L, rescue_max_distance=1000L,
  min_clip_len=50L, match=2, mismatch=-3,
  gap_open=5, gap_extend=2,
  min_identity=90, min_clip_coverage=90, min_aligned_bp=100L,
  partial_identity=95, partial_coverage=80, partial_aligned_bp=200L,
  rescue_min_clips=3L, rescue_min_cluster=3L,
  rescue_min_unique_sequences=2L,
  min_alt_union_bp=200L, alt_end_bp=30L,
  near_full_alt_pct=80,
  min_sequence_entropy=1.2,
  max_align_cells=80000000,
  read_chunk=10000L,
  review_center_bp=100L, review_center_min_peak=4L,
  review_pair_min_peak=3L, review_source_min_peak=2L,
  review_pair_min_separation=20L, review_pair_max_separation=5000L,
  shifted_rescue_window=5000L,
  competitor_score_ratio=0.97, competitor_min_ambiguous_clips=2L
)

ra_chr <- function(chr) {
  chr <- as.character(chr)
  if (length(chr)!=1L || is.na(chr) || !nzchar(chr)) stop('Invalid chromosome')
  if (!grepl('^chr',chr,ignore.case=TRUE)) chr <- paste0('chr',chr)
  sub('^chr','chr',chr,ignore.case=TRUE)
}

ra_clean_seq <- function(x, what='secuencia') {
  x <- toupper(gsub('[[:space:]]', '', as.character(x)))
  if (length(x)!=1L || is.na(x) || !nzchar(x)) stop(what, ' is empty or missing')
  if (grepl('[^ACGTRYSWKMBDHVN]',x)) stop(what, ' contains invalid IUPAC symbols')
  x
}

ra_entropy <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  nuc <- strsplit(toupper(x),'',fixed=TRUE)[[1L]]
  nuc <- nuc[nuc %in% c('A','C','G','T')]
  if (!length(nuc)) return(NA_real_)
  p <- as.numeric(table(nuc))/length(nuc)
  -sum(p*log2(p))
}

ra_empty_clips <- function() {
  data.frame(CLIP_ID=character(),CHR=character(),POS0=numeric(),
    POS1=numeric(),SOFTCLIP_LEN=numeric(),SEQUENCE=character(),
    stringsAsFactors=FALSE)
}

# Supports Megadepth TSV(.gz) and generated FASTA(.gz)
# with |chr:POS1|pos0=...|len=... in the header.
# Retains only query-region sequences in memory.
read_softclip_region <- function(path, chr, pos, window=5000L,
                                 chunk=10000L) {
  if (!file.exists(path)) stop('Soft-clip file not found: ',path)
  chr <- ra_chr(chr); pos <- as.numeric(pos); window <- as.numeric(window)
  if (!is.finite(pos) || !is.finite(window) || window<0) stop('Invalid POS or window')
  con <- if (grepl('\\.gz$',path,ignore.case=TRUE)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con),add=TRUE)
  out <- list(); nout <- 0L
  fasta <- grepl('\\.(fa|fasta)(\\.gz)?$',path,ignore.case=TRUE)
  if (fasta) {
    head <- NULL; seq_lines <- character()
    flush <- function() {
      if (is.null(head)) return(NULL)
      # Cabecera: >CLIP_ID|chr15:29561003|pos0=29561002|len=123
      parts <- strsplit(sub('^>','',head),'\\|')[[1L]]
      loc <- parts[grepl('^chr[^:]+:[0-9]+$',parts)]
      if (!length(loc)) return(NULL)
      p <- strsplit(loc[[1]],':',fixed=TRUE)[[1L]]
      if (ra_chr(p[[1]])!=chr) return(NULL)
      pos1 <- as.numeric(p[[2]])
      if (!is.finite(pos1) || abs(pos1-pos)>window) return(NULL)
      sq <- toupper(paste(seq_lines,collapse=''))
      if (!nzchar(sq)) return(NULL)
      data.frame(CLIP_ID=parts[[1]],CHR=chr,POS0=pos1-1,
        POS1=pos1,SOFTCLIP_LEN=nchar(sq),SEQUENCE=sq,
        stringsAsFactors=FALSE)
    }
    repeat {
      lines <- readLines(con,n=chunk,warn=FALSE)
      if (!length(lines)) break
      for (ln in lines) {
        if (startsWith(ln,'>')) {
          x <- flush()
          if (!is.null(x)) {nout <- nout+1L;out[[nout]] <- x}
          head <- ln;seq_lines <- character()
        } else if (!is.null(head)) seq_lines <- c(seq_lines,trimws(ln))
      }
    }
    x <- flush(); if (!is.null(x)) {nout <- nout+1L;out[[nout]] <- x}
  } else {
    header <- readLines(con,n=1L,warn=FALSE)
    if (!length(header)) stop('Empty soft-clip file: ',path)
    fields <- strsplit(header,'\t',fixed=TRUE)[[1L]]
    needed <- c('CLIP_ID','CHR','POS0','POS1','SOFTCLIP_LEN','SEQUENCE')
    if (!all(needed %in% fields)) stop('Missing soft-clip columns: ',
       paste(setdiff(needed,fields),collapse=','))
    idx <- match(needed,fields)
    repeat {
      lines <- readLines(con,n=chunk,warn=FALSE)
      if (!length(lines)) break
      for (ln in lines) {
        # Keep regional sequences only.
        parts <- strsplit(ln,'\t',fixed=TRUE)[[1L]]
        if (length(parts)<max(idx)) next
        if (ra_chr(parts[[idx[[2]]]])!=chr) next
        pos1 <- suppressWarnings(as.numeric(parts[[idx[[4]]]]))
        if (!is.finite(pos1) || abs(pos1-pos)>window) next
        nout <- nout+1L
        out[[nout]] <- data.frame(
          CLIP_ID=parts[[idx[[1]]]], CHR=chr,
          POS0=suppressWarnings(as.numeric(parts[[idx[[3]]]])),
          POS1=pos1,
          SOFTCLIP_LEN=suppressWarnings(as.numeric(parts[[idx[[5]]]])),
          SEQUENCE=toupper(parts[[idx[[6]]]]),stringsAsFactors=FALSE)
      }
    }
  }
  if (!length(out)) return(ra_empty_clips())
  ans <- do.call(rbind,out)
  rownames(ans) <- NULL
  ans
}

# For multiple loci, loading TSV once is faster but consumes more RAM.
# y pasarla como query_table / source_table a evaluate_clip_rescue().
load_softclip_table <- function(path) {
  if (grepl('\\.(fa|fasta)(\\.gz)?$',path,ignore.case=TRUE))
    stop('Use read_softclip_region for FASTA; whole FASTA is not loaded into RAM')
  con <- if (grepl('\\.gz$',path,ignore.case=TRUE)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con))
  a <- utils::read.delim(con,sep='\t',header=TRUE,quote='',comment.char='',
                         stringsAsFactors=FALSE,check.names=FALSE)
  a
}

ra_select_table <- function(df,chr,pos,window) {
  if (!is.data.frame(df)) stop('query_table/source_table must be a data.frame')
  needed <- c('CLIP_ID','CHR','POS0','POS1','SOFTCLIP_LEN','SEQUENCE')
  if (!all(needed %in% names(df))) stop('Table missing columns: ',paste(setdiff(needed,names(df)),collapse=','))
  chroms <- as.character(df$CHR)
  chroms <- ifelse(grepl('^chr',chroms,ignore.case=TRUE),chroms,paste0('chr',chroms))
  pos1 <- suppressWarnings(as.numeric(df$POS1))
  idx <- which(!is.na(chroms) & chroms==chr & is.finite(pos1) & abs(pos1-pos)<=window)
  if (!length(idx)) return(ra_empty_clips())
  ans <- df[idx,needed,drop=FALSE]
  ans$POS1 <- as.numeric(ans$POS1)
  ans$POS0 <- as.numeric(ans$POS0)
  ans$SOFTCLIP_LEN <- as.numeric(ans$SOFTCLIP_LEN)
  ans$SEQUENCE <- as.character(ans$SEQUENCE)
  ans
}

ra_empty_alignment <- function() {
  data.frame(strand=NA_character_,score=NA_real_,alt_start=NA_integer_,
    alt_end=NA_integer_,clip_start=NA_integer_,clip_end=NA_integer_,
    n_aligned_clip=NA_integer_,n_aligned_alt=NA_integer_,
    n_compared=NA_integer_,n_matches=NA_integer_,n_mismatches=NA_integer_,
    n_gap_clip=NA_integer_,n_gap_alt=NA_integer_,
    identity_pct=NA_real_,clip_coverage_pct=NA_real_,
    alt_coverage_pct=NA_real_,clip_entropy=NA_real_,
    alignment_error=NA_character_,stringsAsFactors=FALSE)
}

# ALT is always the subject; clip always the pattern. Never exchange by length.
# Both clip orientations (forward/reverse complement).
compare_clip_alt <- function(alt_seq,clip_seq,opts=ra_defaults()) {
  result <- ra_empty_alignment()
  alt <- ra_clean_seq(alt_seq,'ALT')
  clip <- ra_clean_seq(clip_seq,'clip')
  result$clip_entropy <- ra_entropy(clip)
  if (as.double(nchar(alt))*as.double(nchar(clip)) > opts$max_align_cells) {
    result$alignment_error <- 'TOO_MANY_DP_CELLS'
    return(result)
  }
  tryCatch({
    m <- pwalign::nucleotideSubstitutionMatrix(
      match=opts$match,mismatch=opts$mismatch,baseOnly=FALSE)
    a <- Biostrings::DNAString(alt)
    b <- Biostrings::DNAString(clip)
    x <- pwalign::pairwiseAlignment(pattern=b,subject=a,type='local',
      substitutionMatrix=m,gapOpening=opts$gap_open,gapExtension=opts$gap_extend)
    y <- pwalign::pairwiseAlignment(pattern=Biostrings::reverseComplement(b),
      subject=a,type='local',substitutionMatrix=m,
      gapOpening=opts$gap_open,gapExtension=opts$gap_extend)
    rev <- as.numeric(pwalign::score(y))>as.numeric(pwalign::score(x))
    aln <- if (rev) y else x
    result$strand <- if (rev) '-' else '+'
    result$score <- as.numeric(pwalign::score(aln))
    pp <- strsplit(as.character(pwalign::alignedPattern(aln)),'',fixed=TRUE)[[1L]]
    ss <- strsplit(as.character(pwalign::alignedSubject(aln)),'',fixed=TRUE)[[1L]]
    if (!length(pp) || length(pp)!=length(ss)) stop('alineamiento vacio')
    valid <- pp!='-' & ss!='-'
    result$n_compared <- sum(valid)
    result$n_matches <- sum(pp==ss & valid)
    result$n_mismatches <- sum(pp!=ss & valid)
    result$n_aligned_clip <- sum(pp!='-')
    result$n_aligned_alt <- sum(ss!='-')
    result$n_gap_clip <- sum(pp=='-' & ss!='-')
    result$n_gap_alt <- sum(ss=='-' & pp!='-')
    result$identity_pct <- if (result$n_compared>0L)
      100*result$n_matches/result$n_compared else NA_real_
    result$clip_coverage_pct <- 100*result$n_aligned_clip/nchar(clip)
    result$alt_coverage_pct <- 100*result$n_aligned_alt/nchar(alt)
    # Coordenadas 1-based del ALT y clip orientado: start/end de XStringViews
    result$alt_start <- as.integer(IRanges::start(pwalign::subject(aln)))
    result$alt_end <- as.integer(IRanges::end(pwalign::subject(aln)))
    cp1 <- as.integer(IRanges::start(pwalign::pattern(aln)))
    cp2 <- as.integer(IRanges::end(pwalign::pattern(aln)))
    if (rev) {
      result$clip_start <- nchar(clip)-cp2+1L
      result$clip_end <- nchar(clip)-cp1+1L
    } else {
      result$clip_start <- cp1;result$clip_end <- cp2
    }
    result
  }, error=function(e) {
    result$alignment_error <- paste0('ALIGNMENT_FAILED: ',conditionMessage(e))
    result
  })
}

ra_add_scores <- function(df,alt_seq,opts,pos) {
  if (!nrow(df)) {
    result <- cbind(df,ra_empty_alignment()[FALSE,,drop=FALSE])
    result$strict <- logical();result$partial_long <- logical()
    result$alt_match <- logical();result$DIST_TO_POS <- numeric()
    result$ALT_END <- logical()
    return(result)
  }
  df <- df[!is.na(df$SEQUENCE) & nzchar(df$SEQUENCE) &
             is.finite(df$POS1) & nchar(df$SEQUENCE)>=opts$min_clip_len,,drop=FALSE]
  if (!nrow(df)) return(ra_add_scores(ra_empty_clips(),alt_seq,opts,pos))
  aligned <- lapply(df$SEQUENCE,function(s) compare_clip_alt(alt_seq,s,opts))
  z <- cbind(df,do.call(rbind,aligned))
  z$DIST_TO_POS <- as.numeric(z$POS1)-as.numeric(pos)
  z$strict <- !is.na(z$identity_pct) &
    z$identity_pct >= opts$min_identity &
    z$clip_coverage_pct >= opts$min_clip_coverage &
    z$n_compared >= min(opts$min_aligned_bp,nchar(alt_seq)) &
    !is.na(z$clip_entropy) & z$clip_entropy>=opts$min_sequence_entropy
  z$partial_long <- !z$strict & !is.na(z$identity_pct) &
    z$identity_pct >= opts$partial_identity &
    z$clip_coverage_pct >= opts$partial_coverage &
    z$n_compared >= min(opts$partial_aligned_bp,nchar(alt_seq)) &
    !is.na(z$clip_entropy) & z$clip_entropy>=opts$min_sequence_entropy
  z$alt_match <- z$strict | z$partial_long
  z$ALT_END <- z$alt_match & !is.na(z$alt_start) & !is.na(z$alt_end) &
     (z$alt_start<=opts$alt_end_bp |
        z$alt_end>=nchar(alt_seq)-opts$alt_end_bp+1L)
  rownames(z) <- NULL
  z
}

# Optional diagnostic: does a different ALT score almost equally well on the
# same clips? This does NOT compare an insertion to the REF allele; for that,
# inspect genome-ALT junctions with QNAME/CIGAR in the original BAM.
ra_competitor_ambiguity <- function(q,competitor_alts,opts) {
  if (is.null(competitor_alts) || !length(competitor_alts) || !nrow(q))
    return(rep(FALSE,nrow(q)))
  comp <- unique(vapply(competitor_alts,ra_clean_seq,character(1L),'competitor ALT'))
  mat <- pwalign::nucleotideSubstitutionMatrix(
      match=opts$match,mismatch=opts$mismatch,baseOnly=FALSE)
  flags <- rep(FALSE,nrow(q))
  good <- which(q$alt_match & !is.na(q$score))
  for (i in good) {
    clip <- Biostrings::DNAString(q$SEQUENCE[[i]])
    rc <- Biostrings::reverseComplement(clip)
    for (seq in comp) {
      if (nchar(seq)*nchar(q$SEQUENCE[[i]])>opts$max_align_cells) next
      alt <- Biostrings::DNAString(seq)
      # Only score is needed, not a full alignment.
      f <- pwalign::pairwiseAlignment(clip,alt,type='local',
          substitutionMatrix=mat,gapOpening=opts$gap_open,
          gapExtension=opts$gap_extend,scoreOnly=TRUE)
      r <- pwalign::pairwiseAlignment(rc,alt,type='local',
          substitutionMatrix=mat,gapOpening=opts$gap_open,
          gapExtension=opts$gap_extend,scoreOnly=TRUE)
      if (max(as.numeric(f),as.numeric(r)) >=
          opts$competitor_score_ratio*q$score[[i]]) {
        flags[[i]] <- TRUE
        break
      }
    }
  }
  flags
}

# Maximum-span clustering, not transitive single-linkage.
ra_clusters <- function(df,pos,opts) {
  output <- data.frame(cluster_id=integer(),start_pos=numeric(),end_pos=numeric(),
   center_pos=numeric(),dist_to_pos=numeric(),n_total=integer(),n_strict=integer(),
   n_partial=integer(),n_alt=integer(),max_exact=integer(),
   median_identity=double(),stringsAsFactors=FALSE)
  if (!nrow(df)) return(output)
  d <- df[order(df$POS1),,drop=FALSE]
  start_idx <- 1L; items <- list(); id <- 0L
  for (i in seq_len(nrow(d)+1L)) {
    split <- i>nrow(d) || (i>start_idx && d$POS1[i]-d$POS1[start_idx]>opts$cluster_bp)
    if (!split) next
    v <- d[start_idx:(i-1L),,drop=FALSE]
    id <- id+1L
    mid <- stats::median(v$POS1)
    vals <- v$identity_pct[v$alt_match & !is.na(v$identity_pct)]
    items[[id]] <- data.frame(cluster_id=id,
       start_pos=min(v$POS1),end_pos=max(v$POS1),center_pos=mid,
       dist_to_pos=mid-pos,n_total=nrow(v),
       n_strict=sum(v$strict),n_partial=sum(v$partial_long),
       n_alt=sum(v$alt_match),max_exact=max(table(v$POS1)),
       median_identity=if(length(vals)) stats::median(vals) else NA_real_)
    start_idx <- i
  }
  result <- do.call(rbind,items)
  rownames(result) <- NULL
  result
}

ra_alt_coverage <- function(df,alt_len) {
  d <- df[df$alt_match & is.finite(df$alt_start) & is.finite(df$alt_end),,drop=FALSE]
  if (!nrow(d)) return(0L)
  starts <- pmax(1L,pmin(alt_len,as.integer(d$alt_start)))
  ends <- pmax(1L,pmin(alt_len,as.integer(d$alt_end)))
  as.integer(sum(IRanges::width(IRanges::reduce(IRanges::IRanges(start=starts,end=ends)))))
}

ra_shared_peaks <- function(cl_q,cl_s,tolerance=25L,alt_only=TRUE) {
  if (alt_only) {
    cl_q <- cl_q[cl_q$n_alt>0,,drop=FALSE]
    cl_s <- cl_s[cl_s$n_alt>0,,drop=FALSE]
  } else {
    cl_q <- cl_q[cl_q$n_total>=2,,drop=FALSE]
    cl_s <- cl_s[cl_s$n_total>=2,,drop=FALSE]
  }
  if (!nrow(cl_q) || !nrow(cl_s)) return(0L)
  sum(vapply(cl_q$center_pos,function(p) any(abs(cl_s$center_pos-p)<=tolerance),logical(1)))
}

# Detailed result: list(summary=one row, clips=table, clusters=table,
#                       source_clips=optional table, source_clusters=optional table).
# RESCUE_CANDIDATE never implies a confirmed genotype.
evaluate_clip_rescue <- function(chr,pos,alt_seq,
                                 query_clips=NULL,source_clips=NULL,
                                 query_table=NULL,source_table=NULL,
                                 opts=list(), competitor_alts=NULL,
                                 sample_id=NA_character_,
                                 source_id=NA_character_,alt_id=NA_character_) {
  cfg <- utils::modifyList(ra_defaults(),opts)
  chr <- ra_chr(chr);pos <- as.numeric(pos)
  if (length(pos)!=1L || !is.finite(pos) || pos<1 || pos!=floor(pos)) stop('POS must be a positive 1-based integer')
  alt <- ra_clean_seq(alt_seq,'ALT')
  if (is.null(query_table) && is.null(query_clips)) stop('Provide query_clips or query_table')
  if (cfg$window<cfg$core_window) stop('window must be >= core_window')
  for (k in c('window','core_window','cluster_bp','shared_peak_bp','rescue_max_distance',
              'min_clip_len','min_aligned_bp','partial_aligned_bp','rescue_min_clips',
              'rescue_min_cluster','rescue_min_unique_sequences',
              'min_alt_union_bp','alt_end_bp','read_chunk',
              'review_center_bp','review_center_min_peak',
              'review_pair_min_peak','review_source_min_peak',
              'review_pair_min_separation','review_pair_max_separation',
              'shifted_rescue_window','competitor_min_ambiguous_clips')) {
    if (!is.finite(cfg[[k]]) || cfg[[k]]<0) stop('Invalid parameter: ',k)
  }
  if (!is.finite(cfg$max_align_cells) || cfg$max_align_cells<1) stop('max_align_cells must be >0')
  input_q <- if (!is.null(query_table)) ra_select_table(query_table,chr,pos,cfg$window) else
    read_softclip_region(query_clips,chr,pos,cfg$window,cfg$read_chunk)
  q <- ra_add_scores(input_q,alt,cfg,pos)
  ambiguous <- ra_competitor_ambiguity(q,competitor_alts,cfg)
  q$competitor_similar_score <- ambiguous
  n_competitor_ambiguous <- sum(ambiguous)
  competitor_blocks <- n_competitor_ambiguous >=cfg$competitor_min_ambiguous_clips
  cq <- ra_clusters(q,pos,cfg)
  s <- ra_add_scores(ra_empty_clips(),alt,cfg,pos)
  cs <- ra_clusters(s,pos,cfg)
  has_source <- !is.null(source_clips) || !is.null(source_table)
  if (has_source) {
    input_s <- if (!is.null(source_table)) ra_select_table(source_table,chr,pos,cfg$window) else
       read_softclip_region(source_clips,chr,pos,cfg$window,cfg$read_chunk)
    s <- ra_add_scores(input_s,alt,cfg,pos)
    cs <- ra_clusters(s,pos,cfg)
  }
  n_strict <- sum(q$strict)
  n_partial <- sum(q$partial_long)
  n_good <- sum(q$alt_match)
  n_unique_alt_seq <- length(unique(q$SEQUENCE[q$alt_match]))
  n_alt_positions <- length(unique(q$POS1[q$alt_match]))
  n_raw_core <- sum(abs(q$DIST_TO_POS)<=cfg$core_window)
  close_cl <- cq[abs(cq$dist_to_pos)<=cfg$rescue_max_distance,,drop=FALSE]
  max_close <- if (nrow(close_cl)) max(close_cl$n_alt) else 0L
  max_cluster <- if (nrow(cq)) max(cq$n_alt) else 0L
  raw_peak <- if (nrow(cq)) max(cq$n_total) else 0L
  n_good_clusters <- sum(cq$n_alt>0)
  covered_bp <- ra_alt_coverage(q,nchar(alt))
  end_support <- sum(q$ALT_END)
  n_near_full_alt <- sum(q$alt_match & !is.na(q$alt_coverage_pct) &
                         q$alt_coverage_pct>=cfg$near_full_alt_pct)
  shared_good <- if(has_source) ra_shared_peaks(cq,cs,cfg$shared_peak_bp,TRUE) else NA_integer_
  shared_raw <- if(has_source) ra_shared_peaks(cq,cs,cfg$shared_peak_bp,FALSE) else NA_integer_
  n_source_good <- if(has_source) sum(s$alt_match) else NA_integer_
  n_failed <- sum(!is.na(q$alignment_error))
  n_source_failed <- if(has_source) sum(!is.na(s$alignment_error)) else NA_integer_
  # Spatial review is deliberately stricter than generic clipping.
  # A single strong near-breakpoint peak, or two peaks on opposite sides.
  # When A is available, displaced pairs must reproduce in A (+/- tolerance).
  matched_to_source <- function(cluster_row, min_source=cfg$review_source_min_peak) {
    if (!has_source || !nrow(cs)) return(FALSE)
    any(cs$n_total>=min_source &
        abs(cs$center_pos-cluster_row$center_pos)<=cfg$shared_peak_bp)
  }
  central_peaks <- cq[cq$n_total>=cfg$review_center_min_peak &
                        abs(cq$dist_to_pos)<=cfg$review_center_bp,,drop=FALSE]
  strong_central <- nrow(central_peaks)>0L
  pair_candidates <- cq[cq$n_total>=cfg$review_pair_min_peak &
                          abs(cq$dist_to_pos)<=cfg$window,,drop=FALSE]
  pair_support <- FALSE
  n_paired_structural_peaks <- 0L
  if (nrow(pair_candidates)) {
    paired <- vapply(seq_len(nrow(pair_candidates)),function(i)
       matched_to_source(pair_candidates[i,,drop=FALSE]),logical(1))
    n_paired_structural_peaks <- sum(paired)
    if (nrow(pair_candidates)>=2L) {
      for (i in seq_len(nrow(pair_candidates)-1L)) {
        for (j in seq.int(i+1L,nrow(pair_candidates))) {
          left <- pair_candidates[i,,drop=FALSE]
          right <- pair_candidates[j,,drop=FALSE]
          separation <- abs(left$center_pos-right$center_pos)
          opposite <- (left$dist_to_pos * right$dist_to_pos < 0)
          spatial_pair <- opposite &&
            separation>=cfg$review_pair_min_separation &&
            separation<=cfg$review_pair_max_separation
          # Without A: two independently prominent clusters required.
          strong_without_source <- !has_source &&
               left$n_total>=max(4,cfg$review_pair_min_peak) &&
               right$n_total>=max(4,cfg$review_pair_min_peak)
          pair_with_source <- has_source && paired[i] && paired[j]
          if (spatial_pair && (strong_without_source || pair_with_source))
            pair_support <- TRUE
        }
      }
    }
  }
  structural_review <- strong_central || pair_support

  # Rescue: alignments against ALT + compact cluster + meaningful ALT span
  # and an ALT end. Shifted clusters (> core-window) are accepted ONLY when
  # a recurrent corresponding peak exists in the source sample.
  enough_alt <- covered_bp>=min(cfg$min_alt_union_bp,nchar(alt))
  clustered <- max_close>=cfg$rescue_min_cluster
  spatial_reconciled <- is.na(shared_good) || n_source_good<3L || shared_good>=1L
  rescue_base <- n_good>=cfg$rescue_min_clips &&
    n_unique_alt_seq>=cfg$rescue_min_unique_sequences && enough_alt &&
    end_support>=1L
  core_rescue <- rescue_base && clustered && spatial_reconciled
  shifted_good <- cq[cq$n_alt>=cfg$rescue_min_cluster &
     abs(cq$dist_to_pos)>cfg$rescue_max_distance &
     abs(cq$dist_to_pos)<=min(cfg$window,cfg$shifted_rescue_window),,drop=FALSE]
  shifted_reproduced <- FALSE
  if (has_source && nrow(shifted_good)>0L)
    shifted_reproduced <- any(vapply(seq_len(nrow(shifted_good)),
      function(i) matched_to_source(shifted_good[i,,drop=FALSE]),logical(1)))
  shifted_rescue <- rescue_base && shifted_reproduced

  if ((core_rescue || shifted_rescue) && !competitor_blocks) {
    status <- 'RESCUE_CANDIDATE'
    reason <- if (core_rescue)
      'ALT_MATCHES_AND_CLUSTER_NEAR_BREAKPOINT'
    else 'ALT_MATCHES_SHIFTED_CLUSTER_REPRODUCED_IN_SOURCE'
  } else if (n_good>0L) {
    status <- 'REVIEW'
    if (competitor_blocks) reason <- 'ALT_MATCH_NOT_SPECIFIC_VS_COMPETITOR'
    else if (n_good<cfg$rescue_min_clips) reason <- 'ALT_MATCHES_TOO_FEW'
    else if (n_unique_alt_seq<cfg$rescue_min_unique_sequences)
      reason <- 'ALT_MATCHES_REDUNDANT_SEQUENCES'
    else if (end_support==0L) reason <- 'ALT_MATCHES_INTERNAL_ONLY'
    else if (!enough_alt) reason <- 'ALT_MATCHES_INSUFFICIENT_SPAN'
    else if (!clustered && !shifted_reproduced)
      reason <- 'ALT_MATCHES_DISPERSED_OR_SHIFTED_UNVERIFIED'
    else if (!spatial_reconciled) reason <- 'SOURCE_QUERY_PEAKS_DISCORDANT'
    else reason <- 'ALT_MATCHES_NOT_SUFFICIENT_FOR_RESCUE'
  } else if (n_failed>0L) {
    # Cannot interpret missing matches if some local alignments were skipped.
    status <- 'REVIEW'
    reason <- 'ALIGNMENTS_UNEVALUATED'
  } else if (strong_central) {
    status <- 'REVIEW';reason <- 'STRONG_CENTRAL_CLIPPING_WITHOUT_ALT_MATCH'
  } else if (pair_support) {
    status <- 'REVIEW';reason <- 'REPRODUCED_TWO_PEAK_STRUCTURE_WITHOUT_ALT_MATCH'
  } else {
    status <- 'NO_CLIP_SUPPORT'
    reason <- if (nrow(q)==0L) 'NO_REGIONAL_SOFTCLIPS' else 'NO_ALT_OR_CLEAR_SPATIAL_STRUCTURE'
  }
  # Never infer ALT absence here; candidate rescues need genome-ALT junction validation.
  summary <- data.frame(
    sample_id=as.character(sample_id),source_id=as.character(source_id),
    alt_id=as.character(alt_id), CHR=chr, POS=pos, alt_length=nchar(alt),
    window_bp=cfg$window,core_window_bp=cfg$core_window,
    n_clips_window=nrow(q),n_clips_core=n_raw_core,
    n_alt_strict=n_strict,n_alt_partial=n_partial,n_alt_total=n_good,
    n_unique_alt_clip_sequences=n_unique_alt_seq,
    n_distinct_alt_clip_positions=n_alt_positions,
    alt_match_fraction=if(nrow(q)) n_good/nrow(q) else NA_real_,
    n_alt_clusters=n_good_clusters,max_alt_cluster=max_cluster,
    max_alt_cluster_near_pos=max_close,max_raw_cluster=raw_peak,
    n_alt_end_clips=end_support,n_near_full_alt=n_near_full_alt,
    alt_union_covered_bp=covered_bp,
    alt_union_covered_pct=100*covered_bp/nchar(alt),
    nearest_alt_clip_bp=if(n_good) min(abs(q$DIST_TO_POS[q$alt_match])) else NA_real_,
    median_alt_identity_pct=if(n_good) stats::median(q$identity_pct[q$alt_match]) else NA_real_,
    median_alt_clip_coverage_pct=if(n_good) stats::median(q$clip_coverage_pct[q$alt_match]) else NA_real_,
    max_alt_compared_bp=if(n_good) max(q$n_compared[q$alt_match]) else NA_real_,
    n_alignment_skipped_or_failed=n_failed,
    n_alt_clips_competitor_ambiguous=n_competitor_ambiguous,
    n_source_clips=if(has_source) nrow(s) else NA_integer_,
    n_source_alt_matches=n_source_good,
    n_source_alignment_failed=n_source_failed,
    n_shared_alt_peaks=shared_good,n_shared_raw_peaks=shared_raw,
    status=status,reason=reason,
    interpretation=if(status=='RESCUE_CANDIDATE')
        'Candidate insertion support; confirm junction, read independence and sequence specificity before assigning genotype'
      else if(status=='REVIEW')
        'Evidence requires review; do not promote to high-confidence differential from this module'
      else 'No supporting clip/ALT evidence detected; biological absence is NOT established',
    n_reproduced_raw_peaks=n_paired_structural_peaks,
    review_strong_central=strong_central, review_two_peaks=pair_support,
    rescue_shifted_verified=shifted_reproduced,
    genotype_confirmed=FALSE,
    stringsAsFactors=FALSE)
  list(summary=summary, clips=q,clusters=cq,source_clips=s,source_clusters=cs,
       parameters=cfg)
}

# Public API: EXACTLY three columns, one row, convenient for bind_rows().
# Set detailed=TRUE for a list(decision, audit, clips, clusters, ...).
retroatlas_clip_classify <- function(chr,pos,alt_seq,
                                    query_clips=NULL,source_clips=NULL,
                                    query_table=NULL,source_table=NULL,
                                    params=list(), competitor_alts=NULL,
                                    detailed=FALSE) {
  full <- evaluate_clip_rescue(chr=chr,pos=pos,alt_seq=alt_seq,
       query_clips=query_clips,source_clips=source_clips,
       query_table=query_table,source_table=source_table,
       opts=params, competitor_alts=competitor_alts)
  decision <- data.frame(
    resolution=as.character(full$summary$status[[1L]]),
    reason=as.character(full$summary$reason[[1L]]),
    n_good_alt_clips=as.integer(full$summary$n_alt_total[[1L]]),
    stringsAsFactors=FALSE
  )
  if (!detailed) return(decision)
  c(list(decision=decision, audit=full$summary),full[c('clips','clusters',
      'source_clips','source_clusters','parameters')])
}

ra_load_alt <- function(alt_seq=NULL,alt_fasta=NULL,alt_rds=NULL,
                        alt_id=NULL,chr=NULL,pos=NULL) {
  n <- sum(!vapply(list(alt_seq,alt_fasta,alt_rds),is.null,logical(1)))
  if(n!=1L) stop('Usa exactamente UNA: --alt-seq, --alt-fasta o --alt-rds')
  if(!is.null(alt_seq)) return(ra_clean_seq(alt_seq))
  if(!is.null(alt_fasta)) {
    if(!file.exists(alt_fasta)) stop('FASTA not found: ',alt_fasta)
    con <- if(grepl('\\.gz$',alt_fasta)) gzfile(alt_fasta,'rt') else file(alt_fasta,'rt')
    on.exit(close(con))
    x <- readLines(con,warn=FALSE)
    h <- which(startsWith(x,'>'))
    if(length(h)!=1L) stop('ALT FASTA must contain exactly ONE sequence')
    return(ra_clean_seq(paste(x[-h],collapse='')))
  }
  obj <- readRDS(alt_rds)
  loci <- if (is.list(obj) && !is.null(obj$loci)) obj$loci else obj
  if(!is.data.frame(loci) || !all(c('CHROM','POS','ALT') %in% names(loci)))
    stop('RDS must contain loci with CHROM, POS, ALT')
  if(!is.null(alt_id)) {
    id_cols <- intersect(c('ID','IDs','event_id','SOURCE_IDS'),names(loci))
    if(!length(id_cols)) stop('RDS has no identifier column')
    index <- Reduce('|', lapply(id_cols, function(nm) {
      vals <- as.character(loci[[nm]])
      !is.na(vals) & vals==as.character(alt_id)
    }))
    loci <- loci[index,,drop=FALSE]
  } else {
    loci <- loci[as.character(loci$CHROM)==ra_chr(chr) &
                    as.numeric(loci$POS)==as.numeric(pos),,drop=FALSE]
  }
  if(nrow(loci)!=1L) stop('ALT lookup is not unique; use --alt-id or check CHR/POS (matches=',nrow(loci),')')
  ra_clean_seq(loci$ALT[[1L]])
}

ra_help <- function() {
  cat('Usage:\n',
' Rscript retroatlas_softclip_rescue_v2.R --chr chr15 --pos 29561003 \\\n',
'   --alt-seq ACGT... --query-clips B.detail.tsv.gz [opciones]\n',
' Rscript retroatlas_softclip_rescue_v2.R --chr chr15 --pos 29561003 \\\n',
'   --alt-rds Twin_0040.rds --alt-id Atlas_5192 \\\n',
'   --query-clips B.detail.tsv.gz --source-clips A.detail.tsv.gz \\\n',
'   --output-dir resultados\n\n',
' ALT: --alt-seq / --alt-fasta (una entrada) / --alt-rds y --alt-id opcional\n',
' Obligatorios: --chr --pos --query-clips y una fuente ALT\n',
' Parameters (defaults):\n',
' --window 5000 --core-window 1000 --cluster-bp 10\n',
' --shared-peak-bp 25 --rescue-max-distance 1000\n',
' --match 2 --mismatch -3 --gap-open 5 --gap-extend 2\n',
' --min-identity 90 --min-clip-coverage 90 --min-aligned-bp 100\n',
' --partial-identity 95 --partial-coverage 80 --partial-aligned-bp 200\n',
' --rescue-min-clips 3 --rescue-min-cluster 3\n',
' --rescue-min-unique-sequences 2\n',

' --min-alt-union-bp 200 --alt-end-bp 30 --near-full-alt-pct 80\n',
' --min-sequence-entropy 1.2 --max-align-cells 80000000\n',
' --review-center-bp 100 --review-center-min-peak 4\n',
' --review-pair-min-peak 3 --review-source-min-peak 2\n',
' --shifted-rescue-window 5000\n',
' --competitor-fasta alternative_alts.fa (optional, multi-FASTA)\n',
' --competitor-score-ratio 0.97 --competitor-min-ambiguous-clips 2\n',
' --output-dir . --sample-id B --source-id A --alt-id Atlas_X\n\n',
' Outputs: *.decision.tsv (3 cols), *.summary.tsv, *.clips.tsv, *.clusters.tsv, *.result.rds;\n',
'            optional *.source_clips.tsv, *.source_clusters.tsv\n',
' In R: source(...); retroatlas_clip_classify(...) -> data.frame (3 cols)\n',
' NOTE: RESCUE_CANDIDATE is not a confirmed GT; NO_CLIP_SUPPORT is not absence.\n',sep='')
}

ra_parse_cli <- function(args) {
  if ('--help' %in% args || '-h' %in% args) {ra_help();return(NULL)}
  out <- list();i <- 1L
  while(i<=length(args)) {
    a <- args[[i]]
    if(!startsWith(a,'--')) stop('Unexpected flag: ',a)
    z <- strsplit(sub('^--','',a),'=',fixed=TRUE)[[1L]]
    key <- z[1L]
    if(length(z)>1L) value <- paste(z[-1L],collapse='=') else {
      if(i==length(args)) stop('Missing value for --',key)
      value <- args[[i+1L]];i <- i+1L
    }
    if(key %in% names(out)) stop('Duplicate argument --',key)
    out[[key]] <- value;i <- i+1L
  }
  out
}

ra_cli <- function(args=commandArgs(trailingOnly=TRUE)) {
  args <- ra_parse_cli(args)
  if(is.null(args)) return(invisible(NULL))
  need <- c('chr','pos','query-clips')
  if(!all(need %in% names(args))) stop('Missing required arguments: ',paste(setdiff(need,names(args)),collapse=','))
  cfg <- ra_defaults()
  names_cli <- gsub('_','-',names(cfg),fixed=TRUE)
  names_cli[names(cfg)=='max_align_cells'] <- 'max-align-cells'
  aliases <- c('chr','pos','query-clips','source-clips','alt-seq','alt-fasta',
               'alt-rds','alt-id','sample-id','source-id','output-dir','competitor-fasta')
  unknown <- setdiff(names(args),c(aliases,names_cli))
  if(length(unknown)) stop('Unknown flags: ',paste(unknown,collapse=','))
  for(j in seq_along(cfg)) {
    label <- names_cli[j]
    if(!is.null(args[[label]])) {
      v <- suppressWarnings(as.numeric(args[[label]]))
      if(length(v)!=1L || !is.finite(v)) stop('Invalid numeric value: --',label)
      cfg[[j]] <- v
    }
  }
  alt <- ra_load_alt(alt_seq=args[['alt-seq']],alt_fasta=args[['alt-fasta']],
      alt_rds=args[['alt-rds']],alt_id=args[['alt-id']],chr=args$chr,pos=args$pos)
  competitor_alts <- NULL
  if(!is.null(args[['competitor-fasta']])) {
    if(!file.exists(args[['competitor-fasta']])) stop('Competitor FASTA not found')
    competitor_alts <- as.character(Biostrings::readDNAStringSet(args[['competitor-fasta']]))
  }
  result <- evaluate_clip_rescue(chr=args$chr,pos=as.numeric(args$pos),
    alt_seq=alt,query_clips=args[['query-clips']],
    source_clips=args[['source-clips']],opts=cfg,
    competitor_alts=competitor_alts,
    sample_id=if(is.null(args[['sample-id']])) NA_character_ else args[['sample-id']],
    source_id=if(is.null(args[['source-id']])) NA_character_ else args[['source-id']],
    alt_id=if(is.null(args[['alt-id']])) NA_character_ else args[['alt-id']])
  dir <- if(is.null(args[['output-dir']])) '.' else args[['output-dir']]
  dir.create(dir,showWarnings=FALSE,recursive=TRUE)
  if(!dir.exists(dir)) stop('Cannot create output-dir: ',dir)
  prefix <- file.path(dir,paste0('softclip_',gsub('[^A-Za-z0-9._-]','_',ra_chr(args$chr)),
                                 '_',as.integer(args$pos)))
  compact <- result$summary[,c('status','reason','n_alt_total'),drop=FALSE]
  names(compact) <- c('resolution','reason','n_good_alt_clips')
  utils::write.table(compact,paste0(prefix,'.decision.tsv'),
     sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
  utils::write.table(result$summary,paste0(prefix,'.summary.tsv'),
     sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
  utils::write.table(result$clips,paste0(prefix,'.clips.tsv'),
     sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
  utils::write.table(result$clusters,paste0(prefix,'.clusters.tsv'),
     sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
  if(!is.null(args[['source-clips']])) {
    utils::write.table(result$source_clips,paste0(prefix,'.source_clips.tsv'),
       sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
    utils::write.table(result$source_clusters,paste0(prefix,'.source_clusters.tsv'),
       sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE,na='NA')
  }
  saveRDS(result,paste0(prefix,'.result.rds'))
  print(compact,row.names=FALSE)
  invisible(result)
}

# source() no ejecuta la CLI; Rscript script.R si.
if (sys.nframe() == 0L) {
  tryCatch(ra_cli(),error=function(e){
    message('ERROR: ',conditionMessage(e))
    quit(save='no',status=1L)
  })
}
