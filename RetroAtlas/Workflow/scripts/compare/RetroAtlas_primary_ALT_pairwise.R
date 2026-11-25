#!/usr/bin/env Rscript
# RetroAtlas: full-length PRIMARY ALT vs PRIMARY ALT sequence concordance.
# Uses the SAME pwalign/Biostrings alignment primitives as
# retroatlas_softclip_rescue_v2.R, but GLOBAL alignment instead of LOCAL.
# This is a standalone QC module; it does NOT change PRIMARY cluster identity.
# Dependencies: pwalign, Biostrings.
#
# Public functions:
#   retroatlas_alt_pairwise(alt_a, alt_b, params=list(), detailed=FALSE)
#   retroatlas_primary_pair_qc(alt_a, alt_b, class_a, class_b, event_a, event_b,
#                             completeness_a, completeness_b, params=list(),
#                             detailed=FALSE)
#
# Percentages are 0-100. Primary ALTs should be fully reconstructed INSERTED
# sequences, not VCF alleles including reference anchor bases or genomic flanks.

ra_alt_pw_defaults <- function() list(
  min_identity = 95,
  min_reciprocal_coverage = 90,
  min_compared_bp = 50L,
  min_informative_fraction = 0.80,
  match = 2,
  mismatch = -3,
  gap_open = 5,
  gap_extend = 2,
  max_align_cells = 80000000,
  check_reverse_complement = TRUE,
  only_on_annotation_conflict = TRUE,
  check_on_class_conflict = TRUE,
  check_on_event_conflict = FALSE,
  check_on_completeness_conflict = FALSE
)

ra_alt_pw_clean <- function(x, label) {
  if (length(x) != 1L || is.na(x)) stop(label, ': expected one ALT sequence')
  x <- toupper(gsub('[[:space:]]', '', as.character(x)))
  if (!nzchar(x) || !grepl('^[ACGTRYSWKMBDHVN]+$', x))
    stop(label, ': missing/invalid IUPAC DNA sequence (or symbolic ALT)')
  x
}

ra_alt_pw_row <- function(id_a, id_b, len_a=NA_integer_, len_b=NA_integer_) {
  data.frame(
    id_a=as.character(id_a), id_b=as.character(id_b),
    len_a=as.integer(len_a), len_b=as.integer(len_b),
    orientation=NA_character_, aligner='pwalign', mode='global',
    score=NA_real_, alignment_columns=NA_integer_,
    paired_bp=NA_integer_, informative_bp=NA_integer_,
    matches=NA_integer_, mismatches=NA_integer_,
    ambiguous_paired_bp=NA_integer_, gap_bases_a=NA_integer_,
    gap_bases_b=NA_integer_, gap_events=NA_integer_,
    identity_pct=NA_real_, identity_with_gaps_pct=NA_real_,
    coverage_a_pct=NA_real_, coverage_b_pct=NA_real_,
    reciprocal_coverage_min_pct=NA_real_, informative_fraction=NA_real_,
    seq_qc='SEQ_NOT_EVALUABLE', qc_reason=NA_character_,
    alignment_error=NA_character_, stringsAsFactors=FALSE
  )
}

ra_alt_pw_gapruns <- function(x) {
  if (!length(x)) return(0L)
  r <- rle(x == '-')
  as.integer(sum(r$values))
}

ra_alt_pw_from_strings <- function(ans, aa, bb, score, orientation, cfg) {
  # aa is ALT A aligned, bb is ALT B aligned (possibly reverse-complemented).
  pa <- strsplit(aa, '', fixed=TRUE)[[1L]]
  pb <- strsplit(bb, '', fixed=TRUE)[[1L]]
  if (length(pa) != length(pb) || !length(pa)) stop('Invalid aligned strings')
  paired <- pa != '-' & pb != '-'
  acgt <- c('A','C','G','T')
  informative <- paired & pa %in% acgt & pb %in% acgt
  n_info <- sum(informative)
  n_match <- sum(informative & pa == pb)
  n_mismatch <- sum(informative & pa != pb)
  n_pair <- sum(paired)
  n_columns <- length(pa)

  ans$orientation <- orientation
  ans$score <- as.numeric(score)
  ans$alignment_columns <- as.integer(n_columns)
  ans$paired_bp <- as.integer(n_pair)
  ans$informative_bp <- as.integer(n_info)
  ans$matches <- as.integer(n_match)
  ans$mismatches <- as.integer(n_mismatch)
  ans$ambiguous_paired_bp <- as.integer(n_pair - n_info)
  # A-gap means bases unique to B; B-gap means bases unique to A.
  ans$gap_bases_a <- as.integer(sum(pa == '-' & pb != '-'))
  ans$gap_bases_b <- as.integer(sum(pb == '-' & pa != '-'))
  ans$gap_events <- as.integer(ra_alt_pw_gapruns(pa) + ra_alt_pw_gapruns(pb))
  ans$identity_pct <- if(n_info) 100*n_match/n_info else NA_real_
  ans$identity_with_gaps_pct <- 100*n_match/n_columns
  # GLOBAL coverage is not tautologically 100%: count only bases paired
  # to another nucleotide, NOT those opposite a gap.
  ans$coverage_a_pct <- 100*n_pair/ans$len_a
  ans$coverage_b_pct <- 100*n_pair/ans$len_b
  ans$reciprocal_coverage_min_pct <- min(ans$coverage_a_pct,ans$coverage_b_pct)
  ans$informative_fraction <- if(n_pair) n_info/n_pair else 0

  enough_info <- n_info >= min(cfg$min_compared_bp, ans$len_a, ans$len_b)
  if (!enough_info || ans$informative_fraction < cfg$min_informative_fraction) {
    ans$seq_qc <- 'SEQ_NOT_EVALUABLE'
    ans$qc_reason <- 'INSUFFICIENT_UNAMBIGUOUS_BASES'
  } else if (ans$identity_pct >= cfg$min_identity &&
             ans$reciprocal_coverage_min_pct >= cfg$min_reciprocal_coverage) {
    ans$seq_qc <- 'SEQ_CONCORDANT'
    ans$qc_reason <- 'HIGH_IDENTITY_AND_RECIPROCAL_COVERAGE'
  } else if (ans$identity_pct >= cfg$min_identity) {
    ans$seq_qc <- 'SEQ_PARTIAL'
    ans$qc_reason <- 'HIGH_IDENTITY_BUT_LOW_RECIPROCAL_COVERAGE'
  } else {
    ans$seq_qc <- 'SEQ_DISCORDANT'
    ans$qc_reason <- 'IDENTITY_BELOW_THRESHOLD'
  }
  ans
}

# Compare fully reconstructed ALT sequences of two PRIMARY observations.
# By default uses GLOBAL (Needleman-Wunsch) alignment in both orientations.
# Does not automatically fuse or split loci. A/B always refer to original ALTs.
retroatlas_alt_pairwise <- function(alt_a, alt_b,
                                    params=list(), id_a=NA_character_,
                                    id_b=NA_character_, detailed=FALSE) {
  cfg <- utils::modifyList(ra_alt_pw_defaults(), params)
  allowed <- names(ra_alt_pw_defaults())
  unknown <- setdiff(names(params),allowed)
  if(length(unknown)) stop('Unknown parameter(s): ',paste(unknown,collapse=', '))
  for(nm in c('min_identity','min_reciprocal_coverage')) {
    if(!is.numeric(cfg[[nm]]) || length(cfg[[nm]])!=1L ||
       !is.finite(cfg[[nm]]) || cfg[[nm]]<0 || cfg[[nm]]>100)
      stop('Invalid ',nm)
  }
  if (!is.finite(cfg$max_align_cells) || cfg$max_align_cells < 1 ||
      !is.finite(cfg$min_compared_bp) || cfg$min_compared_bp < 1 ||
      !is.finite(cfg$min_informative_fraction) || cfg$min_informative_fraction < 0 ||
      cfg$min_informative_fraction>1) stop('Invalid alignment QC limits')
  if (!is.logical(cfg$check_reverse_complement) || length(cfg$check_reverse_complement)!=1L ||
      is.na(cfg$check_reverse_complement)) stop('Invalid check_reverse_complement')

  # Invalid or synthetic ALT is NOT a negative sequence match.
  a <- tryCatch(ra_alt_pw_clean(alt_a,'ALT A'), error=function(e) e)
  b <- tryCatch(ra_alt_pw_clean(alt_b,'ALT B'), error=function(e) e)
  if(inherits(a,'error') || inherits(b,'error')) {
    ans <- ra_alt_pw_row(id_a,id_b)
    ans$qc_reason <- 'ALT_MISSING_OR_INVALID'
    errs <- character()
    if (inherits(a,'error')) errs <- c(errs,conditionMessage(a))
    if (inherits(b,'error')) errs <- c(errs,conditionMessage(b))
    ans$alignment_error <- paste(errs,collapse='; ')
    return(if(detailed) list(summary=ans, aligned_a=NULL, aligned_b=NULL,
                             parameters=cfg) else ans)
  }
  ans <- ra_alt_pw_row(id_a,id_b,nchar(a),nchar(b))
  if(as.double(nchar(a))*as.double(nchar(b)) > cfg$max_align_cells) {
    ans$qc_reason <- 'TOO_MANY_DP_CELLS'
    ans$alignment_error <- 'Alignment skipped by max_align_cells guard'
    return(if(detailed) list(summary=ans,aligned_a=NULL,aligned_b=NULL,
                             parameters=cfg) else ans)
  }
  if(!requireNamespace('pwalign',quietly=TRUE) || !requireNamespace('Biostrings',quietly=TRUE))
    stop('Install Bioconductor packages pwalign and Biostrings')

  aligned_a <- NULL
  aligned_b <- NULL
  orientation <- '+'
  aln_score <- NA_real_
  error <- NULL
  tryCatch({
    rc_b <- if(isTRUE(cfg$check_reverse_complement))
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(b))) else NULL

    # Exact string equality is common inside clusters: avoid quadratic DP.
    if (identical(a,b) || (!is.null(rc_b) && identical(a,rc_b))) {
      orientation <- if(identical(a,b)) '+' else '-'
      aligned_a <- a
      aligned_b <- a
      aln_score <- NA_real_ # no DP run; do not invent score
    } else {
      substitution <- pwalign::nucleotideSubstitutionMatrix(
        match=cfg$match,mismatch=cfg$mismatch,baseOnly=FALSE)
      a_dna <- Biostrings::DNAString(a)
      b_dna <- Biostrings::DNAString(b)
      # Try both orientations. Only store the winning alignment.
      fw <- pwalign::pairwiseAlignment(pattern=b_dna,subject=a_dna,
        type='global',substitutionMatrix=substitution,
        gapOpening=cfg$gap_open,gapExtension=cfg$gap_extend)
      if (!is.null(rc_b)) {
        rv <- pwalign::pairwiseAlignment(pattern=Biostrings::DNAString(rc_b),
          subject=a_dna,type='global',substitutionMatrix=substitution,
          gapOpening=cfg$gap_open,gapExtension=cfg$gap_extend)
        if(as.numeric(pwalign::score(rv)) > as.numeric(pwalign::score(fw))) {
          best <- rv; orientation <- '-'
        } else {
          best <- fw; orientation <- '+'
        }
      } else best <- fw
      aln_score <- as.numeric(pwalign::score(best))
      aligned_a <- as.character(pwalign::alignedSubject(best))
      aligned_b <- as.character(pwalign::alignedPattern(best))
    }
    ans <- ra_alt_pw_from_strings(ans,aligned_a,aligned_b,aln_score,orientation,cfg)
  },error=function(e) { error <<- conditionMessage(e) })
  if (!is.null(error)) {
    ans$seq_qc <- 'SEQ_NOT_EVALUABLE'
    ans$qc_reason <- 'ALIGNMENT_FAILED'
    ans$alignment_error <- error
    aligned_a <- NULL; aligned_b <- NULL
  }
  if(detailed)
    list(summary=ans,aligned_a=aligned_a,aligned_b=aligned_b,parameters=cfg)
  else ans
}

ra_alt_pw_annot_equal <- function(a,b) {
  if(length(a)!=1L || length(b)!=1L || is.na(a) || is.na(b)) return(NA)
  a <- trimws(as.character(a)); b <- trimws(as.character(b))
  if(!nzchar(a) || !nzchar(b)) return(NA)
  identical(toupper(a),toupper(b))
}

ra_alt_pw_complex <- function(event, architecture) {
  x <- paste(event,architecture,collapse=' ')
  grepl('embedded|mosaic|multi|compound|complex|inver|nested|chimeric',
        x,ignore.case=TRUE)
}

# Annotation-conflict-aware wrapper for PRIMARY candidate pairs.
# Geometry (position and size) MUST be screened by the existing merge first.
# A sequence match is supporting evidence, not proof of a single allelic event.
retroatlas_primary_pair_qc <- function(
    alt_a,alt_b,class_a=NA_character_,class_b=NA_character_,
    event_a=NA_character_,event_b=NA_character_,
    completeness_a=NA_character_,completeness_b=NA_character_,
    architecture_a=NA_character_,architecture_b=NA_character_,
    id_a=NA_character_,id_b=NA_character_,params=list(),detailed=FALSE) {
  cfg <- utils::modifyList(ra_alt_pw_defaults(), params)
  c_equal <- ra_alt_pw_annot_equal(class_a,class_b)
  e_equal <- ra_alt_pw_annot_equal(event_a,event_b)
  k_equal <- ra_alt_pw_annot_equal(completeness_a,completeness_b)
  conflicts <- c(class=identical(c_equal,FALSE),
                 event=identical(e_equal,FALSE),
                 completeness=identical(k_equal,FALSE))
  complex <- ra_alt_pw_complex(event_a,architecture_a) ||
             ra_alt_pw_complex(event_b,architecture_b)
  # By default inspect CLASS disagreement; event/completeness are reported
  # in the population consensus, not automatically re-aligned.
  check <- !isTRUE(cfg$only_on_annotation_conflict) ||
    (isTRUE(cfg$check_on_class_conflict) && conflicts[['class']]) ||
    (isTRUE(cfg$check_on_event_conflict) && conflicts[['event']]) ||
    (isTRUE(cfg$check_on_completeness_conflict) && conflicts[['completeness']])
  if (check) {
    x <- retroatlas_alt_pairwise(alt_a,alt_b,params=params,id_a=id_a,
                                 id_b=id_b,detailed=detailed)
    z <- if(detailed) x$summary else x
  } else {
    z <- ra_alt_pw_row(id_a,id_b)
    z$seq_qc <- 'SEQ_NOT_TESTED'
    z$qc_reason <- 'NO_ANNOTATION_CONFLICT'
    x <- if(detailed) list(summary=z,aligned_a=NULL,aligned_b=NULL,
                            parameters=cfg) else z
  }
  z$class_conflict <- conflicts[['class']]
  z$event_conflict <- conflicts[['event']]
  z$completeness_conflict <- conflicts[['completeness']]
  z$complex_architecture <- complex
  z$merge_sequence_decision <- if (complex && any(conflicts))
    'COMPLEX_ARCHITECTURE_REVIEW' else if (!check) 'GEOMETRY_ONLY' else if(z$seq_qc == 'SEQ_CONCORDANT')
    'SEQUENCE_SUPPORTS_GEOMETRIC_MERGE' else if(z$seq_qc == 'SEQ_DISCORDANT')
    'SEQUENCE_DISCORDANCE_REVIEW' else 'SEQUENCE_UNRESOLVED_REVIEW'
  if(detailed) {x$summary <- z;x} else z
}
