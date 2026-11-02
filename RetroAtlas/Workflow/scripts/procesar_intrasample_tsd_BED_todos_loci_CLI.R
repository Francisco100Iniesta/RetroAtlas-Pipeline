#!/usr/bin/env Rscript
# Single-sample finalization: retain every primary callable locus, with or without TSD.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(IRanges)})
args <- commandArgs(trailingOnly=TRUE)
data.table::setDTthreads(1L)
opts <- list(window='2000', min_tsd='4', primary_min_support='3', min_svlen='30', min_depth='0', upper_ratio='3')
allowed <- c('rds','tsd','mapping','validation','excluded_predata','helpers','sample','sex','procedence','assembly','out','window','min_tsd','panel_bed','regions','primary_min_support','min_svlen','min_depth','upper_ratio')
if (length(args) %% 2L) stop('Use --option value pairs')
for (i in seq(1L, length(args), by=2L)) {
  key <- gsub('-', '_', sub('^--', '', args[i]), fixed=TRUE)
  if(!key %in% allowed) stop('Unknown option: ', args[i])
  opts[[key]] <- args[i+1L]
}
for (k in setdiff(allowed, c('window','min_tsd','panel_bed','regions','primary_min_support','min_svlen','min_depth','upper_ratio'))) if(is.null(opts[[k]])) stop('Missing --',k)
window <- as.numeric(opts$window); min_tsd <- as.numeric(opts$min_tsd)
if (!is.finite(window) || window < 0 || window != floor(window) || !is.finite(min_tsd) || min_tsd < 1) stop('Invalid window/min-tsd')
source(opts$helpers)
panel <- read_panel(opts$panel_bed)
dir.create(opts$out, recursive=TRUE, showWarnings=FALSE)
prefix <- file.path(opts$out, opts$sample)
obj <- readRDS(opts$rds)
if(!is.list(obj) || !is.data.frame(obj$loci)) stop('Expected RDS with $loci')
loci <- as.data.frame(rbindlist(list(obj$loci,obj$rescue_candidates),fill=TRUE,use.names=TRUE))
if(!all(c('IDs','CHROM','POS') %in% names(loci)) || anyNA(loci$IDs) || anyDuplicated(loci$IDs)) stop('Missing or duplicate locus IDs')
tsd <- fread(opts$tsd, colClasses=list(character=c('event_id','tsd_sequence')))
required <- c('event_id','tsd_detected','tsd_length','tsd_sequence')
if(!all(required %in% names(tsd))) stop('Invalid TSD columns')
if(anyNA(tsd$event_id) || anyDuplicated(tsd$event_id)) stop('Duplicate/missing TSD IDs')
if(any(!tsd$event_id %in% loci$IDs)) stop('TSD contains unknown event IDs')
idx <- match(loci$IDs, tsd$event_id)
detected <- tolower(as.character(tsd$tsd_detected)) %in% c('true','t','1')
lengths <- suppressWarnings(as.numeric(tsd$tsd_length))
valid <- detected & is.finite(lengths) & lengths >= min_tsd
loci$tsd_assessed <- !is.na(idx)
loci$tsd_detected <- detected[idx]
loci$tsd_length <- ifelse(valid[idx] %in% TRUE, lengths[idx], NA_real_)
loci$tsd_sequence <- ifelse(valid[idx] %in% TRUE, as.character(tsd$tsd_sequence[idx]), NA_character_)
mapping <- read_chromosome_map(opts$mapping)
validation <- jsonlite::fromJSON(opts$validation)
if (any(mapping$ASSEMBLY != opts$assembly) || validation$assembly != opts$assembly) stop('Assembly mismatch between mapping, validation and requested output')
normalized <- standardize_loci(loci, mapping, opts$sample, opts$sex, opts$procedence, opts$assembly, panel)
cov <- if(!is.null(opts$regions)) read_coverage(opts$regions,mapping,as.numeric(opts$min_depth),as.numeric(opts$upper_ratio),opts$sex,panel) else NULL
primary_sets <- primary_partition(normalized$kept,cov,as.numeric(opts$primary_min_support),as.numeric(opts$min_svlen))
loci <- primary_sets$loci
obj$rescue_candidates <- primary_sets$rescue_candidates
obj$nocallablecandidates <- primary_sets$nocallablecandidates
obj$parameters$Primary <- list(min_support=as.numeric(opts$primary_min_support),min_svlen=as.numeric(opts$min_svlen),min_depth=as.numeric(opts$min_depth),upper_ratio=as.numeric(opts$upper_ratio))
# Audit both pre-reconstruction exclusions and final loci, including alternate contigs and female Y.
early <- fread(opts$excluded_predata)
excluded <- rbindlist(list(early, normalized$excluded, primary_sets$rescue_candidates), fill=TRUE, use.names=TRUE)
obj$loci <- loci
obj$excluded_loci <- as.data.frame(excluded)
obj$sample_id <- opts$sample
obj$metadata <- list(assembly=opts$assembly, sex=opts$sex, procedence=opts$procedence,
                     adaptive_sampling=!is.null(panel), target_regions=as.data.frame(panel),
                     reference_validation=validation,
                     coordinate_system='1-based; original assembly; no liftover',
                     final_chromosome_column='FINAL_CHR')
obj$chromosome_mapping <- mapping
obj$counts$n_rescue_candidates <- nrow(obj$rescue_candidates)
obj$counts$n_nocallablecandidates <- nrow(obj$nocallablecandidates)
obj$counts$n_loci <- nrow(loci)
obj$counts$n_excluded <- nrow(excluded)
if (is.data.frame(obj$fragment_evidence) && 'CHROM' %in% names(obj$fragment_evidence)) {
  evidence <- standardize_loci(as.data.frame(obj$fragment_evidence), mapping, opts$sample, opts$sex, opts$procedence, opts$assembly, panel)
  obj$fragment_evidence <- evidence$kept
  obj$excluded_fragment_evidence <- evidence$excluded
}
count_flags <- c(n_population_primary_eligible='POPULATION_PRIMARY_ELIGIBLE',
                 n_population_differential_eligible='POPULATION_DIFFERENTIAL_ELIGIBLE',
                 n_evidence_only_aggregates='POPULATION_EVIDENCE_ONLY',
                 n_shadow='SOURCE_SHADOW',n_complex_group='SOURCE_COMPLEX',
                 n_compound_synthetic='FRAGMENT_COMPOUND_SYNTHETIC',
                 n_compound_anchored='FRAGMENT_RECONSTRUCTION_CAUTION')
for(k in names(count_flags)) if(count_flags[[k]] %in% names(loci))
  obj$counts[[k]] <- sum(loci[[count_flags[[k]]]] %in% TRUE)
if(is.data.frame(obj$fragment_evidence))
  obj$counts$n_fragment_evidence_components <- nrow(obj$fragment_evidence)
obj$parameters$Finalization <- list(window=window, min_tsd=min_tsd)
write_flat_tsv(obj$nocallablecandidates,paste0(prefix,'_nocallablecandidates.tsv'))
write_flat_tsv(obj$rescue_candidates,paste0(prefix,'_rescue_candidates.tsv'))
saveRDS(obj, paste0(prefix, '_Final.rds'))
write_flat_tsv(loci, paste0(prefix, '_Final.tsv'))
write_flat_tsv(excluded, paste0(prefix, '_excluded_loci.tsv'))
# Native chromosome IDs remain in BEDs consumed by tools working on the native reference.
events <- loci[,c('SAMPLE_ID','IDs','CHROM','FINAL_CHR','POS','ASSEMBLY','SEX','PROCEDENCE','tsd_length'),drop=FALSE]
events$tiene_tsd_valida <- !is.na(events$tsd_length)
write_flat_tsv(events, paste0(prefix, '_events.tsv'))
if(nrow(events)) {
  limits <- mapping$expected_length[match(events$CHROM, mapping$CHROM)]
  gr <- GRanges(events$CHROM, IRanges(start=pmax(1,events$POS-window), end=pmin(limits,events$POS+window)))
  reduced <- reduce(gr, ignore.strand=TRUE)
  bed <- data.frame(CHROM=as.character(seqnames(reduced)), start=start(reduced)-1L, end=end(reduced), Id=paste0(opts$sample,'_region_',seq_along(reduced)))
  bed <- as.data.frame(clip_to_panel(bed,panel))
  bed$Id <- paste0(opts$sample,'_region_',seq_len(nrow(bed)))
  reduced <- GRanges(bed$CHROM,IRanges(bed$start+1,bed$end))
  hits <- findOverlaps(gr, reduced, ignore.strand=TRUE)
  membership <- cbind(Id=bed$Id[subjectHits(hits)], events[queryHits(hits),,drop=FALSE])
} else {
  bed <- data.frame(CHROM=character(),start=numeric(),end=numeric(),Id=character())
  membership <- cbind(Id=character(),events)
}
fwrite(bed, paste0(prefix,'_regions.bed'), sep='\t', quote=FALSE, col.names=FALSE)
write_flat_tsv(membership, paste0(prefix,'_events_per_region.tsv'))
message('Final: ', nrow(loci), ' loci; excluded: ', nrow(excluded), '; TSD annotated without filtering loci.')
