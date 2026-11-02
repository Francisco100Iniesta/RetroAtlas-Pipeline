#!/usr/bin/env Rscript
# Canonical coordinates for analysis; native coordinates for indexed BAM tools.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
spec <- jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
source(spec$helpers)
out <- spec$out
 dir.create(out, recursive=TRUE, showWarnings=FALSE)
paths <- function(sid, suffix) file.path(out,paste0(sid,suffix))
union <- data.table(CHROM=character(),start=numeric(),end=numeric())
allmaps <- list(); manifest <- list()
for(s in spec$samples) {
  x <- readRDS(s$rds)
  if(!is.list(x) || is.null(x$metadata) || !identical(x$metadata$assembly,spec$assembly))
    stop('Compare requires finalized Single RDS in the selected assembly: ',s$id)
  if(!identical(x$metadata$sex,s$sex) || !identical(x$metadata$procedence,s$procedence))
    stop('Sample metadata changed since Single; rerun finalization: ',s$id)
  m <- read_chromosome_map(s$mapping)
  panel <- read_panel(s$panel)
  if(!identical(isTRUE(x$metadata$adaptive_sampling),!is.null(panel)))
    stop('Adaptive Sampling changed since Single; rerun finalization: ',s$id)
  if(!is.null(panel) && !isTRUE(all.equal(as.data.frame(panel),
      as.data.frame(x$metadata$target_regions),check.attributes=FALSE)))
    stop('Target panel changed since Single; rerun finalization: ',s$id)
  if(!all(panel_membership(x$loci,panel))) stop('Finalized loci outside target panel: ',s$id)
  if(any(m$ASSEMBLY != spec$assembly)) stop('Chromosome map assembly mismatch')
  if(anyDuplicated(m$FINAL_CHR[m$included])) stop('Ambiguous chromosome mapping')
  allmaps[[s$id]] <- m
  canonical <- function(d) {
    d <- as.data.frame(d)
    if(!nrow(d)) return(d)
    if(!'CHROM' %in% names(d)) stop('CHROM missing')
    j <- match(d$CHROM,m$CHROM)
    if(anyNA(j) || any(!m$included[j])) stop('Unvalidated chromosome in finalized RDS')
    if('FINAL_CHR' %in% names(d) && any(d$FINAL_CHR != m$FINAL_CHR[j])) stop('FINAL_CHR mismatch')
    d$ORIGINAL_CHROM <- d$CHROM
    d$CHROM <- d$FINAL_CHR <- m$FINAL_CHR[j]
    d
  }
  cov <- read_coverage(s$coverage,m,spec$min_depth,spec$upper_ratio,s$sex,panel)
  # Never trust saved PASS flags after thresholds or representative positions change.
  gates <- primary_partition(x$loci,cov,if(is.null(spec$primary_min_support)) 3 else spec$primary_min_support,
                             if(is.null(spec$min_svlen)) 30 else spec$min_svlen)
  x$loci <- canonical(gates$loci)
  pool <- as.data.frame(rbindlist(list(x$rescue_candidates,gates$rescue_candidates),fill=TRUE,use.names=TRUE))
  if(nrow(pool)) {
    pool <- standardize_loci(pool,m,s$id,s$sex,s$procedence,spec$assembly,panel)$kept
    pool <- annotate_coverage(pool,cov)
  }
  x$rescue_candidates <- canonical(pool)
  x$nocallablecandidates <- if(nrow(x$rescue_candidates)) x$rescue_candidates[
    !(x$rescue_candidates$cover_status %in% 'PASS' & x$rescue_candidates$region_status %in% 'PASS'),,drop=FALSE] else x$rescue_candidates
  if(!is.null(x$fragment_evidence) && nrow(x$fragment_evidence)) x$fragment_evidence <- canonical(x$fragment_evidence)
  saveRDS(x,paths(s$id,'.Final.rds'))
  fwrite(cov$regions[,.(CHROM=FINAL_CHR,start,end,depth)],paths(s$id,'.regions.bed'),sep='\t',col.names=FALSE)
  if(file.info(s$bed)$size>0) {
    b <- fread(s$bed,header=FALSE,select=1:3,col.names=c('CHROM','start','end'),
               colClasses=list(character=1L,numeric=2:3))
    j <- match(b$CHROM,m$CHROM)
    if(anyNA(j)||any(!m$included[j])||any(b$start<0|b$end<=b$start|b$end>m$expected_length[j])) stop('Invalid Single BED')
    pad <- if(is.null(spec$search_padding)) 0 else spec$search_padding
    b[,`:=`(start=pmax(0,start-pad),end=pmin(m$expected_length[j],end+pad))]
    b <- clip_to_panel(b,panel)
    j <- match(b$CHROM,m$CHROM)
    b[,CHROM:=m$FINAL_CHR[j]]
    union <- rbind(union,b)
  }
  panel_path <- paths(s$id,'.panel.native.bed')
  if(is.null(panel)) unlink(panel_path) else fwrite(panel,panel_path,sep='\t',col.names=FALSE)
  manifest[[s$id]] <- data.table(sample_id=s$id,sex=s$sex,procedence=s$procedence,origin=s$procedence,assembly=spec$assembly,
    panel_bed=if(is.null(panel)) '' else normalizePath(panel_path),
    final_rds=normalizePath(paths(s$id,'.Final.rds')),quality_regions=normalizePath(paths(s$id,'.regions.bed')))
  # Optional population metadata; the same preparation supports two or N samples.
  for(key in c('sex_karyotype'))
    if(!is.null(s[[key]])) manifest[[s$id]][,(key):=as.character(s[[key]])]
}
merge_bed <- function(b) {
  if(!nrow(b)) return(b)
  setorder(b,CHROM,start,end)
  b[,grp:=cumsum(start>shift(cummax(end),fill=-1)),by=CHROM]
  b[,.(start=min(start),end=max(end)),by=.(CHROM,grp)][,grp:=NULL][]
}
union <- merge_bed(union)
fwrite(union,file.path(out,'regions.bed'),sep='\t',col.names=FALSE)
for(s in spec$samples) {
  m <- allmaps[[s$id]]; valid <- m[m$included,]
  b <- copy(union); j <- match(b$CHROM,valid$FINAL_CHR)
  excluded <- b[is.na(j)]
  excluded[,reason:=rep('TARGET_CHROMOSOME_EXCLUDED_OR_UNAVAILABLE',.N)]
  fwrite(excluded,paths(s$id,'.excluded_regions.tsv'),sep='\t')
  b <- b[!is.na(j)]; j <- match(b$CHROM,valid$FINAL_CHR)
  b[,CHROM:=valid$CHROM[j]]
  b <- clip_to_panel(b,read_panel(s$panel))
  j <- match(b$CHROM,valid$CHROM)
  fwrite(b,paths(s$id,'.regions.native.bed'),sep='\t',col.names=FALSE)
  b[,`:=`(start=pmax(0,start-spec$flank_bp),end=pmin(valid$expected_length[j],end+spec$flank_bp))]
  b <- clip_to_panel(b,read_panel(s$panel))
  fwrite(merge_bed(b),paths(s$id,'.expanded.native.bed'),sep='\t',col.names=FALSE)
  file.copy(s$mapping,paths(s$id,'.mapping.tsv'),overwrite=TRUE)
}
m <- rbindlist(allmaps)[included==TRUE,.(CHROM=FINAL_CHR,FINAL_CHR,ASSEMBLY,expected_length,included=TRUE)]
if(any(m[,uniqueN(expected_length),by=CHROM]$V1!=1)) stop('Incompatible chromosome lengths')
m <- unique(m,by='CHROM'); setorder(m,CHROM)
fwrite(m,file.path(out,'chromosomes.tsv'),sep='\t')
fwrite(rbindlist(manifest),file.path(out,'samples.tsv'),sep='\t')
jsonlite::write_json(spec,file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
