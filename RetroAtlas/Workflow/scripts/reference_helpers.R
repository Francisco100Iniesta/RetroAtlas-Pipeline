# Shared reference, coordinate and callability contract for reconstruction and reports.
read_chromosome_map <- function(path) {
  m <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
  m$included <- tolower(as.character(m$included)) %in% c('true','1','t')
  m
}
read_panel <- function(path=NULL) {
  if(is.null(path) || !length(path) || !nzchar(path)) return(NULL)
  if(!file.exists(path)) stop('Panel BED missing: ',path)
  if(file.info(path)$size == 0) return(data.table::data.table(CHROM=character(),start=numeric(),end=numeric()))
  data.table::fread(path,header=FALSE,select=1:3,col.names=c('CHROM','start','end'),
                    colClasses=list(character=1L,numeric=2:3))
}
clip_to_panel <- function(regions, panel) {
  d <- data.table::as.data.table(regions)
  if(is.null(panel)) return(data.table::copy(d))
  if(!nrow(d) || !nrow(panel)) return(d[0])
  x <- d[,.(CHROM,start=start+1,end)]
  y <- data.table::as.data.table(panel)[,.(CHROM,start=start+1,end)]
  data.table::setkey(y,CHROM,start,end)
  hits <- data.table::foverlaps(x,y,nomatch=0L,which=TRUE)
  z <- data.table::copy(d[hits$xid])
  z[,`:=`(start=pmax(d$start[hits$xid],y$start[hits$yid]-1),
           end=pmin(d$end[hits$xid],y$end[hits$yid]))]
  z[]
}
panel_membership <- function(loci, panel) {
  if(is.null(panel)) return(rep(TRUE,nrow(loci)))
  pts <- data.table::data.table(CHROM=as.character(loci$CHROM),start=as.numeric(loci$POS)-1,
                                end=as.numeric(loci$POS),panel_row=seq_len(nrow(loci)))
  hits <- clip_to_panel(pts,panel)
  seq_len(nrow(loci)) %in% hits$panel_row
}
standardize_loci <- function(loci, mapping, sample, sex, procedence, assembly, panel=NULL) {
  idx <- match(as.character(loci$CHROM), mapping$CHROM)
  loci$FINAL_CHR <- mapping$FINAL_CHR[idx]
  loci$ASSEMBLY <- rep(assembly, nrow(loci))
  loci$SAMPLE_ID <- rep(sample, nrow(loci))
  loci$SEX <- rep(sex, nrow(loci))
  loci$PROCEDENCE <- rep(procedence, nrow(loci))
  reason <- mapping$status[idx]
  reason[is.na(reason)] <- 'UNKNOWN_CONTIG'
  pos <- suppressWarnings(as.numeric(loci$POS))
  badpos <- !is.finite(pos) | pos < 1 | pos != floor(pos) |
    (!is.na(idx) & pos > mapping$expected_length[idx])
  reason[reason == 'PASS' & badpos] <- 'INVALID_POSITION'
  reason[sex == 'F' & !is.na(loci$FINAL_CHR) & loci$FINAL_CHR == 'chrY'] <- 'SEX_F_EXCLUDES_Y'
  eligible <- which(reason == 'PASS')
  reason[eligible[!panel_membership(loci[eligible,,drop=FALSE],panel)]] <- 'OUTSIDE_TARGET_PANEL'
  loci$EXCLUSION_REASON <- ifelse(reason == 'PASS', NA_character_, reason)
  list(kept=loci[is.na(loci$EXCLUSION_REASON),,drop=FALSE],
       excluded=loci[!is.na(loci$EXCLUSION_REASON),,drop=FALSE])
}
weighted_depth_quantiles <- function(depth, bp, probs=c(.25,.5,.75)) {
  if(!length(depth) || sum(bp) <= 0) return(rep(NA_real_,length(probs)))
  ix <- order(depth); v <- depth[ix]; w <- bp[ix]; cumulative <- cumsum(w)
  # Exact type-7 quantiles of per-base regional means without expanding bases.
  rank <- 1+(sum(w)-1)*probs
  at_rank <- function(r) v[vapply(r,function(k) which(cumulative >= k)[1L],integer(1))]
  lo <- floor(rank); hi <- ceiling(rank)
  at_rank(lo)+(rank-lo)*(at_rank(hi)-at_rank(lo))
}
read_coverage <- function(path, mapping, min_depth=0, upper_ratio=3, sex='U', panel=NULL) {
  # fread does not require an external decompressor when R.utils is available.
  d <- data.table::fread(path, header=FALSE, colClasses=c('character','numeric','numeric','numeric'), col.names=c('CHROM','start','end','depth'))
  if (ncol(d) != 4L || !nrow(d)) stop('Expected nonempty four-column regions BED')
  if (any(!is.finite(d$depth) | d$depth < 0 | d$start < 0 | d$end <= d$start)) stop('Invalid coverage interval/depth')
  idx <- match(d$CHROM, mapping$CHROM)
  d$FINAL_CHR <- mapping$FINAL_CHR[idx]
  d <- d[!is.na(idx) & mapping$included[idx]]
  if (!nrow(d)) stop('No coverage for validated primary chromosomes')
  data.table::setorder(d, CHROM, start, end)
  if (any(d[, start < data.table::shift(end, fill=-1), by=CHROM]$V1)) stop('Coverage intervals overlap')
  lim <- mapping$expected_length[match(d$CHROM, mapping$CHROM)]
  if (any(d$end > lim)) stop('Coverage interval outside reference length')
  autosomal <- grepl('^chr([1-9]|1[0-9]|2[0-2])$', d$FINAL_CHR)
  d <- clip_to_panel(d,panel)
  if(sex == 'F') d <- d[FINAL_CHR != 'chrY']
  autosomal <- grepl('^chr([1-9]|1[0-9]|2[0-2])$',d$FINAL_CHR)
  base_rows <- if(any(autosomal)) autosomal else d$FINAL_CHR != 'chrM'
  base_q <- weighted_depth_quantiles(d$depth[base_rows],(d$end-d$start)[base_rows])
  d[,`:=`(baseline_depth=base_q[2],depth_q1=base_q[1],depth_q3=base_q[3])]
  if(sex == 'M') for(ch in c('chrX','chrY')) {
    ix <- which(d$FINAL_CHR == ch)
    q <- weighted_depth_quantiles(d$depth[ix],(d$end-d$start)[ix])
    d[ix,`:=`(baseline_depth=q[2],depth_q1=q[1],depth_q3=q[3])]
  }
  d[,depth_ratio:=ifelse(is.finite(baseline_depth) & baseline_depth>0,depth/baseline_depth,NA_real_)]
  d[,cover_status:=ifelse(depth>0 & depth>=min_depth,'PASS','NO')]
  d[,region_status:=ifelse(is.finite(depth_ratio) & depth_ratio<=upper_ratio,'PASS','NO')]
  # Keep the normalization median unchanged; describe only callable depths.
  ok <- d$cover_status=='PASS' & d$region_status=='PASS'
  quartile_rows <- base_rows & ok
  if(!any(quartile_rows)) {
    quartile_rows <- ok & !(sex=='M' & d$FINAL_CHR %in% c('chrX','chrY'))
    if(any(quartile_rows & d$FINAL_CHR!='chrM')) quartile_rows <- quartile_rows & d$FINAL_CHR!='chrM'
  }
  callable_q <- weighted_depth_quantiles(d$depth[quartile_rows],(d$end-d$start)[quartile_rows])
  d[,`:=`(depth_q1=callable_q[1],depth_q3=callable_q[3])]
  if(sex == 'M') for(ch in c('chrX','chrY')) {
    ix <- which(d$FINAL_CHR == ch & ok)
    q <- weighted_depth_quantiles(d$depth[ix],(d$end-d$start)[ix])
    d[FINAL_CHR == ch,`:=`(depth_q1=q[1],depth_q3=q[3])]
  }
  d[,depth_band:=ifelse(!ok,'NotCallable',ifelse(!is.finite(depth_q1),'Unavailable',
                        ifelse(depth<depth_q1,'LowerIQR',ifelse(depth>depth_q3,'UpperIQR','IQR'))))]
  # Average RGB channels, not category scores: mixed blue/red must not become green.
  d[,`:=`(color_r=255,color_g=255,color_b=255)]
  d[cover_status=='PASS' & !is.finite(depth_ratio),`:=`(color_r=217,color_g=223,color_b=227)]
  d[is.finite(depth_ratio) & depth_ratio>upper_ratio,`:=`(color_r=0,color_g=0,color_b=0)]
  lower <- which(d$depth_band=='LowerIQR')
  upper <- which(d$depth_band=='UpperIQR')
  within <- which(d$depth_band=='IQR')
  d[lower,`:=`(color_r=0,color_g=128,color_b=255)]
  d[within,`:=`(color_r=0,color_g=255,color_b=128)]
  d[upper,`:=`(color_r=255,color_g=0,color_b=0)]
  list(regions=d, baseline=base_q[2], min_depth=min_depth, upper_ratio=upper_ratio,sex=sex,
       panel=!is.null(panel),baseline_source=if(any(autosomal)) 'autosomal' else 'available non-mitochondrial target',
       baselines=unique(d[,.(FINAL_CHR,baseline_depth,depth_q1,depth_q3)]))
}
annotate_coverage <- function(loci, coverage) {
  d <- data.table::as.data.table(loci)
  d[, variant_id := .I]
  # BED [start,end) corresponds to 1-based positions start < POS <= end.
  pts <- d[, .(variant_id, CHROM, start=as.numeric(POS), end=as.numeric(POS))]
  hits <- coverage$regions[pts, on=.(CHROM, start < start, end >= end), nomatch=0,
                           .(variant_id=i.variant_id, depth, cover_status, region_status)]
  if (anyDuplicated(hits$variant_id)) stop('Ambiguous overlapping coverage windows')
  idx <- match(d$variant_id, hits$variant_id)
  d$depth <- hits$depth[idx]
  d$cover_status <- hits$cover_status[idx]
  d$region_status <- hits$region_status[idx]
  d$EXCLUSION_REASON <- ifelse(is.na(idx), 'NO_COVERAGE_WINDOW',
    ifelse(d$cover_status != 'PASS', 'LOW_DEPTH', ifelse(d$region_status != 'PASS', 'HIGH_DEPTH_OR_ZERO_BASELINE', NA_character_)))
  d[, variant_id := NULL]
  as.data.frame(d)
}
coverage_display <- function(coverage, bin_bp=1000000) {
  d <- data.table::copy(coverage$regions)
  d[,`:=`(ok=cover_status=='PASS' & region_status=='PASS',bin=floor(start/bin_bp))]
  crosses <- d$bin != floor((d$end-1)/bin_bp)
  simple <- d[!crosses,.(chr=FINAL_CHR,bin,bp=end-start,ok,depth_band,color_r,color_g,color_b)]
  pieces <- lapply(which(crosses),function(i) {
    bins <- seq(floor(d$start[i]/bin_bp),floor((d$end[i]-1)/bin_bp))
    bp <- pmin(d$end[i],(bins+1)*bin_bp)-pmax(d$start[i],bins*bin_bp)
    data.table::data.table(chr=d$FINAL_CHR[i],bin=bins,bp=bp,ok=d$ok[i],
                           depth_band=d$depth_band[i],color_r=d$color_r[i],color_g=d$color_g[i],color_b=d$color_b[i])
  })
  z <- data.table::rbindlist(c(list(simple),pieces))
  bins <- z[,.(bp=sum(bp),callableBp=sum(bp[ok]),lowerBp=sum(bp[depth_band=='LowerIQR']),
               iqrBp=sum(bp[depth_band=='IQR']),upperBp=sum(bp[depth_band=='UpperIQR']),
               colorR=sum(bp*color_r)/sum(bp),colorG=sum(bp*color_g)/sum(bp),colorB=sum(bp*color_b)/sum(bp)),by=.(chr,bin)]
  bins[,pct:=100*callableBp/bp]
  info <- bins[,.(totalBp=sum(bp),callableBp=sum(callableBp),lowerBp=sum(lowerBp),
                  iqrBp=sum(iqrBp),upperBp=sum(upperBp)),by=chr]
  info[,`:=`(pct=100*callableBp/totalBp,lowerPct=ifelse(callableBp>0,100*lowerBp/callableBp,0),
              iqrPct=ifelse(callableBp>0,100*iqrBp/callableBp,0),upperPct=ifelse(callableBp>0,100*upperBp/callableBp,0))]
  total <- sum(bins$bp)
  callable_total <- sum(bins$callableBp)
  callable_percent <- function(v) if(callable_total>0) 100*sum(v)/callable_total else 0
  percent <- function(v) if(total>0) 100*sum(v)/total else 0
  list(baselineDepth=coverage$baseline,minDepth=coverage$min_depth,
       thresholdDepthRatio=coverage$upper_ratio,sex=coverage$sex,panel=coverage$panel,
       baselineSource=coverage$baseline_source,baselines=as.data.frame(coverage$baselines),
       binBp=bin_bp,totalBp=total,callableBp=sum(bins$callableBp),
       pctCallable=percent(bins$callableBp),lowerPct=callable_percent(bins$lowerBp),
       iqrPct=callable_percent(bins$iqrBp),upperPct=callable_percent(bins$upperBp),regionCount=nrow(d),
       bins=as.data.frame(bins),chromosomes=as.data.frame(info))
}
write_flat_tsv <- function(x, path) {
  # Preserve list-valued audit columns as JSON without losing rows.
  x <- as.data.frame(x)
  for (k in names(x)) if(is.list(x[[k]])) x[[k]] <- vapply(x[[k]], function(v) as.character(jsonlite::toJSON(v, auto_unbox=TRUE, na='null')), character(1))
  data.table::fwrite(x, path, sep='\t', quote='auto', na='NA')
}

# ALT read counts are never summed across callers or reconstructed components.
primary_alt_support <- function(d) {
  n <- nrow(d)
  number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  result <- rep(NA_real_,n)
  for(k in intersect(c('ALT_SUPPORT_MAX','PRIMARY_ALT_SUPPORT','alt_support'),names(d))) {
    v <- number(d[[k]])
    result <- ifelse(is.finite(v),ifelse(is.finite(result),pmax(result,v),v),result)
  }
  if('LR' %in% names(d)) {
    v <- vapply(strsplit(as.character(ifelse(is.na(d$LR),'',as.character(d$LR))),',',fixed=TRUE),
       function(z) if(length(z)) number(z[if(length(z)>1L) 2L else 1L]) else NA_real_,numeric(1))
    result <- ifelse(is.finite(v),ifelse(is.finite(result),pmax(result,v),v),result)
  }
  result
}

primary_partition <- function(loci, coverage=NULL, min_support=3, min_svlen=30) {
  if(!is.finite(min_support)||min_support<1||!is.finite(min_svlen)||min_svlen<1)
    stop('Primary support and size thresholds must be positive')
  d <- as.data.frame(loci)
  if(!is.null(coverage)) d <- annotate_coverage(d,coverage)
  for(k in c('cover_status','region_status')) if(!k %in% names(d)) d[[k]] <- rep('NO',nrow(d))
  d$PRIMARY_ALT_SUPPORT <- primary_alt_support(d)
  callable <- d$cover_status %in% 'PASS' & d$region_status %in% 'PASS'
  len <- if('svlen' %in% names(d)) d$svlen else d$SVLEN
  if(is.null(len)) len <- rep(NA_real_,nrow(d))
  size_ok <- is.finite(as.numeric(len)) & abs(as.numeric(len))>=min_svlen
  support_ok <- is.finite(d$PRIMARY_ALT_SUPPORT) & d$PRIMARY_ALT_SUPPORT>=min_support
  d$PRIMARY_FILTER_REASON <- ifelse(!callable,'NOT_CALLABLE',ifelse(!support_ok,'LOW_OR_MISSING_ALT_SUPPORT',ifelse(!size_ok,'BELOW_MIN_SVLEN',NA_character_)))
  d$PRIMARY_REPORTABLE <- is.na(d$PRIMARY_FILTER_REASON)
  list(loci=d[d$PRIMARY_REPORTABLE,,drop=FALSE],
       rescue_candidates=d[!d$PRIMARY_REPORTABLE,,drop=FALSE],
       nocallablecandidates=d[!callable,,drop=FALSE])
}

# UCSC MANE (BED12+13) and RefSeq Select Curated (bigGenePred BED12+8).
# Keep the source status explicit: RefSeq models are not MANE transcripts.
read_gene_annotation <- function(path) {
  b <- utils::read.table(path,sep='\t',header=FALSE,quote='',comment.char='',
                         stringsAsFactors=FALSE,colClasses='character',fill=TRUE)
  if(ncol(b)==20L) {
    b[[19L]] <- ifelse(nzchar(b[[19L]]),b[[19L]],b[[13L]])
    b[[21L]] <- ''; b[[22L]] <- b[[4L]]; b[[23L]] <- ''; b[[24L]] <- ''
    b[[25L]] <- 'RefSeq Select Curated'
  } else if(ncol(b)!=25L) stop('Gene annotation requires MANE BED25 or RefSeq bigGenePred BED20')
  b
}
