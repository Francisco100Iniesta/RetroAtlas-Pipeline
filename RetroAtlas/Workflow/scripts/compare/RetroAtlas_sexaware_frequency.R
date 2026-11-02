# RetroAtlas sex-aware ploidy and AF. Functions only, requires data.table.
# Sex karyotype metadata are supplied explicitly; NO sex inference from IDs.
# Assumes reference GRCh38 or GRCh37, 1-based inclusive coordinates.
# PAR representation: X PAR counted; Y PAR quarantined pending X/Y normalization.

sa_chr <- function(z) toupper(sub('^CHR','',as.character(z),ignore.case=TRUE))
sa_build <- function(build) {
  b <- toupper(as.character(build))
  if(b %in% c('GRCH38','GRCH38.P14','HG38')) return('GRCh38')
  if(b %in% c('T2T','T2T-CHM13V2.0')) return('T2T-CHM13v2.0')
  if(b %in% c('GRCH37','HG19')) return('GRCh37')
  stop('Sex-aware frequency needs --reference-build GRCh38/hg38 or GRCh37/hg19; got ',build)
}
sa_ranges <- function(build) {
  b <- sa_build(build)
  if(b=='T2T-CHM13v2.0') {
    # Official CHM13 analysis-set PAR BED (0-based half-open), checked in.
    bed <- read.delim(file.path(compare_script_dir,'../../data/reference/T2T.PAR.bed'),header=FALSE)
    return(lapply(setNames(c('chrX','chrY'),c('X','Y')),function(ch) {
      q <- as.matrix(bed[bed[[1]]==ch,2:3]); q[,1] <- q[,1]+1; q
    }))
  }
  if(b=='GRCh38') {
    list(X=matrix(c(10001,2781479,155701383,156030895),ncol=2,byrow=TRUE),
         Y=matrix(c(10001,2781479,56887903,57217415),ncol=2,byrow=TRUE))
  } else {
    list(X=matrix(c(60001,2699520,154931044,155260560),ncol=2,byrow=TRUE),
         Y=matrix(c(10001,2649520,59034050,59363566),ncol=2,byrow=TRUE))
  }
}
sa_region <- function(chr,pos,build) {
  x <- sa_chr(chr); p <- suppressWarnings(as.numeric(pos))
  r <- sa_ranges(build); z <- rep('UNSUPPORTED',length(x))
  z[x %in% as.character(1:22)] <- 'AUTOSOME'
  z[x=='MT' | x=='M'] <- 'MT_UNMODELED'
  for (ch in c('X','Y')) {
    sel <- x==ch & is.finite(p)
    z[sel] <- paste0(ch,'_NONPAR')
    q <- r[[ch]]
    inpar <- sel & ((p>=q[1,1] & p<=q[1,2]) | (p>=q[2,1] & p<=q[2,2]))
    z[inpar] <- paste0(ch,'_PAR')
  }
  z[x %in% c('X','Y') & !is.finite(p)] <- 'UNSUPPORTED'
  z
}
sa_manifest <- function(manifest,enabled=TRUE) {
  m <- data.table::copy(data.table::as.data.table(manifest))
  if(!'sample_id' %in% names(m) || anyDuplicated(m$sample_id))
    stop('Invalid or duplicated sample_id in manifest')
  # 'sex' is a reported metadata code, NOT a measured chromosomal karyotype.
  # For standard XX/XY baseline ploidy only, F->XX and M->XY are ASSUMPTIONS.
  # When provided, sex_karyotype (XX/XY/UNKNOWN/complex) overrides this proxy.
  if(!'sex' %in% names(m)) m[,sex:='U']
  m[,sex:=toupper(trimws(as.character(sex)))]
  m[is.na(sex) | !nzchar(sex),sex:='U']
  m[sex=='FEMALE',sex:='F']
  m[sex=='MALE',sex:='M']
  m[sex=='UNKNOW',sex:='U']
  bad <- setdiff(unique(m$sex),c('F','M','U'))
  if(length(bad)) stop('Invalid sex code: ',paste(bad,collapse=','),'. Allowed F/M/U.')
  if(!'sex_karyotype' %in% names(m)) m[,sex_karyotype:=NA_character_]
  m[,sex_karyotype:=toupper(trimws(as.character(sex_karyotype)))]
  m[is.na(sex_karyotype) | sex_karyotype %in% c('','NA','U','UNK'),sex_karyotype:='UNKNOWN']
  explicit <- m$sex_karyotype!='UNKNOWN'
  valid <- c('XX','XY','UNKNOWN','XXY','X0','XO','XYY','XXX')
  badk <- setdiff(unique(m$sex_karyotype),valid)
  if(length(badk)) stop('Invalid sex_karyotype values: ',paste(badk,collapse=','))
  m[,sex_karyotype_source:=ifelse(explicit,'EXPLICIT','UNDETERMINED')]
  m[sex_karyotype=='UNKNOWN' & sex=='F',`:=`(
     sex_karyotype='XX',sex_karyotype_source='ASSUMED_FROM_F')]
  m[sex_karyotype=='UNKNOWN' & sex=='M',`:=`(
     sex_karyotype='XY',sex_karyotype_source='ASSUMED_FROM_M')]
  # UNKNOWN/U => no X/Y frequency contributions, even if positive evidence.
  # PROCEDENCE is the single source of origin, including older manifests.
  if('procedence' %in% names(m)) m[,origin:=as.character(procedence)]
  if(!'origin' %in% names(m)) m[,origin:='UNKNOW']
  m[,origin:=toupper(trimws(as.character(origin)))]
  m[is.na(origin) | !nzchar(origin) | origin=='U',origin:='UNKNOW']
  obsolete <- intersect(c('ancestry_group','ethnicity'),names(m))
  if(length(obsolete)) m[,(obsolete):=NULL]
  m[]
}

sa_ploidy <- function(region,karyo) {
  reg <- as.character(region); sx <- as.character(karyo)
  p <- rep(NA_integer_,length(reg))
  p[reg=='AUTOSOME'] <- 2L
  standard <- sx %in% c('XX','XY')
  p[reg=='X_PAR' & standard] <- 2L
  p[reg=='X_NONPAR' & sx=='XX'] <- 2L
  p[reg=='X_NONPAR' & sx=='XY'] <- 1L
  p[reg=='Y_NONPAR' & sx=='XX'] <- 0L
  p[reg=='Y_NONPAR' & sx=='XY'] <- 1L
  # Y PAR is a duplicate representation. Do not count unless mapped to X PAR.
  p
}
sa_gt <- function(z) {
  g <- as.character(z)
  g <- sub(':.*$','',g)
  g <- gsub('|','/',g,fixed=TRUE)
  g[is.na(g) | !nzchar(g) | g %in% c('.','./.')] <- NA_character_
  g[g=='1/0'] <- '0/1'
  g
}
sa_adjust_states <- function(states,manifest,primary_members=NULL,reference_build='GRCh38') {
  d <- data.table::copy(data.table::as.data.table(states))
  m <- sa_manifest(manifest,enabled=TRUE)
  mm <- m[,.(sample_id,sex,sex_karyotype,sex_karyotype_source,origin)]
  if(anyDuplicated(mm$sample_id)) stop('Duplicate sample metadata')
  d <- merge(d,mm,by='sample_id',all.x=TRUE,sort=FALSE)
  if(anyNA(d$sex_karyotype)) stop('Missing sex_karyotype for one or more samples')
  d[,`:=`(sa_region=sa_region(CHROM,POS,reference_build),
         sa_gt_original=as.character(genotype_final),
         sa_state_original=as.character(pair_state),
         sa_af_state_original=as.character(af_state))]
  d[,sa_ploidy:=sa_ploidy(sa_region,sex_karyotype)]
  d[,`:=`(sa_qc='OK',sa_raw_primary_gt=NA_character_)]
  if(!is.null(primary_members) && nrow(primary_members)) {
    p <- data.table::as.data.table(primary_members)
    if(all(c('locus_id','sample_id','GT') %in% names(p))) {
      p <- unique(p[,.(locus_id,sample_id,sa_raw_primary_gt=as.character(GT))],by=c('locus_id','sample_id'))
      d[p,on=.(locus_id,sample_id),sa_raw_primary_gt:=i.sa_raw_primary_gt]
    }
  }
  rawgt <- as.character(d$genotype_final)
  prim <- d$af_state=='PRIMARY_GENOTYPED' & !is.na(d$sa_raw_primary_gt) & nzchar(d$sa_raw_primary_gt)
  rawgt[prim] <- d$sa_raw_primary_gt[prim]
  rawgt <- sa_gt(rawgt)
  d[,sa_input_gt:=rawgt]
  # Only observed genotype states can be normalized as GT; softclip candidates
  # retain GT=NA and their separately labeled assumed one-copy model.
  observed <- which(d$af_state %in% c('PRIMARY_GENOTYPED','RESCUE_GENOTYPED'))
  if(length(observed)) {
    pl <- d$sa_ploidy[observed]; g <- d$sa_input_gt[observed]
    dose <- rep(NA_integer_,length(observed))
    norm <- rep(NA_character_,length(observed))
    dip <- pl==2L & !is.na(pl)
    hap <- pl==1L & !is.na(pl)
    dose[dip & g=='0/1'] <- 1L;norm[dip & g=='0/1'] <- '0/1'
    dose[dip & g=='1/1'] <- 2L;norm[dip & g=='1/1'] <- '1/1'
    # Haploid 1/1 is frequently a legacy diploid-caller representation.
    # Explicitly normalize, and label it rather than silently awarding 2 ALT.
    dose[hap & g %in% c('1','1/1')] <- 1L
    norm[hap & g %in% c('1','1/1')] <- '1'
    good <- is.finite(dose)
    if(any(good)) {
      k <- observed[good]
      d[k,`:=`(alt_copies_observed=dose[good],genotype_final=norm[good],GT_final=norm[good])]
      conv <- k[g[good]=='1/1' & pl[good]==1L]
      if(length(conv)) d[conv,sa_qc:='HAPLOID_DIPLOID_1_1_NORMALIZED']
    }
    if(any(!good)) {
      k <- observed[!good]
      primbad <- k[d$af_state[k]=='PRIMARY_GENOTYPED']
      resbad <- setdiff(k,primbad)
      if(length(primbad)) d[primbad,`:=`(
        af_state='PRIMARY_PLOIDY_GT_REVIEW',sa_qc='PLOIDY_GT_CONFLICT',
        alt_copies_observed=NA_integer_,genotype_final=NA_character_,GT_final=NA_character_)]
      if(length(resbad)) d[resbad,`:=`(
        af_state='RESCUE_PLOIDY_GT_REVIEW',pair_state='REVIEW',sa_qc='PLOIDY_GT_CONFLICT',
        alt_copies_observed=NA_integer_,genotype_final=NA_character_,GT_final=NA_character_)]
    }
  }
  # Zero-copy XX-Y observations are NOT reference genotypes; preserve all raw
  # calls in sa_state_original and v7 evidence, but quarantine positive matches.
  zero <- which(!is.na(d$sa_ploidy) & d$sa_ploidy==0L)
  if(length(zero)) {
    positive <- zero[d$sa_state_original[zero] %in% c('PRIMARY','RESCUE')]
    d[zero,`:=`(sa_qc='ZERO_COPY_Y_EXCLUDED',af_state='SEX_PLOIDY0_EXCLUDED',
        pair_state='NOT_COMPARABLE',genotype_final=NA_character_,GT_final=NA_character_,
        alt_copies_observed=NA_integer_,alt_copies_assumed=NA_integer_,
        alt_copies_adjusted=NA_integer_)]
    if(length(positive)) d[positive,sa_qc:='POSSIBLE_SPURIOUS_Y_XX_EVIDENCE']
  }
  unmodeled <- which(is.na(d$sa_ploidy))
  if(length(unmodeled)) {
    status <- ifelse(d$sa_region[unmodeled]=='Y_PAR','Y_PAR_DUPLICATE_UNNORMALIZED',
      ifelse(d$sa_region[unmodeled]=='MT_UNMODELED','MT_PLOIDY_UNMODELED',
      ifelse(d$sex_karyotype[unmodeled]=='UNKNOWN','SEX_KARYOTYPE_UNKNOWN','SEX_PLOIDY_UNMODELED')))
    d[unmodeled,`:=`(sa_qc=status,af_state='SEX_PLOIDY_UNMODELED',
         # Observed PRIMARY remains a PRIMARY; no GT/AF is assigned.
         # Non-primary uncertain states do not count as verified reference.
         pair_state=ifelse(sa_state_original=='PRIMARY','PRIMARY',
                           ifelse(sa_state_original=='RESCUE','REVIEW','NOT_COMPARABLE')),
         genotype_final=NA_character_,GT_final=NA_character_,
         alt_copies_observed=NA_integer_,alt_copies_assumed=NA_integer_,
         alt_copies_adjusted=NA_integer_)]
  }
  # Softclip remains a one-copy assumption even in diploid contexts.
  d[,is_assumed:=af_state=='RESCUE_DOSAGE_ASSUMED' & !is.na(sa_ploidy) & sa_ploidy>0L]
  d[,AF_imputation_flag:=is_assumed]
  d[,is_callable_reference:=af_state=='REF_CALLABLE_CANDIDATE' & !is.na(sa_ploidy) & sa_ploidy>0L]
  d[,is_af_gt_evaluable:=af_state %in% c('PRIMARY_GENOTYPED','RESCUE_GENOTYPED','REF_CALLABLE_CANDIDATE') &
       !is.na(sa_ploidy) & sa_ploidy>0L]
  d[,is_af_adjusted_evaluable:=is_af_gt_evaluable | is_assumed]
  d[,ALT_dosage_observed:=alt_copies_observed]
  d[,ALT_dosage_assumed:=alt_copies_assumed]
  d[,alt_copies_adjusted:=as.integer(ifelse(is_assumed,alt_copies_assumed,alt_copies_observed))]
  if(any(d$is_af_gt_evaluable & d$alt_copies_observed>d$sa_ploidy,na.rm=TRUE))
    stop('Sex-aware AF invariant violated: ALT dose exceeds local ploidy')
  if(anyDuplicated(d[,.(locus_id,sample_id)])) stop('Duplicated locus/sample after sex-aware metadata join')
  d[]
}
sa_frequency <- function(sample_states,sample_ids=NULL,autosomes_only=FALSE,
                         include_assumed=TRUE,required_ref='callable') {
  if(!identical(required_ref,'callable')) stop('Only callable-reference model implemented')
  d <- data.table::copy(data.table::as.data.table(sample_states))
  if(!is.null(sample_ids)) d <- d[sample_id %in% sample_ids]
  if(!nrow(d)) stop('No eligible samples for sex-aware frequency')
  needed <- c('sa_ploidy','sa_region','is_af_gt_evaluable','is_assumed','alt_copies_observed')
  if(length(setdiff(needed,names(d)))) stop('States not sex-aware: missing ',paste(setdiff(needed,names(d)),collapse=','))
  if(anyDuplicated(d[,.(locus_id,sample_id)])) stop('Duplicate locus/sample in frequency input')
  d[,af_eligible_chrom:=!is.na(sa_ploidy) & sa_ploidy>0L &
      (!autosomes_only | sa_region=='AUTOSOME')]
  d[,gt_evaluable:=is_af_gt_evaluable & af_eligible_chrom]
  d[,adjusted_evaluable:=(is_af_gt_evaluable | (include_assumed & is_assumed)) & af_eligible_chrom]
  d[,`:=`(sa_allele_n_gt=as.integer(ifelse(gt_evaluable,sa_ploidy,0L)),
          sa_allele_n_adj=as.integer(ifelse(adjusted_evaluable,sa_ploidy,0L)),
          sa_allele_c_gt=as.integer(ifelse(gt_evaluable,alt_copies_observed,0L)),
          sa_allele_c_adj=as.integer(ifelse(adjusted_evaluable,alt_copies_adjusted,0L)))]
  z <- d[,.(CHROM=CHROM[1L],POS=POS[1L],SVLEN=SVLEN[1L],sa_region=sa_region[1L],
    N_SAMPLES=.N,N_PRIMARY=sum(pair_state=='PRIMARY'),
    N_PRIMARY_GENOTYPED=sum(af_state=='PRIMARY_GENOTYPED'),
    N_PRIMARY_UNTYPED=sum(af_state %in% c('PRIMARY_UNTYPED','PRIMARY_PLOIDY_GT_REVIEW') |
       (pair_state=='PRIMARY' & af_state=='SEX_PLOIDY_UNMODELED')),
    N_RESCUE_GENOTYPED=sum(af_state=='RESCUE_GENOTYPED'),
    N_RESCUE_ASSUMED=sum(is_assumed),N_REVIEW=sum(pair_state=='REVIEW'),
    N_NOT_COMPARABLE=sum(pair_state=='NOT_COMPARABLE'),
    N_REF_CALLABLE=sum(is_callable_reference),
    N_PLOIDY_GT_CONFLICT=sum(sa_qc=='PLOIDY_GT_CONFLICT'),
    N_Y_XX_SPURIOUS=sum(sa_qc=='POSSIBLE_SPURIOUS_Y_XX_EVIDENCE'),
    N_PLOIDY_UNKNOWN=sum(is.na(sa_ploidy)),
    AC_GT=sum(sa_allele_c_gt,na.rm=TRUE),AN_GT=sum(sa_allele_n_gt,na.rm=TRUE),
    AC_ADJUSTED=sum(sa_allele_c_adj,na.rm=TRUE),AN_ADJUSTED=sum(sa_allele_n_adj,na.rm=TRUE),
    N_AF_GT=sum(gt_evaluable),N_AF_ADJUSTED=sum(adjusted_evaluable)),by=locus_id]
  z[,AF_GT:=data.table::fifelse(AN_GT>0L,AC_GT/AN_GT,NA_real_)]
  z[,AF_RESCUE_ADJUSTED:=data.table::fifelse(AN_ADJUSTED>0L,AC_ADJUSTED/AN_ADJUSTED,NA_real_)]
  z[,MAF_GT:=data.table::fifelse(is.finite(AF_GT),pmin(AF_GT,1-AF_GT),NA_real_)]
  z[,MAF_RESCUE_ADJUSTED:=data.table::fifelse(is.finite(AF_RESCUE_ADJUSTED),
      pmin(AF_RESCUE_ADJUSTED,1-AF_RESCUE_ADJUSTED),NA_real_)]
  z[,AF_STATUS:=data.table::fifelse(autosomes_only & sa_region!='AUTOSOME',
       'NON_AUTOSOMAL_DISABLED',data.table::fifelse(sa_region=='Y_PAR',
       'Y_PAR_NEEDS_X_NORMALIZATION',data.table::fifelse(sa_region=='MT_UNMODELED',
       'MT_PLOIDY_UNMODELED',data.table::fifelse(AN_ADJUSTED==0L,
       'NO_EVALUABLE_COPIES','SEX_AWARE_CALLABLE_REF_CANDIDATES'))))]
  if(any(z$AC_GT>z$AN_GT | z$AC_ADJUSTED>z$AN_ADJUSTED,na.rm=TRUE))
    stop('AC>AN after sex-aware model')
  data.table::setorder(z,locus_id)
  z[]
}
sa_frequency_strata <- function(states,include_assumed=TRUE,autosomes_only=FALSE,
                                stratify='sex_origin',min_public_n=5L) {
  d <- data.table::as.data.table(states)
  mode <- tolower(as.character(stratify))
  if(!mode %in% c('none','sex','origin','sex_origin')) stop('Invalid --af-stratify')
  if(mode=='none') return(data.table::data.table())
  cols <- switch(mode, sex=list('sex'), origin=list('origin'),
     sex_origin=list('sex','origin',c('sex','origin')))
  res <- list();k<-0L
  for(nms in cols) {
    groups <- unique(d[,..nms])
    for(i in seq_len(nrow(groups))) {
      mask <- rep(TRUE,nrow(d))
      for(nm in nms) mask <- mask & d[[nm]]==groups[[nm]][i]
      sub <- d[mask]
      if(!nrow(sub)) next
      label <- paste(paste(nms,vapply(nms,function(nm)as.character(groups[[nm]][i]),character(1)),sep='='),collapse=';')
      k<-k+1L
      ff <- sa_frequency(sub,autosomes_only=autosomes_only,include_assumed=include_assumed)
      ff[,`:=`(STRATUM_TYPE=paste(nms,collapse='+'),STRATUM=label,
               N_STRATUM_SAMPLES=data.table::uniqueN(sub$sample_id),
               PUBLIC_ELIGIBLE=data.table::uniqueN(sub$sample_id)>=min_public_n)]
      res[[k]]<-ff
    }
  }
  data.table::rbindlist(res,fill=TRUE)
}
