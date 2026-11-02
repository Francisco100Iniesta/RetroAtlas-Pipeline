# Cohort output helpers. BED coordinates are always zero-based, half-open.
suppressPackageStartupMessages(library(data.table))

aw_pack <- function(d,cols=NULL) {
  d<-as.data.table(d)
  if(!is.null(cols)) d<-d[,intersect(cols,names(d)),with=FALSE]
  # Iterate column vectors directly: avoid a data.table row allocation per cell row.
  rows<-if(!nrow(d) || !ncol(d)) list() else
    unname(do.call(Map,c(list(f=function(...) unname(list(...))),unname(as.list(d)))))
  list(columns=names(d),rows=rows)
}

aw_pieces <- function(d, width) {
  if(!nrow(d)) return(data.table(bin=integer(),observed_bp=double(),callable_bp=double()))
  first<-floor(d$start/width)+1L;last<-floor((d$end-1)/width)+1L
  n<-last-first+1L;idx<-rep(seq_len(nrow(d)),n)
  bins<-rep(first,n)+sequence(n)-1L
  bp<-pmin(d$end[idx],bins*width)-pmax(d$start[idx],(bins-1L)*width)
  data.table(bin=as.integer(bins),observed_bp=bp,callable_bp=bp*(d$ok[idx] %in% TRUE))[
    ,.(observed_bp=sum(observed_bp),callable_bp=sum(callable_bp)),by=bin]
}

# Complete linkage on size, within a fixed genomic window. No chaining.
aw_size_groups <- function(size, similarity, tolerance) {
  if(!length(size)) return(integer())
  idx<-order(size);out<-integer(length(size));g<-0L;minimum<-NA_real_
  for(i in idx) {
    if(g==0L || size[[i]]-minimum>tolerance || minimum/size[[i]]<similarity) {
      g<-g+1L;minimum<-size[[i]]
    }
    out[[i]]<-g
  }
  out
}

aw_summarize_ins <- function(records, settings, cores=1L) {
  if(nrow(records) && cores>1L && uniqueN(records$CHROM)>1L) {
    parts<-aw_parallel(split(records,records$CHROM),function(r) aw_summarize_ins(r,settings),cores)
    return(list(summary=rbindlist(lapply(parts,`[[`,'summary')),
                records=rbindlist(lapply(parts,`[[`,'records'),fill=TRUE)))
  }
  empty<-data.table(CHROM=character(),start=double(),end=double(),INS_GROUP=character(),
    SVLEN_MEDIAN=double(),SVLEN_MIN=double(),SVLEN_MAX=double(),N_SAMPLES=integer(),
    N_CALLABLE_SAMPLES=integer(),N_NOT_CALLABLE_SAMPLES=integer(),N_RECORDS=integer(),
    N_ALT_GT_SAMPLES=integer(),N_UNKNOWN_GT_SAMPLES=integer(),N_REF_GT_SAMPLES=integer(),
    N_CUTESV=integer(),N_SNIFFLES2=integer(),N_SEVERUS=integer(),N_ORIGINAL=integer(),N_RESCUE=integer(),
    N_HET=integer(),N_HOM_ALT=integer(),N_HAP_ALT=integer(),N_GT_UNKNOWN=integer(),N_GT_CONFLICT=integer(),
    AC_OBSERVED=integer(),AN_OBSERVED=integer(),
    N_REVIEW_SAMPLES=integer(),N_NOT_EVALUATED_SAMPLES=integer(),
    N_NEAR_RESOLVED_SAMPLES=integer(),N_OTHER_EVALUATED_SAMPLES=integer())
  if(!nrow(records)) return(list(summary=empty,records=records))
  if(!'candidate_status' %in% names(records)) records<-copy(records)[,candidate_status:='NOT_EVALUATED']
  records<-copy(records)
  records[,INS_GROUP:=NA_character_]
  excluded<-copy(records[window_callable==FALSE])
  r<-copy(records[window_callable==TRUE])
  if(!nrow(r)) return(list(summary=empty,records=records))
  r[,group:=aw_size_groups(SVLEN,settings$min_size_similarity,settings$max_size_difference_bp),by=.(CHROM,start,end)]
  r[,INS_GROUP:=paste(CHROM,start,end,sprintf('INS%03d',group),sep=':')]
  z<-r[,{
    # One size vote per sample, avoiding caller multiplicity in the median.
    votes<-.SD[,.(size=median(SVLEN)),by=sample_id]
    genotypes<-.SD[,{
      good<-is.finite(dosage) & is.finite(ploidy)
      calls<-unique(paste(dosage[good],ploidy[good],sep='/'))
      conflict<-length(calls)>1L || any(genotype_qc %in% c('GT_PLOIDY_CONFLICT','PLOIDY_UNMODELED'))
      list(dose=if(length(calls)==1L && !conflict) dosage[which(good)[1L]] else NA_integer_,
        pl=if(length(calls)==1L && !conflict) ploidy[which(good)[1L]] else NA_integer_,
        conflict=conflict,review=any(candidate_status=='REVIEW'))
    },by=sample_id]
    list(SVLEN_MEDIAN=median(votes$size),SVLEN_MIN=min(SVLEN),SVLEN_MAX=max(SVLEN),
      N_SAMPLES=uniqueN(sample_id),N_CALLABLE_SAMPLES=uniqueN(sample_id[window_callable]),
      N_NOT_CALLABLE_SAMPLES=uniqueN(sample_id[!window_callable]),
      N_RECORDS=uniqueN(record_key),N_ALT_GT_SAMPLES=uniqueN(sample_id[gt_status=='ALT']),
      N_UNKNOWN_GT_SAMPLES=uniqueN(sample_id[gt_status=='UNKNOWN']),N_REF_GT_SAMPLES=uniqueN(sample_id[gt_status=='REF']),
      N_CUTESV=uniqueN(sample_id[caller=='cutesv']),N_SNIFFLES2=uniqueN(sample_id[caller=='sniffles2']),
      N_SEVERUS=uniqueN(sample_id[caller=='severus']),N_ORIGINAL=uniqueN(sample_id[source=='original']),
      N_RESCUE=uniqueN(sample_id[source=='rescue']),
      N_REVIEW_SAMPLES=uniqueN(sample_id[candidate_status=='REVIEW']),
      N_NOT_EVALUATED_SAMPLES=uniqueN(sample_id[candidate_status=='NOT_EVALUATED']),
      N_NEAR_RESOLVED_SAMPLES=uniqueN(sample_id[candidate_status=='NEAR_RESOLVED_LOCUS']),
      N_OTHER_EVALUATED_SAMPLES=uniqueN(sample_id[candidate_status=='EVALUATED_OTHER']),
      N_HET=sum(genotypes$pl==2 & genotypes$dose==1,na.rm=TRUE),
      N_HOM_ALT=sum(genotypes$pl==2 & genotypes$dose==2,na.rm=TRUE),
      N_HAP_ALT=sum(genotypes$pl==1 & genotypes$dose==1,na.rm=TRUE),
      N_GT_UNKNOWN=sum(is.na(genotypes$dose) & !genotypes$conflict),N_GT_CONFLICT=sum(genotypes$conflict),
      AC_OBSERVED=sum(genotypes$dose[!genotypes$review],na.rm=TRUE),
      AN_OBSERVED=sum(genotypes$pl[!genotypes$review],na.rm=TRUE))
  },by=.(CHROM,start,end,INS_GROUP)]
  list(summary=z,records=rbindlist(list(r,excluded),fill=TRUE))
}

# Fork only independent chromosomes; workers use one data.table thread each.
aw_parallel <- function(items, fun, cores=1L) {
  cores<-max(1L,min(as.integer(cores),length(items)))
  task<-function(item) {data.table::setDTthreads(1L);fun(item)}
  if(!length(items)) return(list())
  out<-if(cores>1L && .Platform$OS.type=='unix')
    parallel::mclapply(items,task,mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE) else lapply(items,task)
  if(length(out)!=length(items) || any(vapply(out,function(x)
      is.null(x) || inherits(x,'try-error'),logical(1)))) stop('Atlas export worker failed')
  out
}

# Only exact PRIMARY observations or VCF records explicitly used by a resolved
# rescue are removed. Proximity to a PRIMARY never proves record identity.
aw_candidate_context <- function(records, x) {
  r<-copy(records)
  if(!nrow(r)) return(r)
  r[,candidate_status:='NOT_EVALUATED']
  states<-as.data.table(x$locus_sample_states)
  loci<-as.data.table(x$loci)
  if(nrow(states) && nrow(loci)) {
    radius<-as.numeric(x$parameters$vcf_window)
    if(!length(radius) || !is.finite(radius)) stop('Missing VCF evaluation window in internal catalogue')
    spans<-merge(states[,.(sample_id,locus_id,state)],
      loci[,.(locus_id,CHROM,left=POS_MIN-radius,right=POS_MAX+radius)],by='locus_id',allow.cartesian=TRUE)
    setkey(spans,sample_id,CHROM,left,right)
    q<-r[,.(sample_id,CHROM,left=POS,right=POS,row_id=.I)]
    hits<-foverlaps(q,spans,nomatch=0L)
    if(nrow(hits)) {
      r[unique(hits$row_id),candidate_status:='EVALUATED_OTHER']
      r[unique(hits[state %in% c('PRIMARY','RESCUE_GT','RESCUE_ASSUMED')]$row_id),candidate_status:='NEAR_RESOLVED_LOCUS']
      r[unique(hits[state=='REVIEW']$row_id),candidate_status:='REVIEW']
    }
    used<-states[state=='RESCUE_GT' & !is.na(vcf_evidence_ids) & nzchar(vcf_evidence_ids),
      .(evidence_id=unlist(strsplit(vcf_evidence_ids,';',fixed=TRUE))),by=sample_id]
    if(nrow(used)) {
      keys<-paste(r$sample_id,paste(r$caller,r$CHROM,r$POS,r$SVLEN,sep=':'),sep='|')
      resolved<-keys %in% paste(used$sample_id,used$evidence_id,sep='|')
      r[resolved & candidate_status!='REVIEW',candidate_status:='RESOLVED_RESCUE']
    }
  }
  p<-as.data.table(x$primary_members)
  if(nrow(p)) {
    p<-p[!is.na(ALT_SEQUENCE) & grepl('^[ACGTRYSWKMBDHVN]+$',ALT_SEQUENCE,ignore.case=TRUE)]
    if(nrow(p)) {
      key<-function(sid,ch,pos,size,alt) paste(sid,ch,pos,size,toupper(alt),sep='|')
      represented<-key(r$sample_id,r$CHROM,r$POS,r$SVLEN,r$vcf_alt_seq) %in%
        key(p$sample_id,p$CHROM,p$POS,p$SVLEN,p$ALT_SEQUENCE)
      r[represented & candidate_status!='REVIEW',candidate_status:='RESOLVED_PRIMARY']
    }
  }
  r
}

aw_candidate_evidence <- function(records, detail='vectors') {
  # Aggregate identical evidence, never publish source/sample/read identifiers.
  # Counts across rows can overlap: this is not a sample-linked genotype matrix.
  if(detail!='vectors' || !nrow(records)) return(data.table())
  r<-records[window_callable==TRUE & !is.na(INS_GROUP)]
  if(!nrow(r)) return(data.table())
  result<-r[,.(N_SAMPLES=uniqueN(sample_id),N_RECORDS=uniqueN(record_key),
       SUPPORT_MIN=min(support),SUPPORT_MAX=max(support)),
    by=.(INS_GROUP,CHROM,POS,SVLEN,ALT_SEQUENCE=vcf_alt_seq,caller,source,
         GT,ploidy,dosage,genotype_qc,candidate_status)]
  setorderv(result,c('CHROM','POS','SVLEN','ALT_SEQUENCE','caller','source','GT','candidate_status'))
  result
}

aw_write <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(!ncol(x)) {con<-file(path,'wt');close(con)} else fwrite(x,path,sep='\t',na='NA')
}

aw_run <- function(spec) {
  cores<-if(is.null(spec$cores)) 1L else max(1L,as.integer(spec$cores))
  data.table::setDTthreads(1L)  # No nested threads in forked workers.
  started<-proc.time()[[3L]]
  message('[ATLAS EXPORT] allocated_cores=',cores)
  source(spec$helpers,local=TRUE)
  engine<-new.env(parent=environment())
  sys.source(spec$engine,envir=engine)
  engine$compare_script_dir<-dirname(spec$engine)
  sys.source(spec$sex_module,envir=engine)
  x<-readRDS(spec$internal);catalog<-readRDS(spec$catalog)
  m<-as.data.table(x$manifest);settings<-spec$settings;w<-as.integer(settings$window_bp)
  if(!identical(x$metadata$reference_build,spec$assembly)) stop('Atlas export assembly mismatch')
  mapping<-as.data.table(read_chromosome_map(spec$mapping))[included==TRUE]
  if(any(mapping$ASSEMBLY!=spec$assembly)) stop('Window mapping assembly mismatch')
  grid<-lapply(seq_len(nrow(mapping)),function(j) {
    n<-ceiling(mapping$expected_length[[j]]/w)
    data.table(CHROM=mapping$FINAL_CHR[[j]],start=(seq_len(n)-1)*w,
      end=pmin(seq_len(n)*w,mapping$expected_length[[j]]),N_SAMPLES=nrow(m),
      N_CALLABLE=0L,N_NOT_CALLABLE=0L,N_SEX_EXCLUDED=0L,N_ASSAYED=0L,N_PARTIAL=0L,N_PRIMARY_LOCI=0L)
  });names(grid)<-mapping$FINAL_CHR
  if(anyDuplicated(names(grid))) stop('Duplicate canonical chromosome in window mapping')
  regions<-lapply(grid,function(g) list(
    first=engine$sa_region(g$CHROM,g$start+1,spec$assembly),
    last=engine$sa_region(g$CHROM,g$end,spec$assembly)))
  # Ploidy depends on karyotype and coordinates, not on the individual sample.
  eligibility<-lapply(setNames(unique(m$sex_karyotype),unique(m$sex_karyotype)),function(karyo)
    lapply(regions,function(reg) {
      first<-engine$sa_ploidy(reg$first,karyo);last<-engine$sa_ploidy(reg$last,karyo)
      list(eligible=is.finite(first)&first>0&is.finite(last)&last>0,
           same_region=reg$first==reg$last)
    }))
  all_records<-list();audits<-list()
  for(i in seq_len(nrow(m))) {
    sid<-m$sample_id[[i]]
    cov<-read_coverage(m$quality_regions[[i]],mapping,spec$min_depth,spec$upper_ratio,m$sex[[i]])$regions
    cov[,ok:=cover_status=='PASS' & region_status=='PASS']
    coverage_parts<-aw_parallel(names(grid),function(ch) {
      g<-grid[[ch]];width<-g$end-g$start
      observed<-good<-numeric(nrow(g));pieces<-aw_pieces(cov[FINAL_CHR==ch],w)
      observed[pieces$bin]<-pieces$observed_bp;good[pieces$bin]<-pieces$callable_bp
      el<-eligibility[[m$sex_karyotype[[i]]]][[ch]]
      ok<-el$eligible & el$same_region & good/width>=settings$min_callable_fraction
      list(callable=ok,not_callable=el$eligible & !ok,sex_excluded=!el$eligible,
           assayed=observed>0,partial=observed>0 & observed<width)
    },cores)
    names(coverage_parts)<-names(grid)
    callable<-lapply(coverage_parts,`[[`,'callable')
    for(ch in names(grid)) {
      g<-grid[[ch]];part<-coverage_parts[[ch]]
      g[,`:=`(N_CALLABLE=N_CALLABLE+as.integer(part$callable),
        N_NOT_CALLABLE=N_NOT_CALLABLE+as.integer(part$not_callable),
        N_SEX_EXCLUDED=N_SEX_EXCLUDED+as.integer(part$sex_excluded),
        N_ASSAYED=N_ASSAYED+as.integer(part$assayed),N_PARTIAL=N_PARTIAL+as.integer(part$partial))]
    }
    rm(coverage_parts)
    sources<-spec$sources[[sid]]
    parsed_sources<-aw_parallel(sources,function(s) {
      retained<-data.table()
      cfg<-engine$cfg_defaults();cfg$vcf_require_pass<-settings$require_pass
      parsed<-engine$ra_read_vcf(s$path,sid,s$caller,cfg)
      if(!parsed$audit$status %in% c('OK','EMPTY')) stop('Invalid window VCF: ',s$path)
      d<-parsed$table
      if(nrow(d)) {
        d<-d[(is.na(SVTYPE)|SVTYPE=='INS') & is.finite(support) & support>=s$min_support & SVLEN>=spec$min_svlen]
        d<-d[CHROM %in% names(grid)]
      }
      if(nrow(d)) {
        d[,bin:=as.integer(floor((POS-1)/w)+1L)]
        d<-d[POS>=1 & POS<=mapping$expected_length[match(CHROM,mapping$FINAL_CHR)]]
        d[,`:=`(start=(bin-1)*w,end=pmin(bin*w,mapping$expected_length[match(CHROM,mapping$FINAL_CHR)]),source=s$source)]
        d[,window_callable:=callable[[.BY$CHROM]][bin],by=CHROM]
        d[,gt_status:=vapply(GT,function(gt) {
          if(is.na(gt)||!grepl('^[0-9]+([/|][0-9]+)*$',gt)) return('UNKNOWN')
          if(any(as.integer(strsplit(gt,'[/|]')[[1]])>0)) 'ALT' else 'REF'
        },character(1))]
        pl<-engine$sa_ploidy(engine$sa_region(d$CHROM,d$POS,spec$assembly),m$sex_karyotype[[i]])
        gt<-lapply(seq_len(nrow(d)),function(j) engine$ra_gt_info(d$GT[[j]],pl[[j]]))
        d[,`:=`(ploidy=as.integer(pl),dosage=vapply(gt,`[[`,integer(1),'dosage'),
          genotype_qc=vapply(gt,`[[`,character(1),'qc'))]
        # Same record in raw and rescue VCF must not double N_RECORDS.
        d[,record_key:=paste(sample_id,caller,CHROM,POS,SVLEN,ALT,sep='|')]
        retained<-d[,.(sample_id,caller,source,CHROM,POS,start,end,SVLEN,GT,vcf_alt_seq,
          support,VCF_ID,gt_status,window_callable,ploidy,dosage,genotype_qc,record_key)]
      }
      a<-copy(parsed$audit);a[,`:=`(source=s$source,min_support=s$min_support,n_retained=nrow(d))]
      list(records=retained,audit=a)
    },cores)
    all_records[[length(all_records)+1L]]<-rbindlist(lapply(parsed_sources,`[[`,'records'),fill=TRUE)
    audits[[length(audits)+1L]]<-rbindlist(lapply(parsed_sources,`[[`,'audit'))
    rm(parsed_sources)
    message('[ATLAS WINDOWS] sample ',i,'/',nrow(m),' ',sid)
  }
  # Locus spans include every original PRIMARY position, never a shifted consensus.
  for(ch in names(grid)) {
    g<-grid[[ch]];l<-as.data.table(x$loci)[CHROM==ch]
    if(nrow(l)) {
      first<-floor((l$POS_MIN-1)/w)+1L;last<-floor((l$POS_MAX-1)/w)+1L
      bins<-unlist(Map(seq.int,first,last),use.names=FALSE)
      g[,N_PRIMARY_LOCI:=tabulate(bins,nbins=nrow(g))]
    }
  }
  records<-if(length(all_records)) rbindlist(all_records) else data.table()
  records<-aw_candidate_context(records,x)
  pending<-if(nrow(records)) records[!candidate_status %in% c('RESOLVED_PRIMARY','RESOLVED_RESCUE')] else records
  merged<-aw_summarize_ins(pending,settings,cores);bins<-merged$summary
  setorder(bins,CHROM,start,end,INS_GROUP)
  if(nrow(records) && nrow(merged$records)) {
    records[,INS_GROUP:=NA_character_]
    records[unique(merged$records[,.(record_key,source,INS_GROUP)]),on=.(record_key,source),INS_GROUP:=i.INS_GROUP]
  }
  windows<-rbindlist(grid)
  windows[,`:=`(N_INS_SAMPLES=0L,INS_SUMMARY='',N_INS_GROUPS=0L)]
  if(nrow(bins)) {
    agg<-bins[,.(N_INS_GROUPS=.N,INS_SUMMARY=paste0('INS ',format(SVLEN_MEDIAN,trim=TRUE,scientific=FALSE),':',N_SAMPLES,collapse='; ')),by=.(CHROM,start,end)]
    windows[agg,on=.(CHROM,start,end),`:=`(N_INS_GROUPS=i.N_INS_GROUPS,INS_SUMMARY=i.INS_SUMMARY)]
    people<-pending[window_callable==TRUE,.(N_INS_SAMPLES=uniqueN(sample_id)),by=.(CHROM,start,end)]
    windows[people,on=.(CHROM,start,end),N_INS_SAMPLES:=i.N_INS_SAMPLES]
  }
  # Preserve every mapped window, including no coverage and no VCF calls.
  stopifnot(all(windows$N_CALLABLE+windows$N_NOT_CALLABLE+windows$N_SEX_EXCLUDED==windows$N_SAMPLES))
  windows[,`:=`(ASSEMBLY=spec$assembly,NO_PRIMARY=N_PRIMARY_LOCI==0L,
    EVIDENCE_SCOPE='ORIGINAL_ALL_REGIONS_PLUS_REGIONAL_RESCUE')]
  candidates<-merge(bins,windows[,.(CHROM,start,end,N_CALLABLE,N_NOT_CALLABLE,N_SEX_EXCLUDED,N_PRIMARY_LOCI)],
    by=c('CHROM','start','end'),all.x=TRUE,sort=FALSE)
  candidates[,`:=`(ASSEMBLY=rep(spec$assembly,.N),INTERPRETATION=rep('INS_CANDIDATE_NOT_CONFIRMED_TE_OR_RESCUE',.N))]
  # These tables are a backup, not new population states or allele frequencies.
  public<-spec$public;private<-spec$private
  dir.create(public,recursive=TRUE,showWarnings=FALSE);dir.create(private,recursive=TRUE,showWarnings=FALSE)
  aw_write(records,file.path(private,'window_evidence_PRIVATE.tsv'))
  aw_write(rbindlist(audits),file.path(private,'window_sources_PRIVATE.tsv'))
  summary<-copy(catalog$clusters)
  summary[,`:=`(start=pmax(0,POS_MIN-1),end=POS_MAX,ASSEMBLY=rep(spec$assembly,.N))]
  # TSV is a concise summary; full synchronized vectors remain in the RDS.
  listcols<-names(summary)[vapply(summary,is.list,logical(1))]
  if(length(listcols)) summary[,(listcols):=NULL]
  aw_write(summary,file.path(public,'population_summary.tsv'))
  beds<-summary[,.(CHROM,start,end,locus_id)]
  fwrite(beds,file.path(public,'population_loci.bed'),sep='\t',col.names=FALSE)
  states<-as.data.table(x$locus_sample_states)
  rescues<-states[state %in% c('RESCUE_GT','RESCUE_ASSUMED')]
  rescues<-merge(rescues,as.data.table(x$loci)[,.(locus_id,CHROM,POS,SVLEN)],by='locus_id',all.x=TRUE,sort=FALSE)
  aw_write(rescues,file.path(private,'rescues_PRIVATE.tsv'))
  rescue_summary<-rescues[,.(N_SAMPLES=uniqueN(sample_id)),by=.(locus_id,state,rescue_type)]
  aw_write(rescue_summary,file.path(public,'rescues_summary.tsv'))
  catalog$window_backup<-list(format_version='RetroAtlas_Interpopulation_Backup_v2',settings=settings,source_scope='ORIGINAL_ALL_REGIONS_PLUS_REGIONAL_RESCUE',
    original_min_support=spec$primary_support,rescue_min_support=spec$rescue_support,
    windows=windows,ins_candidates=candidates,
    candidate_evidence=aw_candidate_evidence(merged$records,catalog$metadata$detail),
    detail=catalog$metadata$detail,
    selection='Callable VCF INS excluding exact PRIMARY or explicitly used resolved VCF rescues; REVIEW retained',
    evaluation_context='Nearby population states within VCF search range are context, not proof of allele identity or exhaustive window evaluation',
    callable_only=TRUE,allele_counts_scope='Observed concordant VCF genotypes among candidate samples only, excluding REVIEW; not population AF',
    interpretation='No VCF call is not confirmed reference. Size bins are candidates, not confirmed TE/rescues or population AF. With vectors detail, candidate_evidence preserves exact positions and available inserted sequences, with no cross-row sample linkage.')
  catalog$metadata$coordinates<-'BED0_HALF_OPEN; insertion POS is 1-based'
  catalog$rescue_summary<-rescue_summary
  saveRDS(catalog,file.path(public,'Atlas_merge.rds'),compress='gzip')
  message('[ATLAS EXPORT] complete | seconds=',round(proc.time()[[3L]]-started,2),
    ' | windows=',nrow(windows),' | candidate_groups=',nrow(candidates))
  invisible(list(internal=x,catalog=catalog,windows=windows,candidates=candidates,rescues=rescues))
}
