#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
this <- local({
  z<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(z)) dirname(normalizePath(sub('^--file=','',z[[1L]]))) else getwd()
})
source(file.path(this,'RetroAtlas_population_merge.R'),local=TRUE)
cat('[POP_TEST] Testing complete linkage (no transitive chain)\n')
p<-data.table(uid=c('A::1','B::1','C::1'),sample_id=c('A','B','C'),
 CHROM='chr1',POS=c(1000,1400,1800),SVLEN=100,GT='0/1',alt_support=10,
 RQ_SCORE=5,event_class='Alu',event='single',completeness='full')
ed<-data.table(uid_a=c('A::1','B::1'),uid_b=c('B::1','C::1'),
 sample_a=c('A','B'),sample_b=c('B','C'),Q=c(4L,4L),score=c(4000001,4000000),
 dpos=c(400,400),size_sim=c(1,1),seq_status='NOT_EVALUATED',
 seq_identity=NA_real_,seq_coverage=NA_real_,reason='TEST')
cfg<-pm_cfg_defaults()
c1<-pm_cluster(p,ed,cfg)
stopifnot(nrow(c1$loci)==2L,max(c1$loci$N_PRIMARY)==2L,
  !anyDuplicated(c1$primary[,.(locus_id,sample_id)]))
cat('[POP_TEST] Testing same-sample alternative PRIMARY rejection\n')
p2<-copy(p[1:2]);p2<-rbind(p2,p2[1L])
p2[3L,`:=`(uid='A::2',POS=1005)]
e2<-data.table(uid_a=c('A::1','A::2'),uid_b=c('B::1','B::1'),
 sample_a='A',sample_b='B',Q=c(10L,9L),score=c(10000000,9000000),
 dpos=c(0,5),size_sim=c(1,1),seq_status='NOT_EVALUATED',
 seq_identity=NA_real_,seq_coverage=NA_real_,reason='TEST')
c2<-pm_cluster(p2,e2,cfg)
stopifnot(nrow(c2$loci)==2L,any(c2$edge_decisions$decision=='SAME_SAMPLE_COMPETITION'))
cat('[POP_TEST] Testing one state per locus/sample and no double-counting\n')
loc<-c1$loci
st<-rbindlist(list(pm_empty_state(loc$locus_id[[1L]],'A'),
  pm_empty_state(loc$locus_id[[1L]],'B'),
  pm_empty_state(loc$locus_id[[1L]],'C')))
st[1L,`:=`(state='PRIMARY',ploidy=2L,dosage=2L)]
st[2L,`:=`(state='REF_CALLABLE',ploidy=2L)]
st[3L,`:=`(state='RESCUE_ASSUMED',ploidy=2L,assumption_one_alt=TRUE)]
meta<-data.table(sample_id=c('A','B','C'),sex='F',sex_karyotype='XX',
 ancestry_group='EU',ethnicity='EU',origin='EU')
f<-pm_frequency(st,loc[1L],meta,NULL,cfg)
stopifnot(nrow(f$summary)==1L,
  f$summary$AC_PRIMARY==2L,f$summary$AN_PRIMARY==4L,
  f$summary$AC_GT==2L,f$summary$AN_GT==4L,
  f$summary$AC_ADJUSTED==3L,f$summary$AN_ADJUSTED==6L,
  abs(f$summary$AF_RESCUE_ADJUSTED-0.5)<1e-8)
cat('[POP_TEST] Testing rescue competition across linked loci\n')
d<-rbindlist(list(pm_empty_state('L1','X'),pm_empty_state('L2','X')))
d[,`:=`(state='RESCUE_GT',rescue_type='VCF',dosage=1L,GT='0/1',
  vcf_evidence_ids=c('caller1:chr1:100:100','caller2:chr1:110:100'))]
lk<-data.table(locus_a='L1',locus_b='L2',plausible=TRUE)
g<-pm_resolve_reuse(d,lk)
stopifnot(all(g$states$state=='REVIEW'),nrow(g$competing_rescues)==2L)
cat('[POP_TEST] Testing reused VCF record across two unlinked loci\n')
d[,vcf_evidence_ids:='caller1:chr1:100:100']
g<-pm_resolve_reuse(d,data.table())
stopifnot(all(g$states$state=='REVIEW'),nrow(g$reused_evidence)==1L)
cat('[POP_TEST] ALL CORE TESTS PASS\n')
