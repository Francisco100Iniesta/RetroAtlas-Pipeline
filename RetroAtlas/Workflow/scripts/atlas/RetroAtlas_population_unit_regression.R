#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
this <- local({
  z<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  if(length(z)) dirname(normalizePath(sub('^--file=','',z[[1L]]))) else getwd()
})
source(file.path(this,'RetroAtlas_population_merge.R'),local=TRUE)
e<-new.env(parent=globalenv())
sys.source(file.path(this,'../compare/RetroAtlas_sexaware_frequency.R'),envir=e)
e$compare_script_dir<-normalizePath(file.path(this,'../compare'))
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
f<-pm_frequency(st,loc[1L],meta,e,cfg)
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

cat('[POP_TEST] PAR boundaries, original coordinates and unchanged Compare ploidy\n')
for(build in c('GRCh37','GRCh38','T2T-CHM13v2.0')) {
 r<-e$sa_ranges(build)
 for(ch in c('X','Y')) for(j in 1:2) {
  pp<-c(r[[ch]][j,1]-1,r[[ch]][j,1],r[[ch]][j,2],r[[ch]][j,2]+1)
  stopifnot(identical(e$sa_par(rep(paste0('chr',ch),4),pp,build),
    c('NO_PAR',paste0('PAR',j),paste0('PAR',j),'NO_PAR')))
 }
}
stopifnot(is.na(e$sa_ploidy('Y_PAR','XY')),
 pm_detection_eligible('Y_PAR','XY',NA_integer_),
 !pm_detection_eligible('Y_PAR','XX',NA_integer_),
 !pm_detection_eligible('Y_NONPAR','XX',0L))
loc2<-data.table(locus_id=c('AUTO','XP','YP','XN','YN'),
 CHROM=c('chr1','chrX','chrY','chrX','chrY'),POS=c(1000,10001,10001,3000000,3000000))
m2<-data.table(sample_id=c('M','F','OTHER'),sex=c('M','F','M'),
 sex_karyotype=c('XY','XX','XXY'),origin='EU')
st2<-rbindlist(lapply(loc2$locus_id,function(lid)
 rbindlist(lapply(m2$sample_id,function(sid) pm_empty_state(lid,sid)))))
st2[,state:='PRIMARY'];st2[,GT_RAW:='0/1']
st2[,sex_region:=e$sa_region(loc2$CHROM[match(locus_id,loc2$locus_id)],
 loc2$POS[match(locus_id,loc2$locus_id)],cfg$reference_build)]
st2[,ploidy:=e$sa_ploidy(sex_region,m2$sex_karyotype[match(sample_id,m2$sample_id)])]
st2[,`:=`(GT=NA_character_,dosage=NA_integer_,genotype_qc='GT_PLOIDY_CONFLICT')]
st2[ploidy==2L,`:=`(GT='0/1',dosage=1L,genotype_qc='OK')]
st2[sample_id=='F' & locus_id=='YN',`:=`(state='SEX_EXCLUDED',genotype_qc=NA_character_)]
st2[sample_id=='M' & locus_id=='XN',GT_RAW:='1|0']
before<-copy(st2)
f2<-pm_frequency(st2,loc2,m2,e,cfg)
q<-f2$summary
cols<-grep('^(AC|AN|AF)_',names(q),value=TRUE)
stopifnot(identical(st2,before),all(is.na(unlist(q[PAR!='NO_PAR',..cols]))),
 all(q[locus_id %in% c('XN','YN'),N_PRIMARY_HAPLOID_ASSUMED]==1L),
 all(q[locus_id %in% c('AUTO','XP','YP'),N_PRIMARY_HAPLOID_ASSUMED]==0L),
 q[locus_id=='AUTO',AF_PRIMARY]==0.5,
 q[locus_id=='XN',AC_PRIMARY]==1L,q[locus_id=='XN',AN_PRIMARY]==2L,
 q[locus_id=='XN',AF_PRIMARY_HAPLOID_ASSUMED]==2/3,
 q[locus_id=='YN',AF_PRIMARY_HAPLOID_ASSUMED]==1,
 identical(f2$private_states[PAR!='NO_PAR',GT],f2$private_states[PAR!='NO_PAR',GT_RAW]))
sc<-grep('^(AC|AN|AF)_',names(f2$private_strata),value=TRUE)
stopifnot(all(is.na(unlist(f2$private_strata[PAR!='NO_PAR',..sc]))))
# A genotype that is missing, reference, multiallelic or already resolved cannot be assumed.
for(gt in c(NA_character_,'./.','0/0','1/2','1/1')) {
 ss<-copy(st2);ss[sample_id=='M' & locus_id=='XN',GT_RAW:=gt]
 ff<-pm_frequency(ss,loc2,m2,e,cfg)
 stopifnot(ff$summary[locus_id=='XN',N_PRIMARY_HAPLOID_ASSUMED]==0L)
}
ss<-copy(st2);ss[sample_id=='M' & locus_id=='XN',state:='REVIEW']
stopifnot(pm_frequency(ss,loc2,m2,e,cfg)$summary[locus_id=='XN',N_PRIMARY_HAPLOID_ASSUMED]==0L)
stopifnot(inherits(try(pm_frequency(rbind(st2,st2[1]),loc2,m2,e,cfg),silent=TRUE),'try-error'))
empty<-pm_frequency(st2[0],loc2[0],m2,e,cfg)
stopifnot(nrow(empty$summary)==0L,'PAR' %in% names(empty$summary))
cat('[POP_TEST] ALL PAR AND HAPLOID ASSUMPTION TESTS PASS\n')

cat('[POP_TEST] Real loader and rescue path: Y-PAR retained, QC gate unchanged\n')
sys.source(file.path(this,'../compare/RetroAtlas_compare.R'),envir=e)
sys.source(file.path(this,'../reference_helpers.R'),envir=e)
sys.source(file.path(this,'../compare/RetroAtlas_sexaware_frequency.R'),envir=e)
e$compare_script_dir<-normalizePath(file.path(this,'../compare'))
tmp<-tempfile('atlas-par-');dir.create(tmp)
map<-data.table(CHROM=c('chrX','chrY'),FINAL_CHR=c('chrX','chrY'),
 ASSEMBLY='GRCh38.p14',expected_length=c(156040895,57227415),included=TRUE,status='PASS')
fwrite(map,file.path(tmp,'map.tsv'),sep='\t')
fwrite(data.table(CHROM=c('chrX','chrY'),start=10000,end=30000,depth=20),
 file.path(tmp,'coverage.bed'),sep='\t',col.names=FALSE)
input<-data.table(CHROM=c('chrX','chrY'),POS=20001,svlen=100,
 GT='0|1',ALT_SUPPORT_MAX=6,event_class='Alu',event='single',completeness='full',
 ALT=strrep('A',100),cover_status='PASS',region_status='PASS')
saveRDS(list(loci=input),file.path(tmp,'input.rds'))
cf<-cfg;cf$chromosome_map<-file.path(tmp,'map.tsv')
row<-data.table(sample_id='M',sex='M',sex_karyotype='XY',
 final_rds=file.path(tmp,'input.rds'),quality_regions=file.path(tmp,'coverage.bed'))
loaded<-pm_load_sample(row,e,cf)$loci
stopifnot(nrow(loaded)==2L,all(loaded$primary_eligible),identical(loaded$GT,input$GT))
edge0<-data.table(plausible=logical(),unresolved_nearby=logical(),locus_a=character(),locus_b=character())
source<-as.list(loaded[CHROM=='chrY']);source$seq_valid<-FALSE
locus<-data.table(locus_id='YP',CHROM='chrY',POS=20001)
quality<-data.table(callable=TRUE,quality_status='CALLABLE',frac_ge20=1,mean_depth=20)
vcf<-data.table(CHROM='chrY',POS=20001,SVLEN=100,caller='sniffles2',GT='0|1',support=6,vcf_alt_seq=NA_character_)
cf$clip_enabled<-FALSE;cf$fragment_enabled<-FALSE
res<-pm_one_state(locus,source,'M',row,loaded[0],edge0,quality,vcf,data.table(),NULL,'',e,cf,FALSE)
stopifnot(res$state=='RESCUE_GT',res$GT_RAW=='0|1')
fr<-pm_frequency(res,locus,m2[sample_id=='M'],e,cf)
stopifnot(fr$private_states$GT=='0|1',is.na(fr$private_states$dosage),is.na(fr$summary$AF_GT))
quality[,`:=`(callable=FALSE,quality_status='LOW_DEPTH',mean_depth=2)]
res<-pm_one_state(locus,source,'M',row,loaded[0],edge0,quality,vcf,data.table(),NULL,'',e,cf,FALSE)
stopifnot(res$state=='NOT_CALLABLE',res$quality_gate_skipped)
pm_regional_gate_preflight()
unlink(tmp,recursive=TRUE)
cat('[POP_TEST] LOADER, PAR RESCUE AND QC GATE PASS\n')

# End-to-end small export: public catalogue, candidate backup and report payload.
local({
workflow_dir<-normalizePath(file.path(this,'../..'))
source(file.path(this,'atlas_outputs.R'))
t<-tempfile('atlas-export-');dir.create(t);on.exit(unlink(t,recursive=TRUE),add=TRUE)
map<-data.table(CHROM=c('chrX','chrY'),FINAL_CHR=c('chrX','chrY'),ASSEMBLY='GRCh38.p14',expected_length=30000,included=TRUE,status='PASS')
fwrite(map,file.path(t,'mapping.tsv'),sep='\t')
fwrite(data.table(CHROM=c('chrX','chrY'),start=0,end=30000,depth=20),file.path(t,'coverage.bed'),sep='\t',col.names=FALSE)
p<-data.table(sample_id='M',uid='M::1',locus_id='YP',CHROM='chrY',POS=20001,
 SVLEN=100,GT='0|1',ALT_SEQUENCE=strrep('A',100),ploidy=NA_integer_,PAR='PAR1',
 event_class='Alu',event='single',completeness='full',element_type='Alu',
 polyA=NA_character_,polyA_length=NA_real_,tsd_detected=NA_character_,
 tsd_length=NA_real_,tsd_sequence=NA_character_)
m<-data.table(sample_id='M',sex='M',sex_karyotype='XY',origin='EU',quality_regions=file.path(t,'coverage.bed'))
l<-data.table(locus_id='YP',CHROM='chrY',POS=20001,POS_MIN=20001,POS_MAX=20001,SVLEN=100,N_PRIMARY_SAMPLES=1L,CLASS_MODE='Alu',EVENT_MODE='single',Q_MIN=NA_integer_)
s<-pm_empty_state('YP','M');s[,`:=`(state='PRIMARY',GT_RAW='0|1',GT=NA_character_,sex_region='Y_PAR',callable=TRUE)]
cfg<-pm_cfg_defaults();cfg$reference_build<-'GRCh38.p14';cfg$public_detail<-'vectors'
f<-pm_frequency(s,l,m,e,cfg);catobj<-pm_catalog(f$summary,p,m,f$public_strata,cfg)
x<-list(metadata=list(reference_build=cfg$reference_build,n_samples=1L,n_loci=1L,catalog_id='TEST'),manifest=m,loci=f$summary,locus_sample_states=f$private_states,primary_members=p,parameters=cfg,neighboring_loci=data.table())
saveRDS(x,file.path(t,'internal.rds'));saveRDS(catobj,file.path(t,'catalog.rds'))
vcf<-c('##fileformat=VCFv4.2','#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM',
 paste('chrY',22001,'v1','A',paste0('A',strrep('C',100)),'.','PASS','SVTYPE=INS;SVLEN=100','GT:DV','0|1:6',sep='\t'),
 paste('chrX',22001,'v2','A',paste0('A',strrep('C',100)),'.','PASS','SVTYPE=INS;SVLEN=100','GT:DV','0/1:6',sep='\t'),
 paste('chrX',5001,'v3','A',paste0('A',strrep('C',100)),'.','PASS','SVTYPE=INS;SVLEN=100','GT:DV','1:6',sep='\t'))
writeLines(vcf,file.path(t,'calls.vcf'))
spec<-list(internal=file.path(t,'internal.rds'),catalog=file.path(t,'catalog.rds'),
 helpers=normalizePath(file.path(workflow_dir,'scripts/reference_helpers.R')),engine=normalizePath(file.path(workflow_dir,'scripts/compare/RetroAtlas_compare.R')),sex_module=normalizePath(file.path(workflow_dir,'scripts/compare/RetroAtlas_sexaware_frequency.R')),
 mapping=file.path(t,'mapping.tsv'),assembly=cfg$reference_build,min_depth=10,upper_ratio=3,min_svlen=50,primary_support=3,rescue_support=1,
 settings=list(window_bp=1000L,min_callable_fraction=0.8,require_pass=TRUE,min_size_similarity=0.8,max_size_difference_bp=50),
 sources=list(M=list(list(path=file.path(t,'calls.vcf'),caller='sniffles2',source='original',min_support=3))),
 public=file.path(t,'public'),private=file.path(t,'private'),cores=1,
 module=normalizePath(file.path(workflow_dir,'scripts/atlas/atlas_outputs.R')),template=normalizePath(file.path(workflow_dir,'scripts/atlas/report_template.html')),page_size=50,public_link='public')
r<-aw_run(spec)
stopifnot(nrow(r$candidates)==3L,
 all(is.na(r$candidates[PAR!='NO_PAR',AC_OBSERVED])),all(is.na(r$candidates[PAR!='NO_PAR',AN_OBSERVED])),
 r$candidates[CHROM=='chrY',N_CALLABLE]==1L,
 r$candidates[PAR=='NO_PAR',AC_OBSERVED]==1L,
 identical(r$catalog$window_backup$candidate_evidence[CHROM=='chrY',GT],'0|1'),
 all(is.na(r$catalog$window_backup$candidate_evidence[PAR!='NO_PAR',dosage])))
jsonlite::write_json(spec,file.path(t,'spec.json'),auto_unbox=TRUE)
status<-system2('Rscript',shQuote(c(file.path(workflow_dir,'scripts/atlas/export_atlas.R'),file.path(t,'spec.json'))))
stopifnot(status==0L,file.exists(file.path(t,'private/atlas_report.html')))
html<-paste(readLines(file.path(t,'private/atlas_report.html')),collapse='\n')
stopifnot(grepl('AF_PRIMARY_HAPLOID_ASSUMED',html,fixed=TRUE),!grepl('PAYLOAD_BASE64',html,fixed=TRUE))
validator<-file.path(this,'RetroAtlas_population_validate.R')
args<-shQuote(c(validator,file.path(t,'internal.rds'),file.path(t,'catalog.rds')))
vlog<-file.path(t,'validation.log')
stopifnot(system2('Rscript',args,stdout=vlog,stderr=vlog)==0L)
bad<-readRDS(file.path(t,'catalog.rds'));bad$clusters$AC_PRIMARY<-0L
saveRDS(bad,file.path(t,'catalog.rds'))
stopifnot(suppressWarnings(system2('Rscript',args,stdout=vlog,stderr=vlog))!=0L,
 any(grepl('POPULATION_VALIDATION_FAILED',readLines(vlog),fixed=TRUE)))
cat('EXPORT, PAR BACKUP, ORIGINAL GT, HTML AND VALIDATOR PASS\n')

})
