#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
argv<-commandArgs(trailingOnly=TRUE)
if(length(argv)!=2L) stop('Usage: Rscript RetroAtlas_population_validate.R INTERNAL.rds Atlas_merge.rds')
i<-readRDS(argv[[1L]]);c<-readRDS(argv[[2L]])
p<-as.data.table(i$primary_members); s<-as.data.table(i$locus_sample_states)
l<-as.data.table(i$loci); pub<-as.data.table(c$clusters)
ns<-nrow(i$manifest);nl<-nrow(l)
checks<-data.table(check=c(
 'state matrix N x L','unique locus+sample','unique PRIMARY UID',
 'one PRIMARY per sample/locus','PRIMARY member count equals state count',
 'catalog locus identities','all catalog AF in [0,1]',
 'frequency allele counts satisfy AC <= AN',
 'one state belongs to supported categories',
 'no sample IDs in public field names','public metadata privacy warning'),
 pass=c(
 nrow(s)==ns*nl,
 anyDuplicated(s[,.(locus_id,sample_id)])==0L,
 anyDuplicated(p$uid)==0L,
 anyDuplicated(p[,.(locus_id,sample_id)])==0L,
 nrow(p)==sum(s$state=='PRIMARY'),
 setequal(l$locus_id,pub$locus_id),
 all(unlist(pub[,.(AF_PRIMARY,AF_GT,AF_RESCUE_ADJUSTED)],use.names=FALSE) %in% c(NA_real_,0) |
       (unlist(pub[,.(AF_PRIMARY,AF_GT,AF_RESCUE_ADJUSTED)],use.names=FALSE)>=0 &
        unlist(pub[,.(AF_PRIMARY,AF_GT,AF_RESCUE_ADJUSTED)],use.names=FALSE)<=1),na.rm=TRUE),
 all(l[PAR=='NO_PAR',AC_PRIMARY<=AN_PRIMARY & AC_GT<=AN_GT & AC_ADJUSTED<=AN_ADJUSTED]),
 all(s$state %in% c('PRIMARY','RESCUE_GT','RESCUE_ASSUMED','REVIEW',
                    'NOT_CALLABLE','REF_CALLABLE','SEX_EXCLUDED')),
 !any(grepl('sample_id|row_uid|source_ids|VCF_ID|CLIP_ID',names(pub),ignore.case=TRUE)),
 identical(c$metadata$privacy,'DEIDENTIFIED_NOT_ANONYMOUS')))
allele_cols<-grep('^(AC|AN|AF)_',names(l),value=TRUE)
add_check<-function(name,ok) checks<<-rbind(checks,data.table(check=name,pass=isTRUE(ok)))
add_check('PAR annotation complete',all(l$PAR %in% c('PAR1','PAR2','NO_PAR')))
add_check('PAR has no population allele statistics',all(is.na(unlist(l[PAR!='NO_PAR',..allele_cols]))))
add_check('Public PAR has no population allele statistics',all(is.na(unlist(pub[PAR!='NO_PAR',..allele_cols]))))
strata<-as.data.table(c$strata)
acols<-grep('^(AC|AN|AF)_',names(strata),value=TRUE)
add_check('PAR strata have no population allele statistics',all(is.na(unlist(strata[PAR!='NO_PAR',..acols]))))
expected<-s[,.(conflicts=sum(state=='PRIMARY' & genotype_qc %in% 'GT_PLOIDY_CONFLICT'),
  assumed=sum(haploid_alt_assumed)),by=locus_id]
e<-expected[match(l$locus_id,locus_id)]
add_check('Conflict and assumption counts match states',
 all(l$N_PRIMARY_GT_PLOIDY_CONFLICT==e$conflicts & l$N_PRIMARY_HAPLOID_ASSUMED==e$assumed))
q<-l[PAR=='NO_PAR']
add_check('Alternative allele counts add exactly one per eligible PRIMARY',
 all(q$AC_PRIMARY_HAPLOID_ASSUMED==q$AC_PRIMARY+q$N_PRIMARY_HAPLOID_ASSUMED &
     q$AN_PRIMARY_HAPLOID_ASSUMED==q$AN_PRIMARY+q$N_PRIMARY_HAPLOID_ASSUMED &
     q$AC_PRIMARY_HAPLOID_ASSUMED<=q$AN_PRIMARY_HAPLOID_ASSUMED))
expected_af<-ifelse(q$AN_PRIMARY_HAPLOID_ASSUMED>0,
 q$AC_PRIMARY_HAPLOID_ASSUMED/q$AN_PRIMARY_HAPLOID_ASSUMED,NA_real_)
add_check('Alternative AF matches alternative AC/AN',isTRUE(all.equal(q$AF_PRIMARY_HAPLOID_ASSUMED,expected_af)))
a<-s[haploid_alt_assumed==TRUE];mi<-match(a$sample_id,i$manifest$sample_id)
add_check('Assumption restricted to XY male non-PAR conflicting heterozygous PRIMARY',
 all(a$state=='PRIMARY' & a$PAR=='NO_PAR' & a$ploidy==1L &
 a$sex_region %in% c('X_NONPAR','Y_NONPAR') & a$genotype_qc=='GT_PLOIDY_CONFLICT' &
 is.na(a$dosage) & a$GT_RAW %in% c('0/1','1/0','0|1','1|0') &
 i$manifest$sex[mi]=='M' & i$manifest$sex_karyotype[mi]=='XY'))
j<-match(paste(p$locus_id,p$sample_id),paste(s$locus_id,s$sample_id))
add_check('PRIMARY input GT preserved',identical(as.character(p$GT),as.character(s$GT_RAW[j])))
add_check('PAR state dosage unavailable',all(is.na(s[PAR!='NO_PAR',dosage])))
print(checks)
cat('Samples:',ns,'Loci:',nl,'PRIMARY:',sum(s$state=='PRIMARY'),
 'RESCUE_GT:',sum(s$state=='RESCUE_GT'),
 'RESCUE_ASSUMED:',sum(s$state=='RESCUE_ASSUMED'),
 'REVIEW:',sum(s$state=='REVIEW'),'\n')
if(any(!checks$pass)) stop('POPULATION_VALIDATION_FAILED: ',paste(checks$check[!checks$pass],collapse=';'))
cat('POPULATION_VALIDATION_OK\n')
