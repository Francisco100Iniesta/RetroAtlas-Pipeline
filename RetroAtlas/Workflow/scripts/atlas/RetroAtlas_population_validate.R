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
 all(l$AC_PRIMARY<=l$AN_PRIMARY & l$AC_GT<=l$AN_GT & l$AC_ADJUSTED<=l$AN_ADJUSTED),
 all(s$state %in% c('PRIMARY','RESCUE_GT','RESCUE_ASSUMED','REVIEW',
                    'NOT_CALLABLE','REF_CALLABLE','SEX_EXCLUDED')),
 !any(grepl('sample_id|row_uid|source_ids|VCF_ID|CLIP_ID',names(pub),ignore.case=TRUE)),
 identical(c$metadata$privacy,'DEIDENTIFIED_NOT_ANONYMOUS')))
print(checks)
cat('Samples:',ns,'Loci:',nl,'PRIMARY:',sum(s$state=='PRIMARY'),
 'RESCUE_GT:',sum(s$state=='RESCUE_GT'),
 'RESCUE_ASSUMED:',sum(s$state=='RESCUE_ASSUMED'),
 'REVIEW:',sum(s$state=='REVIEW'),'\n')
if(any(!checks$pass)) stop('POPULATION_VALIDATION_FAILED: ',paste(checks$check[!checks$pass],collapse=';'))
cat('POPULATION_VALIDATION_OK\n')
