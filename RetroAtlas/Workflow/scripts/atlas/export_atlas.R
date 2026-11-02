#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
spec<-jsonlite::fromJSON(commandArgs(TRUE)[[1]],simplifyVector=FALSE)
source(spec$module)
result<-aw_run(spec)
x<-result$internal
# Offline report: no ALT/read sequences, embedded data and no remote dependencies.
pack<-aw_pack
state_cols<-c('locus_id','sample_id','state','reason','GT','dosage','ploidy','callable','rescue_type','source_sample','linked_loci')
primary_cols<-c('locus_id','sample_id','uid','CHROM','POS','SVLEN','GT','event_class','event','completeness','alt_support')
payload<-list(meta=x$metadata,page_size=spec$page_size,
 loci=pack(x$loci,c('locus_id','CHROM','POS','POS_MIN','POS_MAX','SVLEN','CLASS_MODE','N_PRIMARY','N_RESCUE_GT','N_RESCUE_ASSUMED','N_REVIEW','N_NOT_CALLABLE','N_REF_CALLABLE','AF_PRIMARY','AF_GT','AF_RESCUE_ADJUSTED','locus_interpretation')),
 states=pack(x$locus_sample_states,state_cols),primary=pack(x$primary_members,primary_cols),
 links=pack(x$neighboring_loci),rescues=pack(result$rescues,c(state_cols,'CHROM','POS','SVLEN')),
 candidates=pack(result$candidates),samples=as.character(x$manifest$sample_id),
 windows=list(n=nrow(result$windows),without_primary=sum(result$windows$NO_PRIMARY),
   with_ins=sum(result$windows$N_INS_SAMPLES>0),fully_callable=sum(result$windows$N_CALLABLE==nrow(x$manifest)),
   bp=spec$settings$window_bp),public_link=spec$public_link)
json<-jsonlite::toJSON(payload,auto_unbox=TRUE,na='null',null='null',digits=15)
b64<-jsonlite::base64_enc(charToRaw(enc2utf8(json)))
html<-paste(readLines(spec$template,warn=FALSE),collapse='\n')
html<-sub('PAYLOAD_BASE64',b64,html,fixed=TRUE)
writeLines(html,file.path(spec$private,'atlas_report.html'),useBytes=TRUE)
message('[ATLAS OUTPUTS] Public catalogue, window backup and private report complete')
