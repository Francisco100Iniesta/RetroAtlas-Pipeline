#!/usr/bin/env Rscript
# ATLAS Intra-sample Insertion Explorer
# A self-contained HTML and a standards-oriented symbolic insertion VCF.
# Dependency: jsonlite. No web services are required.
# Usage:
# Rscript report_final.R --rds sample.rds --mane mane.bed \
#   --regions-bed sample.regions.bed.gz \
#   --out atlas_report_final.html --vcf atlas_report_final.vcf

parse_args <- function(argv) {
  o <- list(rds = NULL, mane = NULL, regions_bed = NULL,
            chromosome_map=NULL, chromosome_lengths=NULL, catalog=NULL, helpers=NULL, template=NULL,
            tsv=NULL, mane_assembly=NULL, primary_min_support='3', min_svlen='30', min_depth='0', upper_depth_ratio='3', panel_bed=NULL,
            out = 'atlas_report_final.html', vcf = NULL, readme = NULL, run_config='{}',
            reference = 'unknown', promoter_bp = '2000')
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (a %in% c('-h','--help')) {
      cat('Usage: Rscript report_final.R --rds input.rds --mane mane.bed [--regions-bed regions.bed.gz] [--out report.html] [--vcf report.vcf] [--readme additional.txt] [--reference GRCh38] [--promoter-bp 2000]\n')
      quit(save='no',status=0)
    }
    if (!startsWith(a,'--')) stop('Unexpected option: ',a)
    if (grepl('=',a,fixed=TRUE)) {
      j <- regexpr('=',a,fixed=TRUE)[1L]; key <- substring(a,3L,j-1L)
      v <- substring(a,j+1L)
    } else {
      key <- substring(a,3L)
      if (i == length(argv)) stop('Missing value: --',key)
      i <- i + 1L; v <- argv[[i]]
    }
    key <- gsub('-','_',key,fixed=TRUE)
    if (!key %in% names(o)) stop('Unknown option: --',key)
    o[[key]] <- v; i <- i + 1L
  }
  for (key in c('rds','chromosome_map','chromosome_lengths','catalog','helpers','template'))
    if(is.null(o[[key]])) stop('Missing --', key)
  if (!is.null(o$mane) && (is.null(o$mane_assembly) || o$mane_assembly != o$reference)) stop('MANE assembly must match the reference; renaming chromosomes is not liftover')
  o$min_depth <- as.numeric(o$min_depth); o$upper_depth_ratio <- as.numeric(o$upper_depth_ratio)
  if (!is.finite(o$min_depth) || o$min_depth < 0 || !is.finite(o$upper_depth_ratio) || o$upper_depth_ratio <= 0) stop('Invalid coverage limits')
  for (k in c('rds','mane','regions_bed','readme','chromosome_map','chromosome_lengths','catalog','helpers','template','panel_bed'))
    if (!is.null(o[[k]]) && !file.exists(o[[k]])) stop('File missing: ',o[[k]])
  o$promoter_bp <- as.integer(o$promoter_bp)
  if (is.na(o$promoter_bp)||o$promoter_bp<0L) stop('Invalid --promoter-bp')
  if (is.null(o$vcf)) o$vcf <- paste0(tools::file_path_sans_ext(o$out),'.vcf')
  if (is.null(o$tsv)) o$tsv <- paste0(tools::file_path_sans_ext(o$out),'.tsv')
  o
}
if (!requireNamespace('jsonlite',quietly=TRUE))
  stop("Install jsonlite first: install.packages('jsonlite')")
opts <- parse_args(commandArgs(trailingOnly=TRUE))
message('[ATLAS] Reading RDS: ',opts$rds)
input <- readRDS(opts$rds)
if (is.data.frame(input)) {
  loci <- input; sample_id <- 'sample'; format_version <- 'data.frame RDS'
} else if (is.list(input) && is.data.frame(input$loci)) {
  loci <- input$loci
  sample_id <- if (length(input$sample_id)) as.character(input$sample_id[[1L]]) else 'sample'
  format_version <- if (length(input$format_version)) as.character(input$format_version[[1L]]) else 'loci'
} else stop('RDS must contain a data.frame or a list with $loci')
if (!all(c('CHROM','POS') %in% names(loci))) stop('Expected CHROM and POS')

source(opts$helpers)
data.table::setDTthreads(1L)
mapping <- read_chromosome_map(opts$chromosome_map)
chromosome_lengths <- read.delim(opts$chromosome_lengths, stringsAsFactors=FALSE)
chroms <- chromosome_lengths$FINAL_CHR
if (!length(chroms) || anyDuplicated(chroms) || any(chromosome_lengths$length <= 0)) stop('Invalid chromosome lengths')
if (is.null(input$metadata) || input$metadata$assembly != opts$reference) stop('RDS reference mismatch')
if (!'FINAL_CHR' %in% names(loci) || any(!loci$FINAL_CHR %in% chroms)) stop('Expected validated primary FINAL_CHR')
if (input$metadata$sex == 'F' && any(loci$FINAL_CHR == 'chrY')) stop('Female sample contains chrY loci')
report_panel <- read_panel(opts$panel_bed)
if(is.null(report_panel) && isTRUE(input$metadata$adaptive_sampling))
  report_panel <- data.table::as.data.table(input$metadata$target_regions)
if(!all(panel_membership(loci,report_panel))) stop('Finalized loci outside report target panel')
# Recheck current thresholds even for older finalized bundles.
report_coverage <- if(!is.null(opts$regions_bed)) read_coverage(opts$regions_bed,mapping,
  opts$min_depth,opts$upper_depth_ratio,input$metadata$sex,report_panel) else NULL
report_primary <- primary_partition(loci,report_coverage,as.numeric(opts$primary_min_support),as.numeric(opts$min_svlen))
loci <- report_primary$loci
write_flat_tsv(report_primary$rescue_candidates,paste0(opts$out,'.excluded_primary.tsv'))
native_chroms <- mapping$CHROM[match(chroms, mapping$FINAL_CHR)]
if(anyNA(native_chroms) || anyDuplicated(native_chroms) ||
   any(!loci$CHROM %in% native_chroms)) stop('Invalid native VCF chromosome mapping')
chr_norm <- function(x) toupper(sub('^chr','',as.character(x),ignore.case=TRUE))
catalog <- read.delim(opts$catalog, stringsAsFactors=FALSE, check.names=FALSE)
aliases <- setNames(rep(catalog$FINAL_CHR, lengths(strsplit(catalog$aliases,';',fixed=TRUE))), unlist(strsplit(catalog$aliases,';',fixed=TRUE)))
annotation_available <- !is.null(opts$mane)

message('[ATLAS] Reading assembly-specific gene annotation')
mane <- if(annotation_available) read_gene_annotation(opts$mane) else
  as.data.frame(matrix(character(), nrow=0L, ncol=25L))
names(mane) <- c('chrom','chromStart','chromEnd','transcriptName','score','strand',
  'thickStart','thickEnd','itemRgb','blockCount','blockSizes','chromStarts',
  'transcriptID','cdsStartStat','cdsEndStat','exonFrames','featureType','geneID',
  'geneSymbol','geneType','proteinID','refseqTranscriptID','geneIDNCBI',
  'refseqProteinID','maneStatus')
original_mane_n <- nrow(mane)
mane$chrom <- unname(aliases[mane$chrom])
mane <- mane[!is.na(mane$chrom) & mane$chrom %in% chroms,,drop=FALSE]
if(annotation_available && !nrow(mane)) stop('MANE has no chromosomes matching this assembly')
if (nrow(mane)) {
  limits <- chromosome_lengths$length[match(mane$chrom, chroms)]
  if (any(!is.finite(as.numeric(mane$chromStart)) | as.numeric(mane$chromStart) < 0 |
          as.numeric(mane$chromEnd) > limits | as.numeric(mane$chromEnd) <= as.numeric(mane$chromStart))) stop('MANE BED coordinates outside reference')
}

mane$maneStatus <- trimws(mane$maneStatus)
mane <- mane[mane$maneStatus %in% c('MANE Select','MANE Plus Clinical','RefSeq Select Curated') &
               !is.na(mane$geneSymbol)&nzchar(mane$geneSymbol),,drop=FALSE]
mane$priority <- ifelse(mane$maneStatus=='MANE Select',0L,1L)
mane <- mane[order(mane$priority,chr_norm(mane$chrom),mane$geneSymbol,
                   as.numeric(mane$chromStart)),,drop=FALSE]
mane <- mane[!duplicated(paste(chr_norm(mane$chrom),mane$geneSymbol,sep='||')),,drop=FALSE]
rownames(mane) <- NULL
mane$start <- as.numeric(mane$chromStart)+1
mane$end <- as.numeric(mane$chromEnd)
mane$thickStart <- as.numeric(mane$thickStart)
mane$thickEnd <- as.numeric(mane$thickEnd)
mane$chrkey <- chr_norm(mane$chrom)
mane$promStart <- ifelse(mane$strand=='+',pmax(1,mane$start-opts$promoter_bp),mane$end+1)
mane$promEnd <- ifelse(mane$strand=='+',mane$start-1,mane$end+opts$promoter_bp)
chr_index <- split(seq_len(nrow(mane)),mane$chrkey)
parse_bed_int <- function(x) { y <- strsplit(as.character(x),',',fixed=TRUE)[[1L]]; as.numeric(y[nzchar(y)]) }
annotate_one <- function(chrom,pos) {
  if (!annotation_available) return(list(Gene=NA_character_,Location='unannotated',Region='Unannotated'))
  nohit <- list(Gene='NO',Location='intergenic',Region='Intergenic')
  p <- suppressWarnings(as.numeric(pos))
  if (is.na(chrom)||!is.finite(p)) return(nohit)
  ind <- chr_index[[chr_norm(chrom)]]
  if (is.null(ind)) return(nohit)
  in_gene <- !is.na(mane$start[ind]) & !is.na(mane$end[ind]) &
    p>=mane$start[ind] & p<=mane$end[ind]
  in_prom <- !is.na(mane$promStart[ind]) & !is.na(mane$promEnd[ind]) &
    p>=mane$promStart[ind] & p<=mane$promEnd[ind]
  candidates <- ind[which(in_gene|in_prom)]
  if (!length(candidates)) return(nohit)
  best <- nohit; best_score <- 0L
  for(j in candidates) {
    tx <- mane[j,,drop=FALSE]; strand <- tx$strand[[1L]]
    if(is.na(strand)||!strand %in% c('+','-'))next
    score<-0L; gene<-tx$geneSymbol[[1L]];location<-'intergenic';region<-'Intergenic'
    if(p<tx$start[[1L]]||p>tx$end[[1L]]){
      if(p>=tx$promStart[[1L]]&&p<=tx$promEnd[[1L]]){
        score<-1L;location<-'promoter';region<-'Promoter'
      }
    } else {
      offsets<-parse_bed_int(tx$chromStarts[[1L]])
      sizes<-parse_bed_int(tx$blockSizes[[1L]])
      expected<-suppressWarnings(as.integer(tx$blockCount[[1L]]))
      if(is.na(expected)||length(offsets)!=expected||length(sizes)!=expected||
         anyNA(offsets)||anyNA(sizes))next
      xs<-as.numeric(tx$chromStart[[1L]])+offsets+1
      xe<-as.numeric(tx$chromStart[[1L]])+offsets+sizes
      oo<-order(xs,decreasing=identical(strand,'-'));xs<-xs[oo];xe<-xe[oo]
      exon_hit<-which(p>=xs & p<=xe)
      if(length(exon_hit)){
        k<-exon_hit[[1L]];location<-paste0('exon',k)
        region<-'Exon'
        score<-4L
      }else if(length(xs)>1L){
        for(k in seq_len(length(xs)-1L)){
          intr_start<-min(xe[k],xe[k+1L])+1
          intr_end<-max(xs[k],xs[k+1L])-1
          if(p>=intr_start&&p<=intr_end){
            location<-paste0('intron',k)
            region<-'Intron'
            score<-2L
            break
          }
        }
      }
    }
    if(score>best_score){best<-list(Gene=gene,Location=location,Region=region);best_score<-score}
    if(best_score==4L)break
  }
  best
}
message('[ATLAS] Annotating ',nrow(loci),' loci')
ann <- lapply(seq_len(nrow(loci)),function(i)annotate_one(loci$FINAL_CHR[[i]],loci$POS[[i]]))
ann_gene<-vapply(ann,`[[`,character(1),'Gene');ann_loc<-vapply(ann,`[[`,character(1),'Location')
ann_reg<-vapply(ann,`[[`,character(1),'Region')
chrom_key<-loci$FINAL_CHR; ranks<-match(chrom_key,chroms)
ranks[is.na(ranks)]<-length(chroms)+1L
positions<-suppressWarnings(as.numeric(loci$POS))
oo<-order(ranks,chrom_key,-positions,na.last=TRUE)
loci<-loci[oo,,drop=FALSE];chrom_key<-chrom_key[oo];positions<-positions[oo]
ann_gene<-ann_gene[oo];ann_loc<-ann_loc[oo];ann_reg<-ann_reg[oo]
n<-nrow(loci)
field<-function(k){if(k %in% names(loci))loci[[k]] else rep(NA,n)}
num<-function(k) { v<-suppressWarnings(as.numeric(field(k)));v[!is.finite(v)]<-NA_real_;v }
chars<-function(k){ v<-as.character(field(k));v[is.na(v)|!nzchar(trimws(v))]<-NA_character_;v }
flag<-function(k){toupper(trimws(chars(k))) %in% c('TRUE','T','1','YES')}
polylen<-num('polyA_length');polylen[is.na(polylen)|polylen<=0]<-NA_real_
tsdlen<-num('tsd_length');tsdlen[is.na(tsdlen)|tsdlen<=0]<-NA_real_
polyflag<-flag('polyA');bio <- polyflag & !is.na(polylen) & !is.na(tsdlen)
status1<-chars('cover_status');status2<-chars('region_status')
callstatus<-ifelse(is.na(status1)|is.na(status2),'NA',ifelse(status1=='PASS'&status2=='PASS','PASS','NO'))
source_ref<-chars('REF');source_alt<-chars('ALT')
fragq<-num('FRAGMENT_Q');fragq[!is.na(fragq)]<-round(fragq[!is.na(fragq)])
rows<-data.frame(
  ID=as.character(loci$IDs),CHROM=chrom_key,FINAL_CHR=chrom_key,ORIGINAL_CHROM=as.character(loci$CHROM),POS=positions,
  ASSEMBLY=chars('ASSEMBLY'),SEX=chars('SEX'),PROCEDENCE=chars('PROCEDENCE'),
  Gene=ann_gene,Location=ann_loc,Region=ann_reg,
  SVLEN=num('svlen'),GT=chars('GT'),Class=chars('event_class'),
  Elements=chars('repeat_elements'),Event=chars('event'),
  Completeness=chars('completeness'),PolyA=as.logical(polyflag),
  PolyA_length=polylen,TSD=tsdlen,TSD_sequence=chars('tsd_sequence'),
  BIO=as.logical(bio),Fragment_Q=fragq,VCF_QUAL=chars('QUAL'),
  VCF_FILTER=chars('FILTER'),Cover_status=status1,Region_status=status2,
  Call_status=callstatus,
  PolyT5=num('length_poly_T_5'),PolyT5_purity=num('purity_poly_T_5'),
  PolyA3=num('length_poly_A_3'),PolyA3_purity=num('purity_poly_A_3'),
  Gap_min=num('FRAGMENT_MIN_GAP_BP'),Gap_max=num('FRAGMENT_MAX_GAP_BP'),
  Fragment_span=num('FRAGMENT_COMPONENT_SPAN_BP'),
  Event_coverage=num('event_coverage_pct'),Model_coverage=num('model_coverage_pct'),
  RepeatMasker_coverage=num('repeatmasker_coverage_pct'),
  REF=source_ref,ALT=source_alt,stringsAsFactors=FALSE,check.names=FALSE)

# Callability uses the exact same regions-only policy as reconstruction.
callability <- list()
if(!is.null(opts$regions_bed)) {
  callability <- coverage_display(report_coverage)
}
meta<-list(chromosomes=as.list(chroms), chromosomeLengths=as.list(chromosome_lengths$length),
  vcfChromosomes=as.list(native_chroms), vcfChromosomeLengths=as.list(chromosome_lengths$length),
  sex=input$metadata$sex, procedence=input$metadata$procedence, annotationAvailable=annotation_available,
  warnings=as.list(c(input$metadata$reference_validation$errors, input$metadata$reference_validation$warnings,
    if(!annotation_available) 'Gene annotation unavailable: no MANE BED in this assembly. Coordinates have not been lifted over.')),
  sample=sample_id,version=format_version,rdsFile=basename(opts$rds),
  maneFile=if(annotation_available) basename(opts$mane) else 'Not supplied for this assembly',promoterBp=opts$promoter_bp,
  maneOriginal=original_mane_n,maneSelect=sum(mane$maneStatus=='MANE Select'),
  manePlus=sum(mane$maneStatus=='MANE Plus Clinical'),refseqSelect=sum(mane$maneStatus=='RefSeq Select Curated'),reference=opts$reference,
  coverageFile=if(is.null(opts$regions_bed)) 'Not supplied' else basename(opts$regions_bed))
meta$runConfiguration <- jsonlite::fromJSON(opts$run_config,simplifyVector=FALSE)
meta$rdsParameters <- input$parameters
meta$reportParameters <- list(primary_min_support=as.numeric(opts$primary_min_support),
  min_svlen=as.numeric(opts$min_svlen),min_depth=opts$min_depth,
  upper_depth_ratio=opts$upper_depth_ratio,promoter_bp=opts$promoter_bp,reference=opts$reference)
extra_readme<-if(is.null(opts$readme)) '' else
  paste(readLines(opts$readme,warn=FALSE,encoding='UTF-8'),collapse='\n')

# VCF 4.3. GT is written in FORMAT, not INFO.
# TSD is missing (.) if no TSD sequence evaluation result is available.
vcf_safe <- function(x) {
  if(length(x)==0L||is.na(x)||!nzchar(as.character(x)))return('.')
  x<-as.character(x)
  x<-gsub('%','%25',x,fixed=TRUE)
  x<-gsub(';','%3B',x,fixed=TRUE);x<-gsub(',','%2C',x,fixed=TRUE)
  x<-gsub('=','%3D',x,fixed=TRUE);x<-gsub(' ','%20',x,fixed=TRUE)
  x<-gsub('\t','%09',x,fixed=TRUE);x<-gsub('\n','%0A',x,fixed=TRUE)
  x
}
info_defs<-list(
  c('ORIGINAL_CHROM','1','String','Original reference chromosome identifier'),
  c('FINAL_CHR','1','String','Canonical chromosome label used in the report'),
  c('ASSEMBLY','1','String','Reference assembly'),
  c('SEX','1','String','Sample sex M F or U'),
  c('PROCEDENCE','1','String','Sample provenance'),
  c('SVTYPE','1','String','Structural variant type'),
  c('SVLEN','1','Integer','Length of inserted ALT sequence or summary length for symbolic INS'),
  c('RPT_SVLEN','1','Float','Original locus-level summary insertion length'),
  c('CLASS','1','String','Repeat element class'),
  c('ELEMENTS','1','String','RepeatMasker element annotation'),
  c('EVENT','1','String','Structural event architecture'),
  c('COMPLETENESS','1','String','Repeat model completeness'),
  c('POLYA','1','Integer','PolyA or polyT signal detected: 1 yes, 0 no'),
  c('POLYA_LEN','1','Integer','PolyA or polyT tail length'),
  c('TSD','1','Integer','TSD presence: 1 detected, 0 assessed absent, . unavailable'),
  c('TSD_LEN','1','Integer','Target-site duplication length in bp'),
  c('BIO','1','Integer','PolyA plus TSD double biological evidence: 1 yes, 0 no'),
  c('FRAGMENT_Q','1','Integer','Fragmentation quality score from 1 to 5'),
  c('EVENT_COV','1','Float','Event sequence coverage percent'),
  c('MODEL_COV','1','Float','Reference repeat model coverage percent'),
  c('RM_COV','1','Float','RepeatMasker query coverage percent'))
vcf_header<-c('##fileformat=VCFv4.3','##source=ATLAS_intrasample',
  paste0('##reference=',opts$reference),
  paste0('##contig=<ID=',native_chroms,',length=',chromosome_lengths$length,'>'),
  '##ALT=<ID=INS,Description="Insertion with unresolved alternate sequence">',
  '##FILTER=<ID=PASS,Description="Input call passed its original filter">',
  vapply(info_defs,function(x)paste0('##INFO=<ID=',x[1L],',Number=',x[2L],
    ',Type=',x[3L],',Description="',x[4L],'">'),character(1)),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste(c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',
    gsub('[^A-Za-z0-9_.-]','_',sample_id)),collapse='\t'))
make_vcf_row<-function(i){
  r<-rows[i,,drop=FALSE]
  ref<-if(!is.na(r$REF)&&grepl('^[ACGTN]$',r$REF,ignore.case=TRUE))toupper(r$REF) else 'N'
  sequence_alt<-if(!is.na(r$ALT) && grepl('^[ACGTN]+$',r$ALT,ignore.case=TRUE) &&
                  nchar(r$ALT)>nchar(ref) && startsWith(toupper(r$ALT),ref)) toupper(r$ALT) else '<INS>'
  insertion_length<-if(sequence_alt!='<INS>') nchar(sequence_alt)-nchar(ref) else
    if(is.na(r$SVLEN))NA else floor(r$SVLEN+0.5)
  raw_info<-c(ORIGINAL_CHROM=r$ORIGINAL_CHROM,FINAL_CHR=r$CHROM,ASSEMBLY=r$ASSEMBLY,SEX=r$SEX,PROCEDENCE=r$PROCEDENCE,SVTYPE='INS',SVLEN=insertion_length,RPT_SVLEN=r$SVLEN,CLASS=r$Class,ELEMENTS=r$Elements,
    EVENT=r$Event,COMPLETENESS=r$Completeness,POLYA=as.integer(r$PolyA),
    POLYA_LEN=r$PolyA_length,TSD=if(is.na(r$TSD))NA else 1L,
    TSD_LEN=r$TSD,BIO=as.integer(r$BIO),FRAGMENT_Q=r$Fragment_Q,
    EVENT_COV=r$Event_coverage,MODEL_COV=r$Model_coverage,
    RM_COV=r$RepeatMasker_coverage)
  full<-paste(paste(names(raw_info)[!is.na(raw_info)],
                  vapply(raw_info[!is.na(raw_info)],vcf_safe,character(1)),sep='='),collapse=';')
  qual<-if(!is.na(r$VCF_QUAL)&&grepl('^[0-9]+(\\.[0-9]+)?$',r$VCF_QUAL))r$VCF_QUAL else '.'
  gt<-if(!is.na(r$GT)&&grepl('^[0-9.]+[/|][0-9.]+$',r$GT))r$GT else './.'
  paste(c(r$ORIGINAL_CHROM,as.integer(r$POS),r$ID,ref,sequence_alt,qual,
    ifelse(is.na(r$VCF_FILTER),'PASS',r$VCF_FILTER),
    if(nzchar(full))full else '.','GT',gt),collapse='\t')
}
dir.create(dirname(opts$tsv),showWarnings=FALSE,recursive=TRUE)
write_flat_tsv(rows, opts$tsv)
dir.create(dirname(opts$vcf),showWarnings=FALSE,recursive=TRUE)
message('[ATLAS] Writing VCF: ',opts$vcf)
con<-file(opts$vcf,open='wt');writeLines(vcf_header,con,useBytes=TRUE)
for (i in seq_len(n)) writeLines(make_vcf_row(i),con,useBytes=TRUE)
close(con)

encode<-function(x,...){
  ans<-as.character(jsonlite::toJSON(x,auto_unbox=TRUE,na='null',null='null',digits=NA,...))
  gsub('<','\\u003c',ans,fixed=TRUE)
}
replace_token<-function(x,token,value){if(!grepl(token,x,fixed=TRUE))stop('Missing HTML token: ',token)
  paste(strsplit(x,token,fixed=TRUE)[[1L]],collapse=value)}

# Editable, standalone template; reference metadata is injected below.
html_template <- paste(readLines(opts$template, warn=FALSE), collapse='\n')
html <- replace_token(html_template,'@@ATLAS_DATA@@',encode(rows,dataframe='rows'))
html <- replace_token(html,'@@ATLAS_META@@',encode(meta))
html <- replace_token(html,'@@ATLAS_CALL@@',encode(callability))
html <- replace_token(html,'@@ATLAS_README@@',encode(extra_readme))
dir.create(dirname(opts$out),recursive=TRUE,showWarnings=FALSE)
writeLines(html,con=opts$out,useBytes=TRUE)
message('[ATLAS] HTML report: ',normalizePath(opts$out,winslash='/',mustWork=TRUE))
message('[ATLAS] VCF: ',normalizePath(opts$vcf,winslash='/',mustWork=TRUE))
message('[ATLAS] Loci: ',n,' | BIO: ',sum(rows$BIO,na.rm=TRUE))
