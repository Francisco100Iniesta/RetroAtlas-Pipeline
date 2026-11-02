library(Biostrings)
library(dplyr)
library(stringr)

log_connection <- file(snakemake@log[[1]], "wt")
sink(log_connection)
sink(log_connection, type = "message")
min_svlen <- as.integer(snakemake@params[["min_svlen"]])
stopifnot(!is.na(min_svlen), min_svlen > 0)

vcf_names <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
               "INFO", "FORMAT", "GENOTYPE")
data <- as.data.frame(setNames(rep(list(character()), 10), vcf_names))
data$POS <- integer()
data$svtype <- character()
data$svlen <- numeric()
data$ID2 <- integer()
data$header <- character()

files <- sort(list.files(snakemake@input[["results"]], pattern = "\\.vcf$", full.names = TRUE))
for (vcf in files) {
    lines <- readLines(vcf, warn = FALSE)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    x <- read.table(text = paste(lines, collapse = "\n"), sep = "\t",
                    header = FALSE, quote = "", comment.char = "", colClasses = "character")
    if (ncol(x) != 10) stop("Se esperan 10 columnas en ", vcf)
    names(x) <- vcf_names
    x$POS <- as.integer(x$POS)
    x$svtype <- str_match(x$INFO, "(?:^|;)SVTYPE=([^;]+)")[, 2]
    x$svlen <- suppressWarnings(as.numeric(str_match(x$INFO, "(?:^|;)SVLEN=([^;]+)")[, 2]))
    x <- filter(x, svtype == "INS", svlen >= min_svlen)
    if (!nrow(x)) next
    # Solo ALT resueltos como secuencia pueden convertirse a DNAStringSet.
    if (any(!grepl("^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$", x$ALT))) {
        stop("ALT no resuelto como secuencia en ", vcf)
    }
    x$ID <- tools::file_path_sans_ext(basename(vcf))
    x$ID2 <- seq_len(nrow(x))
    x$header <- paste(x$ID, x$ID2, sep = ";")
    data <- bind_rows(data, x)
}

fasta <- DNAStringSet(data$ALT)
names(fasta) <- data$header
if (length(fasta)) {
    writeXStringSet(fasta, snakemake@output[["fasta"]], format = "fasta")
} else {
    file.create(snakemake@output[["fasta"]])
}
saveRDS(data, snakemake@output[["rds"]])
message("Inserciones exportadas: ", nrow(data), "; min_svlen: ", min_svlen)
sink(type = "message")
sink()
close(log_connection)
