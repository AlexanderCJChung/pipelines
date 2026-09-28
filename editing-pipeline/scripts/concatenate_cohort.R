#!/usr/bin/env Rscript
# Concatenate the per-chromosome filtered matrices and emit ANNOVAR input.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo).
#
# Upstream had `rat_files` holding the COVERAGE paths and `cov_files` holding the RATIO
# paths, then swapped them back at the read - functionally correct, actively misleading.
# It also hardcoded the chromosome list here as well as in the Snakefile, so the two
# could drift. The list now arrives as --chroms from config, and the names say what
# they hold.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inDir'), help = 'directory holding {chrom}_{coverage,ratio}.tsv.gz', default = ''),
    make_option(c('--chroms'), help = 'comma-separated chromosome list', default = ''),
    make_option(c('--outCov'), help = 'output concatenated coverage matrix', default = ''),
    make_option(c('--outRat'), help = 'output concatenated ratio matrix', default = ''),
    make_option(c('--outAvinput'), help = 'output ANNOVAR avinput', default = ''),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
setDTthreads(opt$threads)

chroms <- trimws(strsplit(opt$chroms, ",", fixed = TRUE)[[1]])
chroms <- chroms[nzchar(chroms)]
stopifnot(length(chroms) > 0)

read_all <- function(suffix) {
    paths <- file.path(opt$inDir, paste0(chroms, "_", suffix, ".tsv.gz"))
    missing <- !file.exists(paths)
    if (any(missing)) stop("missing input: ", paste(paths[missing], collapse = ", "))
    # A chromosome can legitimately contribute zero surviving sites; rbind with
    # use.names guards against any column-order difference between shards.
    rbindlist(lapply(paths, fread, sep = "\t", showProgress = FALSE), use.names = TRUE)
}

message(" * reading coverage shards")
all_cov <- read_all("coverage")
message(" * reading ratio shards")
all_rat <- read_all("ratio")

stopifnot(identical(all_cov$ESid, all_rat$ESid))
message(" * ", nrow(all_cov), " sites concatenated across ", length(chroms), " chromosomes")

message(" * writing coverage")
fwrite(all_cov, opt$outCov, sep = "\t", compress = "gzip")
message(" * writing ratios")
fwrite(all_rat, opt$outRat, sep = "\t", compress = "gzip")

# ANNOVAR avinput: Chr Start End Ref Alt, 1-based and inclusive. Editing sites are
# SNVs, so Start == End == the ESid position.
message(" * writing ANNOVAR avinput")
av <- data.table(ESid = all_cov$ESid)
av[, c("Chr", "Pos", "Ref", "Alt") := tstrsplit(ESid, ":", fixed = TRUE)]
av <- av[, .(Chr, Start = Pos, End = Pos, Ref, Alt)]
fwrite(av, opt$outAvinput, sep = "\t", col.names = FALSE)
message(" * done")
