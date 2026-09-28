#!/usr/bin/env Rscript
# Merge per-sample discovery calls for ONE chromosome into sites x samples matrices.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo).
#
# Two changes from upstream, both needed to reach 583 samples:
#
#   * Upstream builds the matrices with `reduce(cbind)` over N per-sample data.frames,
#     which reallocates the whole matrix once per sample - O(N^2) copying, and the
#     reason the rule was provisioned 100 GB of himem. Here the long table is filled
#     into a preallocated matrix by (row, col) index in one pass.
#
#   * Upstream discovers samples by globbing the output directory. That silently picks
#     up stale files from earlier runs and makes column order depend on the filesystem.
#     The sample list and its order now come from the metadata sheet, so every
#     chromosome produces identically-ordered columns and concatenate_cohort can rbind
#     them safely.
#
# The awk prefilter is kept from upstream: each sample's .filtAll is genome-wide, and
# one Snakemake job runs per chromosome, so pushing the chromosome filter into awk
# avoids reading every whole file 24 times.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inDir'), help = 'project directory holding {sample}/{subDir}/{sample}.filtAll', default = '.'),
    make_option(c('--subDir'), help = 'per-tool subdirectory under {sample}/ (may be empty)', default = ''),
    make_option(c('--metadata'), help = 'TSV with a `sample` column', default = ''),
    make_option(c('--chr'), help = 'chromosome to work on', default = ''),
    make_option(c('--outCov'), help = 'output coverage matrix (.tsv.gz)', default = ''),
    make_option(c('--outRat'), help = 'output editing-rate matrix (.tsv.gz)', default = ''),
    make_option(c('--minSamples'), type = 'integer', default = 2L,
                help = 'drop sites seen in fewer than this many samples [default %default]'),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
setDTthreads(opt$threads)

chrom <- opt$chr
meta <- fread(opt$metadata, sep = "\t")
stopifnot("sample" %in% names(meta))
sample_ids <- as.character(meta$sample)

sample_dir <- if (nzchar(opt$subDir)) file.path(opt$inDir, sample_ids, opt$subDir) else file.path(opt$inDir, sample_ids)
files <- file.path(sample_dir, paste0(sample_ids, ".filtAll"))
missing <- !file.exists(files)
if (any(missing)) {
    stop("missing .filtAll for ", sum(missing), " sample(s), e.g. ",
         paste(head(basename(files[missing]), 3), collapse = ", "))
}
message(" * ", length(files), " samples, chromosome ", chrom)

# Column 2 of .filtAll is `chr` (see parse_jacusa.R). NR==1 keeps the header.
read_chrom <- function(f) {
    cmd <- sprintf("awk 'NR == 1 || $2 == \"%s\"' %s", chrom, shQuote(f))
    fread(cmd = cmd, showProgress = FALSE,
          select = c("ESid", "total_cov", "edit_rate"),
          colClasses = list(character = "ESid"))
}

message(" * reading files")
dat <- rbindlist(lapply(seq_along(files), function(i) {
    d <- read_chrom(files[i])
    if (nrow(d)) d[, sample_idx := i]
    d
}), fill = TRUE)

if (nrow(dat) == 0L) {
    message(" * no sites on ", chrom, " in any sample - writing empty matrices")
    empty <- data.table(ESid = character())
    for (s in sample_ids) empty[, (s) := numeric()]
    fwrite(empty, opt$outCov, sep = "\t", compress = "gzip")
    fwrite(empty, opt$outRat, sep = "\t", compress = "gzip")
    quit(save = "no", status = 0)
}

# Drop singletons - a site seen in only one sample is not analysable across a cohort
# and dominates the matrix row count. Upstream hardcoded this as `duplicated()`, i.e.
# minSamples = 2; it is a parameter here.
n_per_site <- dat[, .N, by = ESid]
keep <- n_per_site[N >= opt$minSamples, ESid]
message(" * ", nrow(n_per_site), " sites seen; ", length(keep),
        " in >= ", opt$minSamples, " samples")
dat <- dat[ESid %chin% keep]

sites <- sort(unique(dat$ESid))
dat[, row_idx := match(ESid, sites)]

message(" * filling ", length(sites), " x ", length(sample_ids), " matrices")
idx <- cbind(dat$row_idx, dat$sample_idx)

cov_mat <- matrix(NA_integer_, nrow = length(sites), ncol = length(sample_ids))
cov_mat[idx] <- as.integer(dat$total_cov)

rat_mat <- matrix(NA_real_, nrow = length(sites), ncol = length(sample_ids))
rat_mat[idx] <- as.numeric(dat$edit_rate)

write_matrix <- function(m, path) {
    out <- as.data.table(m)
    setnames(out, sample_ids)
    out[, ESid := sites]
    setcolorder(out, "ESid")
    fwrite(out, path, sep = "\t", compress = "gzip")
}

write_matrix(cov_mat, opt$outCov)
write_matrix(rat_mat, opt$outRat)
message(" * done")
