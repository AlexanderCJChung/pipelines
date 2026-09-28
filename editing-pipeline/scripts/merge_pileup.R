#!/usr/bin/env Rscript
# Merge the per-sample pileup back-fill into the final cohort matrices.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo), rewritten.
#
# ============================ DELIBERATE DIVERGENCE FROM UPSTREAM ============================
# Allele selection happens HERE, against the annotation, rather than per-sample in
# parse_pileup.R. See the header of parse_pileup.R for the full argument. The short
# version: the pileup runs against a fixed site list whose alternative allele the cohort
# has already decided, so every sample must be asked for the SAME allele at the SAME
# locus. Doing it per-sample (upstream) means a sample with zero edited reads emits no
# row and becomes NA - a true zero turned into missing data - and samples whose
# top-covered alt differs from the cohort's fail to join at all.
#
# Consequence: a locus a sample covered at >= --minCov reads but did not edit is
# recorded as total_cov = n, edit_rate = 0. NA now means only "not enough coverage to
# say", which is what the missingness filters below assume.
# ============================================================================================
#
# Also fixed: upstream assigned the flipped ESids positionally
#   anno_DF <- anno_df[anno_df$ESid %in% coverage_DF$ESid,]
#   coverage_DF$ESid <- anno_DF$ESid2
# which is correct only while both frames happen to share row order. Here it is a match().

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inDir'), help = 'project directory holding {sample}/{subDir}/{sample}.pileup.parsed.txt', default = ''),
    make_option(c('--subDir'), help = 'per-tool subdirectory under {sample}/ (may be empty)', default = ''),
    make_option(c('--metadata'), help = 'TSV with a `sample` column', default = ''),
    make_option(c('--anno'), help = 'editing_annotation.txt from filter_annovar.R', default = ''),
    make_option(c('--outCov'), help = 'output coverage matrix', default = ''),
    make_option(c('--outRat'), help = 'output editing-rate matrix', default = ''),
    make_option(c('--outAnno'), help = 'output annotation', default = ''),
    make_option(c('--outDtu'), help = 'output ref/alt count matrix for DTU-style testing', default = ''),
    make_option(c('--minCov'), type = 'integer', default = 5L,
                help = 'reads required at a locus before an editing rate is trusted; below this the cell is NA [default %default]'),
    make_option(c('--missing_sites'), type = 'double', default = 0.75,
                help = 'keep sites with non-NA coverage in >= this proportion of samples [default %default]'),
    make_option(c('--missing_samples'), type = 'double', default = 0.5,
                help = 'keep samples with non-NA coverage at >= this proportion of sites [default %default]'),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
setDTthreads(opt$threads)

BASES <- c("A", "C", "G", "T")

meta <- fread(opt$metadata, sep = "\t")
stopifnot("sample" %in% names(meta))
sample_ids <- as.character(meta$sample)

anno <- fread(opt$anno, sep = "\t", showProgress = FALSE, colClasses = list(character = c("Chr", "Ref", "Alt", "ESid")))
stopifnot(!anyDuplicated(anno$ESid))
message(" * ", nrow(anno), " annotated sites")

anno[, chrpos := paste0(Chr, ":", Start)]
ref_col <- match(anno$Ref, BASES)
alt_col <- match(anno$Alt, BASES)
if (anyNA(ref_col) || anyNA(alt_col)) {
    stop("annotation contains non-ACGT alleles at ", sum(is.na(ref_col) | is.na(alt_col)), " site(s)")
}

n_sites <- nrow(anno)
n_samp  <- length(sample_ids)

cov_m <- matrix(NA_integer_, n_sites, n_samp)
ref_m <- matrix(NA_integer_, n_sites, n_samp)
alt_m <- matrix(NA_integer_, n_sites, n_samp)

message(" * reading ", n_samp, " pileup files")
for (i in seq_len(n_samp)) {
    sdir <- if (nzchar(opt$subDir)) file.path(opt$inDir, sample_ids[i], opt$subDir) else file.path(opt$inDir, sample_ids[i])
    f <- file.path(sdir, paste0(sample_ids[i], ".pileup.parsed.txt"))
    if (!file.exists(f)) stop("missing pileup for sample ", sample_ids[i], ": ", f)
    pu <- fread(f, sep = "\t", showProgress = FALSE, colClasses = list(character = "chrpos"))
    if (nrow(pu) == 0L) next
    idx <- match(anno$chrpos, pu$chrpos)
    base_m <- as.matrix(pu[, ..BASES])
    cov_m[, i] <- pu$total_cov[idx]
    ref_m[, i] <- base_m[cbind(idx, ref_col)]
    alt_m[, i] <- base_m[cbind(idx, alt_col)]
    if (i %% 50L == 0L) message("   - ", i, "/", n_samp)
}

# Below minCov the locus is not measured well enough to distinguish "unedited" from
# "undersampled", so it is NA everywhere. At or above minCov, alt_cov = 0 is a real zero.
under <- is.na(cov_m) | cov_m < opt$minCov
cov_m[under] <- NA_integer_
ref_m[under] <- NA_integer_
alt_m[under] <- NA_integer_

rat_m <- alt_m / cov_m
message(" * ", format(sum(!is.na(rat_m))), " observed site x sample cells; of those ",
        format(sum(rat_m == 0, na.rm = TRUE)), " are true zeros")

# --- missingness -------------------------------------------------------------------
observed <- !is.na(cov_m)
site_keep <- rowSums(observed) >= opt$missing_sites * n_samp
message(" * keeping ", sum(site_keep), " of ", n_sites, " sites with coverage in >= ",
        opt$missing_sites * 100, "% of samples")
if (!any(site_keep)) stop("no sites survived the site-missingness filter")

cov_m <- cov_m[site_keep, , drop = FALSE]
ref_m <- ref_m[site_keep, , drop = FALSE]
alt_m <- alt_m[site_keep, , drop = FALSE]
rat_m <- rat_m[site_keep, , drop = FALSE]
anno  <- anno[site_keep]

samp_keep <- colSums(!is.na(cov_m)) >= opt$missing_samples * nrow(cov_m)
message(" * keeping ", sum(samp_keep), " of ", n_samp, " samples with coverage at >= ",
        opt$missing_samples * 100, "% of sites")
if (!any(samp_keep)) stop("no samples survived the sample-missingness filter")

cov_m <- cov_m[, samp_keep, drop = FALSE]
ref_m <- ref_m[, samp_keep, drop = FALSE]
alt_m <- alt_m[, samp_keep, drop = FALSE]
rat_m <- rat_m[, samp_keep, drop = FALSE]
kept_samples <- sample_ids[samp_keep]

# --- orient alleles to the gene ----------------------------------------------------
# JACUSA reports everything on the reference strand. For a minus-strand gene the
# biological A>I edit appears as T>C, so complement the alleles. The COUNTS are
# unchanged - they are counts of the same reads, relabelled to the gene's strand.
message(" * flipping alleles on minus-strand genes")
comp <- c(A = "T", C = "G", G = "C", T = "A")
flip <- !is.na(anno$strand) & anno$strand == "-"
anno[, Ref_gene := ifelse(flip, unname(comp[Ref]), Ref)]
anno[, Alt_gene := ifelse(flip, unname(comp[Alt]), Alt)]
anno[, ESid_gene := paste0(Chr, ":", Start, ":", Ref_gene, ":", Alt_gene)]
message("   - ", sum(flip), " of ", nrow(anno), " sites flipped")

esid <- anno$ESid_gene

write_matrix <- function(m, path) {
    out <- as.data.table(m)
    setnames(out, kept_samples)
    out[, ESid := esid]
    setcolorder(out, "ESid")
    fwrite(out, path, sep = "\t", compress = "gzip")
}

message(" * writing matrices")
write_matrix(cov_m, opt$outCov)
write_matrix(rat_m, opt$outRat)

# DTU-style matrix: two rows per site (reference and edited allele counts), which is
# the shape DRIMSeq/satuRn style tests expect.
dtu <- rbind(
    data.table(ESid = esid, allele = paste0(esid, ":ref"), as.data.table(ref_m)),
    data.table(ESid = esid, allele = paste0(esid, ":alt"), as.data.table(alt_m))
)
setnames(dtu, c("ESid", "allele", kept_samples))
setorder(dtu, ESid, allele)
fwrite(dtu, opt$outDtu, sep = "\t", compress = "gzip")

anno_out <- anno[, .SD, .SDcols = setdiff(names(anno), c("chrpos", "ESid"))]
setnames(anno_out, "ESid_gene", "ESid")
setcolorder(anno_out, "ESid")
fwrite(anno_out, opt$outAnno, sep = "\t", compress = "gzip")
message(" * done")
