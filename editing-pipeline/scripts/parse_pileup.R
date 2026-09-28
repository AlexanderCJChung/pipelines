#!/usr/bin/env Rscript
# Parse a JACUSA2 pileup (call-1 -A -b sites.bed) into one row per LOCUS.
#
# ============================ DELIBERATE DIVERGENCE FROM UPSTREAM ============================
# Upstream (RajLabMSSM/editing-pipeline) resolves the alternative allele here, in the
# per-sample parse, exactly as the discovery parser does. That is wrong for the pileup
# pass, in two compounding ways:
#
#   1. Commit 91d3898 ("manuscript version") added `alt_cov > 0` to the filter. A sample
#      with 200 clean reads and zero edited reads then emits NO row for the locus, so
#      merge_pileup's left_join yields NA - indistinguishable from a sample that had no
#      coverage at all. True zeros become missing data. Editing rates are biased upward
#      and the missingness filters discard well-measured unedited samples.
#      The pre-91d3898 version knew this: "rows where alt_cov=0 are kept here but
#      considered by missingness filter later on".
#
#   2. Even with zeros kept, per-sample allele resolution is the wrong operation. The
#      pileup runs against a FIXED site list whose alternative allele is already decided
#      by the cohort discovery pass. If sample X's top alt at chr1:100 is G and sample Y's
#      is C, upstream builds different ESids for the same locus and the join drops both.
#      Restoring the old `alt_cov > 0` behaviour alone does not fix this - with all-zero
#      alt counts the "top" allele is whichever sorts first, i.e. arbitrary.
#
# So this script does no allele selection and no filtering at all. It emits the full
# A/C/G/T count vector per locus and lets merge_pileup.R read off the count for the
# allele the annotation specifies. Zero is then a real, observed zero, and every sample
# is queried for the same allele at the same locus.
# ============================================================================================
#
# JACUSA2 BED6-extended output: line 1 is a track line, line 2 is the header
#   #contig start end name score strand bases11 info filter ref

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--input'),  help = 'JACUSA2 pileup result file', default = ''),
    make_option(c('--output'), help = 'output TSV', default = ''),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
setDTthreads(opt$threads)

if (!file.exists(opt$input)) stop("File ", opt$input, " does not exist")
message("Loading JACUSA2 pileup from ", opt$input)

dt <- fread(opt$input, skip = 1L, header = TRUE, showProgress = FALSE)

out_cols <- c("chrpos", "ref", "A", "C", "G", "T", "total_cov")
if (nrow(dt) == 0L || !"bases11" %in% names(dt)) {
    empty <- setNames(
        data.table(character(), character(), integer(), integer(), integer(),
                   integer(), integer()),
        out_cols)
    fwrite(empty, opt$output, sep = "\t")
    message(" * 0 loci in input - wrote empty output")
    quit(save = "no", status = 0)
}

dt[, c("A", "C", "G", "T") := tstrsplit(bases11, ",", fixed = TRUE, type.convert = TRUE)]
setnames(dt, "#contig", "chr")
dt[, chrpos := paste0(chr, ":", end)]
dt[, total_cov := A + C + G + T]

# A locus can appear once per BED interval; the site list is unique by construction,
# but guard anyway so the downstream match() is unambiguous.
dt <- dt[!duplicated(chrpos)]

message(" * ", nrow(dt), " loci parsed")
fwrite(dt[, ..out_cols], opt$output, sep = "\t")
