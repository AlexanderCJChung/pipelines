#!/usr/bin/env Rscript
# Apply cohort-level filters to one chromosome's sites x samples matrices.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo), rewritten.
# Three fixes, all documented in the README:
#
#   * ESids are no longer round-tripped through make.names. Upstream transposed the
#     matrix with `data.frame(t(...))`, which mangles "chr1:12345:A:G" into
#     "chr1.12345.A.G", then converted back with gsub("\\.", ":"). That is lossy for
#     any contig whose name contains a dot, and it cost two full transposes of a
#     sites x samples matrix. Here the matrix is filtered in place with rowSums/rowMeans
#     over group column indices; ESids are never used as R names.
#
#   * The metadata schema is no longer hardcoded. Upstream did
#     `select(group_df, -c(sample, bam_path, library))` and treated whatever remained
#     as the grouping variable, so any extra column silently corrupted the grouping.
#     The `group` column is now read by name, and is OPTIONAL - absent means one
#     cohort-wide group.
#
#   * Grouping semantics are explicit. Upstream evaluated the coverage criterion and
#     the editing-rate criterion in separate passes, combined each across groups with
#     full_join, and intersected the results - so a site could be kept because it had
#     coverage in group A and an adequate editing rate in group B. Here both criteria
#     must hold WITHIN THE SAME group, and --groupLogic decides whether one qualifying
#     group is enough ("any", the closest match to upstream) or all groups must qualify
#     ("all", what the upstream config comment actually described).

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inCov'), help = 'input coverage matrix (.tsv.gz)', default = ''),
    make_option(c('--inRat'), help = 'input editing-rate matrix (.tsv.gz)', default = ''),
    make_option(c('--outCov'), help = 'output coverage matrix (.tsv.gz)', default = ''),
    make_option(c('--outRat'), help = 'output editing-rate matrix (.tsv.gz)', default = ''),
    make_option(c('--metadata'), help = 'TSV with `sample` and optionally `group`', default = ''),
    make_option(c('--minSamples'), type = 'integer', default = 40L,
                help = 'site must have coverage in >= this many samples of a group [default %default]'),
    make_option(c('--minER'), type = 'double', default = 0.1,
                help = 'minimum mean editing rate within a group [default %default]'),
    make_option(c('--groupLogic'), default = 'any',
                help = '"any" = one qualifying group is enough, "all" = every group must qualify [default %default]'),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(opt$groupLogic %in% c('any', 'all'))
setDTthreads(opt$threads)

meta <- fread(opt$metadata, sep = "\t")
stopifnot("sample" %in% names(meta))
sample_ids <- as.character(meta$sample)

if ("group" %in% names(meta)) {
    groups <- as.character(meta$group)
} else {
    message(" * no `group` column in metadata - treating the cohort as a single group")
    groups <- rep("all", length(sample_ids))
}

message(" * reading coverage matrix")
cov_dt <- fread(opt$inCov, sep = "\t", showProgress = FALSE)
message(" * reading editing-rate matrix")
rat_dt <- fread(opt$inRat, sep = "\t", showProgress = FALSE)

write_empty <- function() {
    empty <- data.table(ESid = character())
    for (s in sample_ids) empty[, (s) := numeric()]
    fwrite(empty, opt$outCov, sep = "\t", compress = "gzip")
    fwrite(empty, opt$outRat, sep = "\t", compress = "gzip")
}

if (nrow(cov_dt) == 0L) {
    message(" * empty input - writing empty output")
    write_empty()
    quit(save = "no", status = 0)
}

stopifnot(identical(cov_dt$ESid, rat_dt$ESid))
esid <- cov_dt$ESid

# merge_jacusa.R writes columns in metadata order; assert rather than assume, because a
# silent column/sample mismatch would misassign every value in the matrix.
stopifnot(identical(setdiff(names(cov_dt), "ESid"), sample_ids))
stopifnot(identical(setdiff(names(rat_dt), "ESid"), sample_ids))

cov_m <- as.matrix(cov_dt[, ..sample_ids])
rat_m <- as.matrix(rat_dt[, ..sample_ids])
rm(cov_dt, rat_dt); invisible(gc())

group_levels <- unique(groups)
message(" * ", length(group_levels), " group(s): ", paste(group_levels, collapse = ", "))

pass <- NULL
for (g in group_levels) {
    cols <- which(groups == g)
    cov_ok <- rowSums(!is.na(cov_m[, cols, drop = FALSE])) >= opt$minSamples
    # All-NA rows give NaN from rowMeans; NaN >= x is NA, so coerce to FALSE.
    er <- rowMeans(rat_m[, cols, drop = FALSE], na.rm = TRUE)
    er_ok <- !is.na(er) & er >= opt$minER
    pass_g <- cov_ok & er_ok
    message("   - ", g, " (n=", length(cols), "): ", sum(pass_g), " sites qualify")
    pass <- if (is.null(pass)) pass_g else if (opt$groupLogic == 'any') pass | pass_g else pass & pass_g
}

message(" * keeping ", sum(pass), " of ", length(pass), " sites (groupLogic = ", opt$groupLogic, ")")

if (!any(pass)) {
    write_empty()
    quit(save = "no", status = 0)
}

write_matrix <- function(m, path) {
    out <- as.data.table(m[pass, , drop = FALSE])
    setnames(out, sample_ids)
    out[, ESid := esid[pass]]
    setcolorder(out, "ESid")
    fwrite(out, path, sep = "\t", compress = "gzip")
}

write_matrix(cov_m, opt$outCov)
write_matrix(rat_m, opt$outRat)
message(" * done")
