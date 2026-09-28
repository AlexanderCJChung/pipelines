#!/usr/bin/env Rscript
# Release checks and QC summary for the editing pipeline. Fails the release if:
#   * matrix columns are not the sample sheet (minus samples dropped by the missingness
#     filter) in sheet order;
#   * A>G (gene-oriented) is under --minAG of the final sites - anything else means the
#     strand flip or the SNP filters are wrong;
#   * there are no true zeros (the back-fill must separate unedited from unmeasured);
#   * an AEI sample is missing.
# Run-specific numbers go to --summary / --spectrum instead of the methods doc.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--samples"), make_option("--editing", default = ""),
    make_option("--annotation", default = ""), make_option("--aei", default = ""),
    make_option("--minAG", type = "double", default = 0.8),
    make_option("--summary"), make_option("--spectrum"), make_option("--passed")
)))

sheet <- fread(opt$samples, sep = "\t")$sample
summary <- list(samples_in_sheet = length(sheet))
spectrum <- data.table(substitution = character(), sites = integer(), fraction = numeric())
fail <- character()

if (nzchar(opt$editing)) {
    ed <- fread(opt$editing)
    anno <- fread(opt$annotation, select = c("ESid", "Ref_gene", "Alt_gene", "REDIportal_info"))
    cols <- setdiff(names(ed), "ESid")
    if (!identical(cols, sheet[sheet %in% cols]))
        fail <- c(fail, "matrix columns are not the sample sheet in sheet order")
    if (!identical(ed$ESid, anno$ESid)) fail <- c(fail, "editing and annotation site order differ")

    m <- as.matrix(ed[, ..cols])
    spectrum <- anno[, .(sites = .N), by = .(substitution = paste0(Ref_gene, ">", Alt_gene))]
    spectrum[, fraction := sites / sum(sites)]
    setorder(spectrum, -sites)
    ag <- spectrum[substitution == "A>G", sum(fraction)]
    if (ag < opt$minAG) fail <- c(fail, sprintf("A>G fraction %.3f < %.2f", ag, opt$minAG))
    true_zero <- sum(m == 0, na.rm = TRUE)
    if (true_zero == 0) fail <- c(fail, "no true zeros in the editing matrix")

    dropped <- setdiff(sheet, cols)
    summary <- c(summary, list(
        sites = nrow(ed), samples = length(cols),
        samples_dropped_missingness = length(dropped),
        dropped = paste(dropped, collapse = ","),
        a_to_g_fraction = round(ag, 4),
        rediportal_fraction = round(mean(!anno$REDIportal_info %in% c("", ".", NA)), 4),
        cells = length(m), na_cells = sum(is.na(m)), true_zero_cells = true_zero,
        median_editing_rate = round(median(m, na.rm = TRUE), 4)
    ))
}

if (nzchar(opt$aei)) {
    aei <- fread(opt$aei)
    missing <- setdiff(sheet, aei$sample)
    if (length(missing)) fail <- c(fail, paste("AEI missing for", length(missing), "samples"))
    summary <- c(summary, list(
        aei_samples = uniqueN(aei$sample),
        aei_a2g_median = round(median(aei$A2GEditingIndex), 4),
        aei_c2t_median = round(median(aei$C2TEditingIndex), 4)
    ))
}

fwrite(data.table(metric = names(summary), value = vapply(summary, as.character, "")),
       opt$summary, sep = "\t")
fwrite(spectrum, opt$spectrum, sep = "\t")
if (length(fail)) stop("validation failed:\n  ", paste(fail, collapse = "\n  "))
writeLines(c("passed", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), opt$passed)
message(" * validation passed")
