#!/usr/bin/env Rscript
# Union the de novo calls with the REDIportal known-site calls for one sample.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo). Unchanged in
# behaviour; data.table instead of dplyr::union for speed, and it tolerates either
# input being empty (upstream's read_tsv errored on a header-only file).
#
# The union is what recovers known editing sites that fall below the stricter de novo
# depth threshold - the known-site pass uses --siteDepth 5 against the de novo 10.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--deNovoIN'), help = 'parsed de novo calls (.filt)', default = ''),
    make_option(c('--knownIN'),  help = 'parsed known-site calls (.knownES.filt)', default = ''),
    make_option(c('--output'),   help = 'output TSV (.filtAll)', default = '')
)

opt <- parse_args(OptionParser(option_list = option_list))

for (f in c(opt$deNovoIN, opt$knownIN)) {
    if (!file.exists(f)) stop("File ", f, " does not exist")
}

de_novo <- fread(opt$deNovoIN, sep = "\t", showProgress = FALSE)
known   <- fread(opt$knownIN,  sep = "\t", showProgress = FALSE)

all_sites <- unique(rbindlist(list(de_novo, known), use.names = TRUE, fill = TRUE))
message(" * ", nrow(de_novo), " de novo + ", nrow(known), " known -> ",
        nrow(all_sites), " unique sites")

fwrite(all_sites, opt$output, sep = "\t")
