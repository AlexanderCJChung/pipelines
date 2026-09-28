#!/usr/bin/env Rscript
# Merge all per-sample Alu Editing Index outputs into one table.
#
# Ported from RajLabMSSM/editing-pipeline (Jack Humphrey, 2021).
#
# Two fixes: upstream called write_tsv(path = ...), and `path` was removed in readr 2.0
# so the script errored on any modern R; and no Snakemake rule ever invoked it, so the
# merged table was never produced (see the merge_aei rule in the Snakefile).
#
# RNAEditingIndex writes one EditingIndex.csv per sample directory. The sample name is
# not reliably inside the file, so it is taken from the directory the file was found in.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inFolder', '-i'), help = 'directory containing per-sample subfolders', default = ''),
    make_option(c('--output', '-o'), help = 'output TSV', default = '')
)

opt <- parse_args(OptionParser(option_list = option_list))

all_files <- list.files(opt$inFolder, pattern = "^EditingIndex\\.csv$",
                        recursive = TRUE, full.names = TRUE)
message(" * found ", length(all_files), " AEI files")
if (length(all_files) == 0L) stop("no EditingIndex.csv found under ", opt$inFolder)

# .../{outDir}/{sample}/{...}/EditingIndex.csv - walk up to the child of inFolder.
sample_of <- function(path) {
    rel <- sub(paste0("^", normalizePath(opt$inFolder), "/?"), "",
               normalizePath(path))
    strsplit(rel, .Platform$file.sep, fixed = TRUE)[[1]][1]
}

all_res <- rbindlist(lapply(all_files, function(f) {
    d <- fread(f, showProgress = FALSE)
    d[, sample := sample_of(f)]
    d
}), use.names = TRUE, fill = TRUE)

setcolorder(all_res, "sample")
message(" * merged ", nrow(all_res), " rows across ", uniqueN(all_res$sample), " samples")

fwrite(all_res, opt$output, sep = "\t")
