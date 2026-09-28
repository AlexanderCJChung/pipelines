#!/usr/bin/env Rscript
# Parse a JACUSA2 call-1 result into one row per candidate editing site.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo).
# Two upstream scripts, parse_jacusa.R and parse_knownSites.R, were byte-identical
# apart from the default --siteDepth (10 vs 5), so they are merged here into one
# script and the threshold comes from config. See README "Divergences from upstream".
#
# Reimplemented in data.table. The upstream tidyverse version does a pivot_longer
# over a genome-wide JACUSA output (tens of millions of rows) once per sample; at
# 583 samples that was the dominant cost of the discovery pass.
#
# JACUSA2 BED6-extended output: line 1 is a track line, line 2 is the header
#   #contig start end name score strand bases11 info filter ref
# and bases11 is a comma-separated A,C,G,T count vector for the site.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--input'),  help = 'JACUSA2 call-1 result file', default = ''),
    make_option(c('--output'), help = 'output TSV', default = ''),
    make_option(c('--siteDepth'), type = 'integer', default = 10L,
                help = 'minimum total read coverage of the locus [default %default]'),
    make_option(c('--altDepth'), type = 'integer', default = 3L,
                help = 'minimum reads supporting the alternative allele [default %default]'),
    make_option(c('--minRate'), type = 'double', default = 0.01,
                help = 'minimum editing rate [default %default]'),
    make_option(c('--multiallelic'), default = 'drop',
                help = paste('how to resolve loci with >1 non-reference allele:',
                             '"drop" discards the locus (upstream 91d3898 behaviour),',
                             '"top" keeps the allele with highest coverage',
                             '(pre-91d3898 behaviour) [default %default]')),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(opt$multiallelic %in% c('drop', 'top'))
setDTthreads(opt$threads)

if (!file.exists(opt$input)) stop("File ", opt$input, " does not exist")
message("Loading JACUSA2 output from ", opt$input)

dt <- fread(opt$input, skip = 1L, header = TRUE, showProgress = FALSE)

# An empty result is legitimate - a sample can contribute no sites on a small
# contig - but fread returns a 0-column table for a header-only file, so the
# column references below would error. Emit a correctly-typed empty file instead.
empty_out <- function() {
    fwrite(data.table(ESid = character(), chr = character(), pos = integer(),
                      score = numeric(), total_cov = integer(), ref_cov = integer(),
                      alt_cov = integer(), edit_rate = numeric()),
           opt$output, sep = "\t")
    message(" * 0 sites in input - wrote empty output")
    quit(save = "no", status = 0)
}
if (nrow(dt) == 0L || !"bases11" %in% names(dt)) empty_out()

message(" * splitting counts by base")
dt[, c("A", "C", "G", "T") := tstrsplit(bases11, ",", fixed = TRUE, type.convert = TRUE)]

# `end` is the 1-based position (the BED interval is [end-1, end)).
setnames(dt, "#contig", "chr")
dt <- dt[, .(chr, pos = end, score, ref, A, C, G, T)]
dt[, total_cov := A + C + G + T]
dt[, ref_cov := fcase(ref == "A", A, ref == "C", C, ref == "G", G, ref == "T", T,
                      default = NA_integer_)]
# JACUSA emits the reference base from the FASTA; anything not in ACGT (N, masked)
# cannot yield an interpretable editing call.
dt <- dt[!is.na(ref_cov)]
dt[, chrpos := paste0(chr, ":", pos)]

message(" * ", nrow(dt), " sites loaded")

long <- melt(dt, id.vars = c("chrpos", "chr", "pos", "score", "ref", "ref_cov", "total_cov"),
             measure.vars = c("A", "C", "G", "T"),
             variable.name = "alt", value.name = "alt_cov",
             variable.factor = FALSE)
long <- long[ref != alt & alt_cov > 0L]

# Order by locus, then by descending alt coverage, so "top" can take the first row.
setorder(long, chrpos, -alt_cov)
if (opt$multiallelic == "drop") {
    long[, n_alt := .N, by = chrpos]
    long <- long[n_alt == 1L][, n_alt := NULL]
} else {
    long <- long[!duplicated(chrpos)]
}

out <- long[total_cov >= opt$siteDepth & alt_cov >= opt$altDepth]
out[, edit_rate := alt_cov / total_cov]
out <- out[edit_rate >= opt$minRate]
out[, ESid := paste0(chrpos, ":", ref, ":", alt)]
out <- out[, .(ESid, chr, pos, score, total_cov, ref_cov, alt_cov, edit_rate)]

message(" * ", nrow(out), " sites pass thresholds")
fwrite(out, opt$output, sep = "\t")
