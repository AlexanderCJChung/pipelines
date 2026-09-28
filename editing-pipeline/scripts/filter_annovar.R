#!/usr/bin/env Rscript
# Filter the ANNOVAR annotation down to the final editing-site set and emit the BED
# that drives the pileup back-fill.
#
# Ported from RajLabMSSM/editing-pipeline (Cuddleston, Humphrey, Seo).
#
# Two changes from upstream:
#
#   * The gnomAD cutoff is a parameter. Upstream required `AF == "."`, i.e. the site
#     must be absent from gnomAD entirely, which also discards genuine editing sites
#     that happen to carry a singleton call in a population database. --maxGnomadAF
#     keeps the upstream behaviour at its default of 0.
#
#   * Multi-gene annotations are resolved instead of dropped. ANNOVAR writes
#     Gene.refGene as "GENE1;GENE2" for sites spanning two genes and as
#     "NONE(dist=NONE),LINC01128(dist=1000)" for intergenic sites with --neargene.
#     Upstream joined the raw string to GENCODE gene_name, so every one of those failed
#     to match and was dropped by the !is.na(strand) filter. Here the string is split
#     and the first token matching a GENCODE gene is used.
#
# Sites with no resolvable gene are still dropped: strand is what lets merge_pileup
# orient A>G against T>C, and an unstranded library gives no other source for it.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

option_list <- list(
    make_option(c('--inAnno'), help = 'ANNOVAR *_multianno.txt', default = ''),
    make_option(c('--gencode'), help = 'GENCODE gene_strand TSV (ensembl_id, gene_name, strand)', default = ''),
    make_option(c('--outAnno'), help = 'output annotation TSV', default = ''),
    make_option(c('--outBed'),  help = 'output BED of surviving sites', default = ''),
    make_option(c('--maxGnomadAF'), type = 'double', default = 0,
                help = paste('drop sites with gnomAD AF above this. 0 reproduces upstream,',
                             'which required absence from gnomAD entirely [default %default]')),
    make_option(c('--threads'), type = 'integer', default = 1L,
                help = 'data.table threads [default %default]')
)

opt <- parse_args(OptionParser(option_list = option_list))
setDTthreads(opt$threads)

anno <- fread(opt$inAnno, sep = "\t", showProgress = FALSE, colClasses = "character")
message(" * read in ", nrow(anno), " sites")

# --- common SNP filter -------------------------------------------------------------
stopifnot(all(c("dbsnp153CommonSNV", "Chr", "Start", "End", "Ref", "Alt") %in% names(anno)))
anno <- anno[dbsnp153CommonSNV == "."]
message(" * ", nrow(anno), " sites not overlapping common dbSNP variants")

# `AF` is the first gnomAD column table_annovar emits for gnomad30_genome. It is absent
# when gnomad30_genome is dropped from annovar_protocols to avoid the 43 GB scan - in
# that case the common-dbSNP filter above is the only SNP filter, which is looser.
if ("AF" %in% names(anno)) {
    af <- suppressWarnings(as.numeric(anno$AF))   # "." -> NA, meaning absent from gnomAD
    anno <- anno[is.na(af) | af <= opt$maxGnomadAF]
    message(" * ", nrow(anno), " sites with gnomAD AF <= ", opt$maxGnomadAF, " or absent")
} else {
    message(" * NOTE: no `AF` column - gnomad30_genome was not in annovar_protocols, ",
            "so no population-frequency filter was applied")
}

if (nrow(anno) == 0L) stop("no sites survived the SNP filters - check the ANNOVAR protocols")

anno[, ESid := paste0(Chr, ":", Start, ":", Ref, ":", Alt)]

# --- gene strand from GENCODE ------------------------------------------------------
message(" * adding GENCODE gene strand")
gencode <- fread(opt$gencode, sep = "\t", showProgress = FALSE)
stopifnot(all(c("ensembl_id", "gene_name", "strand") %in% names(gencode)))
gencode <- gencode[!duplicated(gene_name)]
gene_strand <- setNames(gencode$strand, gencode$gene_name)
gene_ens    <- setNames(gencode$ensembl_id, gencode$gene_name)

# "GENE1;GENE2" and "NONE(dist=NONE),LINC01128(dist=1000)" both need splitting, and the
# parenthesised distance suffix stripping, before anything will match a GENCODE name.
resolve_gene <- function(x) {
    toks <- strsplit(x, "[;,]")
    vapply(toks, function(tk) {
        tk <- sub("\\(.*$", "", tk)
        hit <- tk[tk %in% names(gene_strand)]
        if (length(hit)) hit[1] else NA_character_
    }, character(1))
}

anno[, gene_resolved := resolve_gene(Gene.refGene)]
anno[, strand := unname(gene_strand[gene_resolved])]
anno[, ensembl_id := unname(gene_ens[gene_resolved])]

n_before <- nrow(anno)
anno <- anno[!is.na(strand)]
message(" * ", nrow(anno), " sites mapped to a GENCODE gene (dropped ",
        n_before - nrow(anno), " with no resolvable gene)")

if (nrow(anno) == 0L) stop("no sites mapped to genes - check the GENCODE table")

# --- outputs -----------------------------------------------------------------------
# BED is 0-based half-open, so start = Start - 1. Column 4 carries the ESid so the
# pileup output can be traced back to the exact site and allele.
bed <- anno[, .(Chr, start = as.integer(Start) - 1L, end = as.integer(End), ESid)]
setorder(bed, Chr, start)
fwrite(bed, opt$outBed, sep = "\t", col.names = FALSE)

fwrite(anno, opt$outAnno, sep = "\t")
message(" * done")
