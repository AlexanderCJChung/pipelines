#!/usr/bin/env Rscript
#
# Build the three objects APAlyzer counts against: UTRdbraw, dfIPA, and dfLE.
#
# The source is either an official APAlyzer reference RData file (PolyA_DB
# based; objects refUTRraw_<name>, dfIPA_<name>, dfLE_<name>) or a GTF, from
# which PAS2GEF() derives poly(A) sites as annotated transcript ends. GTF input
# must follow Ensembl conventions (no "chr" prefix; see
# prepare_apalyzer_gtf.sh). REF4PAS() is applied in both cases.

suppressPackageStartupMessages({
  library(APAlyzer)
  # PAS2GEF's internals call queryHits() without importing it from S4Vectors
  # (APAlyzer 1.20.0); attaching GenomicRanges puts it on the search path.
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: prepare_apalyzer_reference.R <source.RData> <reference_name> <output.RData>")
}

source_file <- args[1]
reference_name <- args[2]
output_file <- args[3]

if (!file.exists(source_file)) {
  stop("APAlyzer reference not found: ", source_file)
}

reference_env <- new.env(parent = emptyenv())
if (grepl("\\.gtf(\\.gz)?$", source_file, ignore.case = TRUE)) {
  message("[prepare_apalyzer_reference] PAS2GEF on ", source_file)
  raw <- APAlyzer::PAS2GEF(source_file, AnnoMethod = "V2")
  assign("refUTRraw", raw$refUTRraw, envir = reference_env)
  assign("dfIPAraw", raw$dfIPA, envir = reference_env)
  assign("dfLEraw", raw$dfLE, envir = reference_env)
} else {
  load(source_file, envir = reference_env)
}

find_object <- function(candidates, label) {
  found <- candidates[candidates %in% ls(reference_env)]
  if (length(found) == 0L) {
    stop(
      "Could not find ", label, " in ", source_file, ". Objects present: ",
      paste(ls(reference_env), collapse = ", ")
    )
  }
  get(found[1], envir = reference_env)
}

suffix <- paste0("_", reference_name)
ref_utr <- find_object(
  c(paste0("refUTRraw", suffix), "refUTRraw", "UTRdbraw"),
  "3' UTR reference"
)
ref_ipa <- find_object(
  c(paste0("dfIPA", suffix), "dfIPAraw", "dfIPA"),
  "intronic poly(A) reference"
)
ref_le <- find_object(
  c(paste0("dfLE", suffix), "dfLEraw", "dfLE"),
  "last-exon reference"
)

# Official reference files contain raw regions. REF4PAS applies APAlyzer's
# required region preparation. Already prepared files use the standard names.
if (all(c("UTRdbraw", "dfIPA", "dfLE") %in% ls(reference_env))) {
  UTRdbraw <- ref_utr
  dfIPA <- ref_ipa
  dfLE <- ref_le
} else {
  prepared <- APAlyzer::REF4PAS(ref_utr, ref_ipa, ref_le)
  UTRdbraw <- prepared$UTRdbraw
  dfIPA <- prepared$dfIPA
  dfLE <- prepared$dfLE
}

# PASEXP_IPA computes abs(TES - LEstart). Some reference files store these as
# text, which makes that subtraction fail only after the 3' UTR step has already
# written its output. Convert once here instead of per sample.
for (column in c("LEstart", "TES")) {
  if (column %in% names(dfLE) && !is.numeric(dfLE[[column]])) {
    converted <- suppressWarnings(as.numeric(as.character(dfLE[[column]])))
    lost <- sum(is.na(converted) & !is.na(dfLE[[column]]))
    if (lost > 0L) {
      stop("dfLE$", column, " has ", lost, " non-numeric value(s)")
    }
    dfLE[[column]] <- converted
    message("[prepare_apalyzer_reference] converted dfLE$", column, " to numeric")
  }
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
save(UTRdbraw, dfIPA, dfLE, file = output_file)

message("[prepare_apalyzer_reference] reference: ", reference_name)
utr_regions <- unlist(UTRdbraw)
message("[prepare_apalyzer_reference] 3' UTR regions (aUTR + cUTR): ", length(utr_regions),
        " in ", length(unique(as.character(utr_regions$gene_symbol))), " genes")
message("[prepare_apalyzer_reference] IPA rows: ", nrow(dfIPA))
message("[prepare_apalyzer_reference] last-exon rows: ", nrow(dfLE))
