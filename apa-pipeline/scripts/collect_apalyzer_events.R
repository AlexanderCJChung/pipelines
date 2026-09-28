#!/usr/bin/env Rscript
#
# Stack the per-sample APAlyzer tables without calculating new statistics.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop(
    "usage: collect_apalyzer_events.R <3utr_output> <ipa_output> ",
    "<3utr_files...> --ipa <ipa_files...>"
  )
}

separator <- match("--ipa", args)
if (is.na(separator) || separator < 4L || separator == length(args)) {
  stop("Expected 3' UTR files followed by --ipa and IPA files")
}

utr_output <- args[1]
ipa_output <- args[2]
utr_files <- args[3:(separator - 1L)]
ipa_files <- args[(separator + 1L):length(args)]

sample_from_path <- function(path) {
  basename(dirname(dirname(dirname(path))))
}

read_event_files <- function(files) {
  dplyr::bind_rows(lapply(files, function(path) {
    sample_id <- sample_from_path(path)
    event_table <- readr::read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    )

    # APAlyzer prefixes metric columns with the sample ID. Remove that prefix
    # because the sample is stored in its own column.
    prefix <- paste0(sample_id, "_")
    prefixed <- startsWith(names(event_table), prefix)
    names(event_table)[prefixed] <- substring(
      names(event_table)[prefixed],
      nchar(prefix) + 1L
    )
    event_table$sample <- sample_id
    event_table
  }))
}

utr_samples <- vapply(utr_files, sample_from_path, character(1))
ipa_samples <- vapply(ipa_files, sample_from_path, character(1))
if (anyDuplicated(utr_samples) || anyDuplicated(ipa_samples)) {
  stop("Each sample must have one 3' UTR file and one IPA file")
}
if (!setequal(utr_samples, ipa_samples)) {
  stop("3' UTR and IPA files do not contain the same samples")
}

dir.create(dirname(utr_output), recursive = TRUE, showWarnings = FALSE)
utr_events <- read_event_files(utr_files)
ipa_events <- read_event_files(ipa_files)

readr::write_tsv(utr_events, utr_output)
readr::write_tsv(ipa_events, ipa_output)

message("[collect_apalyzer_events] 3' UTR rows: ", nrow(utr_events))
message("[collect_apalyzer_events] IPA rows: ", nrow(ipa_events))
