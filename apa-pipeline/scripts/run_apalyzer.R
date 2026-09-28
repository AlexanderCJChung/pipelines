#!/usr/bin/env Rscript
#
# Quantify 3' UTR and intronic APA events for one sample.
#
#   Rscript run_apalyzer.R <refs.RData> <bam> <sample> <outdir> <strandtype> <seqtype> <threads>
#
# Writes <outdir>/3utr.tsv and <outdir>/ipa.tsv.
#
suppressPackageStartupMessages({
  library(APAlyzer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: run_apalyzer.R <refs.RData> <bam> <sample> <outdir> [strandtype] [seqtype] [threads]")
}
ref_file   <- args[1]
bam_file   <- args[2]
sample_id  <- args[3]
out_dir    <- args[4]
strandtype <- if (length(args) >= 5) args[5] else "NONE"
seqtype    <- if (length(args) >= 6) args[6] else "SingleEnd"
threads    <- if (length(args) >= 7) as.integer(args[7]) else 1L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(bam_file)) stop("BAM not found: ", bam_file)
bai <- c(paste0(bam_file, ".bai"), sub("\\.bam$", ".bai", bam_file))
if (!any(file.exists(bai))) stop("no BAM index (.bai) beside: ", bam_file)

# Rsubread writes a temporary annotation file into the working directory, which
# is the pipeline checkout under Snakemake. Resolve the paths, then work from a
# per-process temporary directory.
ref_file <- normalizePath(ref_file, mustWork = TRUE)
bam_file <- normalizePath(bam_file, mustWork = TRUE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
setwd(tempdir())

if (!strandtype %in% c("forward", "invert", "NONE")) {
  stop("Strandtype must be one of forward / invert / NONE, got: ", strandtype)
}
if (!seqtype %in% c("SingleEnd", "ThreeMostPairEnd")) {
  stop("SeqType must be SingleEnd or ThreeMostPairEnd, got: ", seqtype)
}

# Check the BAM layout before counting.
.bf <- Rsamtools::BamFile(bam_file, yieldSize = 10000L)
open(.bf)
.flags <- Rsamtools::scanBam(.bf, param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag
close(.bf)
if (length(.flags) == 0L) stop("BAM appears to contain no reads: ", bam_file)
paired_frac <- mean(bitwAnd(.flags, 1L) > 0L)
message("[apalyzer] library check: ", round(100 * paired_frac, 1),
        "% of the first ", length(.flags), " reads are paired")

if (paired_frac > 0.5 && seqtype != "ThreeMostPairEnd") {
  stop("BAM is PAIRED-END (", round(100 * paired_frac, 1), "% of sampled reads) but ",
       "apalyzer_seqtype is '", seqtype,
       "'. Set apalyzer_seqtype: \"ThreeMostPairEnd\" in config.yaml.")
}
if (paired_frac <= 0.5 && seqtype != "SingleEnd") {
  stop("BAM is SINGLE-END (", round(100 * paired_frac, 1), "% of sampled reads paired) ",
       "but apalyzer_seqtype is '", seqtype, "'. Set apalyzer_seqtype: \"SingleEnd\".")
}

message("[apalyzer] sample=", sample_id, " strand=", strandtype,
        " seq_type=", seqtype, " threads=", threads)

load(ref_file)   # UTRdbraw, dfIPA, dfLE
for (obj in c("UTRdbraw", "dfIPA", "dfLE")) {
  if (!exists(obj)) stop("reference file is missing object '", obj, "': ", ref_file)
}

flS <- c(bam_file)
names(flS) <- sample_id

# Stop if a reference mismatch produces zero counts.
assert_counted <- function(df, pattern, what) {
  cols <- grep(pattern, names(df), value = TRUE)
  if (length(cols) == 0L) {
    stop("no ", what, " read columns found. Columns: ", paste(names(df), collapse = ", "))
  }
  total <- sum(as.numeric(unlist(df[, cols, drop = FALSE])), na.rm = TRUE)
  message("[apalyzer] ", what, " reads counted: ", format(total, big.mark = ","))
  if (total == 0) {
    stop("all ", what, " counts are zero for ", sample_id,
         ". Check the reference, BAM contig names, and strand setting.")
  }
}

message("[apalyzer] PASEXP_3UTR")
utr <- APAlyzer::PASEXP_3UTR(UTRdbraw, flS, Strandtype = strandtype)
utr$sample <- sample_id
if (nrow(utr) == 0L) stop("PASEXP_3UTR returned no rows for ", sample_id)
assert_counted(utr, "(areads|creads)$", "3'UTR")
write.table(utr, file.path(out_dir, "3utr.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[apalyzer] 3'UTR rows: ", nrow(utr))

message("[apalyzer] PASEXP_IPA")
ipa <- APAlyzer::PASEXP_IPA(dfIPA, dfLE, flS,
                            Strandtype = strandtype, nts = threads, SeqType = seqtype)
ipa$sample <- sample_id
assert_counted(ipa, "(UPreads|DNreads|LEreads)$", "IPA")
write.table(ipa, file.path(out_dir, "ipa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[apalyzer] IPA rows: ", nrow(ipa))

message("[apalyzer] OK")
