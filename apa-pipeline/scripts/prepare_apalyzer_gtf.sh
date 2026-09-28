#!/usr/bin/env bash
#
# Prepare a GENCODE-style GTF for APAlyzer's PAS2GEF().
#
# Usage: prepare_apalyzer_gtf.sh <annotation.gtf[.gz]> <out.gtf.gz>
#
# PAS2GEF was written for Ensembl GTFs: it prepends "chr" to every seqname and
# uses gene_name as a merge key between its PAS, intron, and last-exon tables.
# This script therefore:
#   - keeps chr1-22, chrX, chrY, chrM and strips the "chr" prefix
#   - drops *_PAR_Y transcripts
#   - adds gene_name to transcript/exon/CDS records that lack it (isoMiGA
#     PacBio transcripts), taken from the gene record, or gene_id when the
#     gene record has none (novelGene_*)

set -euo pipefail

gtf="$1"
out="$2"
[ -f "${gtf}" ] || { echo "Missing: ${gtf}" >&2; exit 1; }

reader=cat
case "${gtf}" in *.gz) reader="gzip -dc" ;; esac

# Pass 1 collects gene names; pass 2 rewrites.
${reader} "${gtf}" | awk -F'\t' -v OFS='\t' '
  function attr(key,    m) {
    if (match($9, key " \"[^\"]+\"")) return substr($9, RSTART + length(key) + 2, RLENGTH - length(key) - 3)
    return ""
  }
  NR == FNR {
    if ($0 !~ /^#/ && $3 == "gene") { g = attr("gene_id"); n = attr("gene_name"); gname[g] = (n != "") ? n : g }
    next
  }
  /^#/ { print; next }
  $1 !~ /^chr([0-9]+|X|Y|M)$/ { next }
  $9 ~ /_PAR_Y/ { next }
  {
    sub(/^chr/, "", $1)
    if ($9 !~ /gene_name "/) {
      g = attr("gene_id"); n = (g in gname) ? gname[g] : g
      $9 = $9 " gene_name \"" n "\";"
    }
    print
  }' <(${reader} "${gtf}") - | gzip > "${out}"

echo "$(gzip -dc "${out}" | grep -vc '^#') records -> ${out}"
