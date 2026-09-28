#!/usr/bin/env bash
#
# Build a QAPA 3' UTR library from a GTF with qapa build.
#
# Usage: build_qapa_utr_library.sh <annotation.gtf[.gz]> <identifiers.tsv> \
#          <gencode_polya.txt.gz> <polyasite_clusters.bed.gz> <out_utrs.bed>
#
# gencode_polya.txt.gz   UCSC table dump wgEncodeGencodePolyaV* for the genome
# polyasite_clusters     PolyASite 2.0 atlas clusters BED
#
# The genePred name2 column is set to the gene symbol so event IDs match the
# precompiled QAPA libraries. qapa quant parses event IDs on underscores, so
# underscores in symbols are replaced by hyphens. Gene IDs and names in the
# event table still come from the identifier table.
#
# QAPA_ENV  QAPA conda environment (qapa, gtfToGenePred, bedtools)
# THREADS   qapa build cores (default: 4)

set -euo pipefail

gtf="$1"
identifiers="$2"
gencode_polya="$3"
polyasite="$4"
out_bed="$5"

qapa_env="${QAPA_ENV:-/sc/arion/projects/pd-omics/alex/conda/qapa}"
threads="${THREADS:-4}"

for f in "${gtf}" "${identifiers}" "${gencode_polya}" "${polyasite}"; do
  [ -f "${f}" ] || { echo "Missing: ${f}" >&2; exit 1; }
done
export PATH="${qapa_env}/bin:${PATH}"

work="$(mktemp -d "${TMPDIR:-/tmp}/qapa_build.XXXXXX")"
trap 'rm -rf "${work}"' EXIT

echo "[1/3] genePred"
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons "${gtf}" "${work}/raw.genePred"
# name2 -> gene symbol; drop PAR_Y duplicates.
awk -F'\t' -v OFS='\t' 'NR == FNR { if (FNR > 1) sym[$1] = $5; next }
  $1 ~ /_PAR_Y$/ { next }
  { g = $12; if (g ~ /^ENS/) sub(/\.[0-9]+$/, "", g)
    s = (g in sym) ? sym[g] : $12; gsub(/_/, "-", s); $12 = s; print }' \
  "${identifiers}" "${work}/raw.genePred" > "${work}/annotation.genePred"
echo "  $(wc -l < "${work}/annotation.genePred") transcripts"

echo "[2/3] poly(A) site annotations"
# UCSC genePred dump -> BED6 with the site class as name.
gzip -dc "${gencode_polya}" | awk -F'\t' -v OFS='\t' '{ print $3, $5, $6, $13, $12, $4 }' \
  > "${work}/gencode_polya.bed"
gzip -dc "${polyasite}" > "${work}/polyasite.bed"
echo "  $(grep -c polyA_site "${work}/gencode_polya.bed") GENCODE polyA_site, $(wc -l < "${work}/polyasite.bed") PolyASite clusters"

echo "[3/3] qapa build"
qapa build --db "${identifiers}" -g "${work}/gencode_polya.bed" -p "${work}/polyasite.bed" \
  -t "${work}" -c "${threads}" "${work}/annotation.genePred" > "${out_bed}"
echo "  $(wc -l < "${out_bed}") UTR entries -> ${out_bed}"
