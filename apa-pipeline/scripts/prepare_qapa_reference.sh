#!/usr/bin/env bash
#
# Build the Salmon index and stage the identifier table for one QAPA library.
#
# Usage: prepare_qapa_reference.sh <utrs.bed[.gz]> <identifiers.tsv> <genome.fa> <out_dir>
#
# Writes to out_dir:
#   qapa_3utrs.bed            UTR library
#   ensembl_identifiers.tsv   identifier table (config: identifiers)
#   qapa_3utrs.fa             3' UTR sequences (qapa fasta)
#   salmon_index/             Salmon index (config: salmon_index)
#
# QAPA_ENV  QAPA conda environment (qapa, bedtools)
# SALMON    salmon binary
# THREADS   salmon index threads (default: 4)

set -euo pipefail

utr_bed="$1"
identifiers_in="$2"
genome="$3"
out_dir="$4"

qapa_env="${QAPA_ENV:-/sc/arion/projects/pd-omics/alex/conda/qapa}"
salmon="${SALMON:-/sc/arion/projects/ad-omics/data/software/salmon-1.9.0_linux_x86_64/bin/salmon}"
threads="${THREADS:-4}"

for f in "${utr_bed}" "${identifiers_in}" "${genome}"; do
  [ -f "${f}" ] || { echo "Missing: ${f}" >&2; exit 1; }
done
[ -x "${salmon}" ] || { echo "Salmon not found: ${salmon}" >&2; exit 1; }

mkdir -p "${out_dir}"
export PATH="${qapa_env}/bin:${PATH}"

bed="${out_dir}/qapa_3utrs.bed"
identifiers="${out_dir}/ensembl_identifiers.tsv"
fasta="${out_dir}/qapa_3utrs.fa"
index="${out_dir}/salmon_index"

echo "[1/3] UTR library and identifier table"
case "${utr_bed}" in
  *.gz) gzip -dc "${utr_bed}" > "${bed}" ;;
  *) [ "$(readlink -f "${utr_bed}")" = "$(readlink -f "${bed}")" ] || cp "${utr_bed}" "${bed}" ;;
esac
[ "$(readlink -f "${identifiers_in}")" = "$(readlink -f "${identifiers}")" ] || cp "${identifiers_in}" "${identifiers}"
echo "  $(wc -l < "${bed}") UTR entries, $(($(wc -l < "${identifiers}") - 1)) transcripts"

# Every transcript in the BED must be in the identifier table.
missing="$(cut -f4 "${bed}" | cut -d_ -f1 | tr ',' '\n' | sort -u |
  comm -23 - <(tail -n +2 "${identifiers}" | cut -f2 | sort -u) | wc -l)"
[ "${missing}" -eq 0 ] || {
  echo "  ${missing} BED transcripts are absent from the identifier table" >&2
  exit 1
}
echo "  all BED transcripts present in identifier table"

echo "[2/3] qapa fasta"
qapa fasta -f "${genome}" "${bed}" "${fasta}"
echo "  $(grep -c '^>' "${fasta}") sequences"

echo "[3/3] salmon index"
rm -rf "${index}"
"${salmon}" index -t "${fasta}" -i "${index}" -p "${threads}" > "${out_dir}/salmon_index.log" 2>&1
echo "  $(basename "${index}") written"

echo
echo "Done. Config paths:"
echo "  salmon_index: ${index}"
echo "  identifiers:  ${identifiers}"
