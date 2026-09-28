#!/usr/bin/env bash
#
# Derive the five-column identifier table QAPA expects (--db) from a GTF.
#
# Usage: qapa_identifiers_from_gtf.sh <annotation.gtf[.gz]> > ensembl_identifiers.tsv
#
# Equivalent to the BioMart export in the QAPA documentation for GENCODE
# transcripts. Versions are stripped from ENSG/ENST IDs because qapa build
# strips them before joining. PAR_Y transcripts are dropped.
#
# Transcripts without gene_type/transcript_type attributes (isoMiGA MSTRG
# transcripts) inherit gene type and gene name from their gene record. Their
# transcript type is protein_coding when a CDS is annotated, otherwise
# processed_transcript. Genes without a gene record type (novelGene_*) are
# protein_coding when any transcript has a CDS. qapa build keeps only
# transcripts where both types are protein_coding.

set -euo pipefail

gtf="$1"
[ -f "${gtf}" ] || { echo "Missing: ${gtf}" >&2; exit 1; }

reader=cat
case "${gtf}" in *.gz) reader="gzip -dc" ;; esac

${reader} "${gtf}" | awk -F'\t' '
  function attr(key,    m) {
    if (match($9, key " \"[^\"]+\"")) {
      m = substr($9, RSTART + length(key) + 2, RLENGTH - length(key) - 3)
      return m
    }
    return ""
  }
  function strip(id) { if (id ~ /^ENS/) sub(/\.[0-9]+$/, "", id); return id }
  BEGIN { OFS = "\t" }
  /^#/ { next }
  $3 == "gene" {
    g = attr("gene_id"); gtype[g] = attr("gene_type"); gname[g] = attr("gene_name")
    next
  }
  $3 == "transcript" {
    t = attr("transcript_id"); if (t ~ /_PAR_Y$/) next
    n++; tid[n] = t; tgene[n] = attr("gene_id"); ttype[n] = attr("transcript_type")
    next
  }
  $3 == "CDS" { has_cds[attr("transcript_id")] = 1; g = attr("gene_id"); gene_cds[g] = 1 }
  END {
    print "Gene stable ID", "Transcript stable ID", "Gene type", "Transcript type", "Gene name"
    for (i = 1; i <= n; i++) {
      t = tid[i]; g = tgene[i]
      gt = (g in gtype && gtype[g] != "") ? gtype[g] : ((g in gene_cds) ? "protein_coding" : "unknown")
      tt = ttype[i] != "" ? ttype[i] : ((t in has_cds) ? "protein_coding" : "processed_transcript")
      gn = (g in gname && gname[g] != "") ? gname[g] : g
      print strip(g), strip(t), gt, tt, gn
    }
  }'
