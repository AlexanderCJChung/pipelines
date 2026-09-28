# QAPA reference build log

Two QAPA 3′ UTR libraries were built on 2026-09-10 under
`/sc/arion/projects/pd-omics/Mega_MyND/references/apa/<reference_name>/qapa/`.
Every command below was run as written, from the pipeline checkout, on a
Minerva interactive node. The pinned inputs in this directory are the exact
files used.

## Software

| Tool | Version | Location |
|---|---|---|
| QAPA | 1.4.1 | `/sc/arion/projects/pd-omics/alex/conda/qapa` |
| bedtools (qapa fasta) | 2.29.2 | same env |
| gtfToGenePred (UCSC) | conda `ucsc-gtftogenepred` | same env |
| pandas / pybedtools | 2.3.3 / 0.10.0 | same env |
| Salmon | 1.9.0 | `/sc/arion/projects/ad-omics/data/software/salmon-1.9.0_linux_x86_64/bin/salmon` |

Salmon 1.9.0 is the version `qapa_salmon_quant` uses, so index and
quantification match.

## Inputs

| File | Source | md5 |
|---|---|---|
| `qapa_3utrs.gencode_V31.hg38.bed.gz` | https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_V31.hg38.bed.gz (QAPA's precompiled hg38 library, GENCODE V31) | `b9b5ec2a495836473bfb93d488e0cf63` |
| `wgEncodeGencodePolyaV38.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodePolyaV38.txt.gz (GENCODE V38 poly(A) track, UCSC table dump) | `f613fe83197fae97e9d4a56bd055d57c` |
| `atlas.clusters.2.0.GRCh38.96.bed.gz` | https://www.polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz (PolyASite 2.0) | `d78001b7c10b8fe7e2c847aa32e86425` |
| `gencode.v31.annotation.gtf.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz — kept next to the hg38 build, not in this directory (41 MB) | `0aac86e8a6c9e5fa5afe2a286325da81` |
| `isomiga_gencode.gtf` | `/sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/isomiga_gencode.gtf` (isoMiGA v38: GENCODE v38 + 35,261 microglia PacBio transcripts; 272,340 transcripts) | project file, 2026-08-08 |
| `GRCh38.primary_assembly.genome.fa` | `/sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/GRCh38.primary_assembly.genome.fa` — the genome the cohort BAMs are aligned to | project file |

## What QAPA needs

Following https://github.com/morrislab/qapa:

1. A 3′ UTR library in BED format. Either QAPA's precompiled library or one
   built with `qapa build` from a genePred annotation plus poly(A) site
   annotations.
2. `qapa fasta`: 3′ UTR sequences extracted from the genome.
3. `salmon index` over those sequences.
4. An identifier table with columns `Gene stable ID`, `Transcript stable ID`,
   `Gene type`, `Transcript type`, `Gene name`, used by `qapa build` (to keep
   protein-coding transcripts) and `qapa quant` (to annotate events).

The three scripts under `scripts/` implement these steps:

| Script | Step |
|---|---|
| `scripts/qapa_identifiers_from_gtf.sh` | identifier table from a GTF (4) |
| `scripts/build_qapa_utr_library.sh` | genePred + poly(A) annotations → `qapa build` (1) |
| `scripts/prepare_qapa_reference.sh` | `qapa fasta` + `salmon index`, stages BED and table (2, 3) |

## Build 1: `hg38` — QAPA precompiled GENCODE V31 library

```bash
cd /sc/arion/projects/pd-omics/alex/pipelines/apa-pipeline
OUT=/sc/arion/projects/pd-omics/Mega_MyND/references/apa/hg38/qapa
mkdir -p "$OUT"

# Identifier table from the GENCODE release the BED was built from.
# QAPA documents a BioMart export; the GTF has the same five fields and is
# reproducible offline. Versions are stripped to match the BED's unversioned
# ENST IDs. PAR_Y duplicates are dropped.
curl -sSL -o "$OUT/gencode.v31.annotation.gtf.gz" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
bash scripts/qapa_identifiers_from_gtf.sh "$OUT/gencode.v31.annotation.gtf.gz" \
  > "$OUT/ensembl_identifiers.tsv"

# Sequences and Salmon index.
THREADS=4 bash scripts/prepare_qapa_reference.sh \
  references/qapa_3utrs.gencode_V31.hg38.bed.gz \
  "$OUT/ensembl_identifiers.tsv" \
  /sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/GRCh38.primary_assembly.genome.fa \
  "$OUT"
```

`prepare_qapa_reference.sh` runs, in order:

```bash
gzip -dc references/qapa_3utrs.gencode_V31.hg38.bed.gz > $OUT/qapa_3utrs.bed
# check: every transcript ID in the BED is present in the identifier table
qapa fasta -f GRCh38.primary_assembly.genome.fa $OUT/qapa_3utrs.bed $OUT/qapa_3utrs.fa
salmon index -t $OUT/qapa_3utrs.fa -i $OUT/salmon_index -p 4
```

Result: 42,034 BED entries (16,886 genes, 9,017 with more than one 3′ UTR
isoform); 41,856 sequences after `qapa fasta` skipped 178 (shorter than 100 nt
or duplicate sequence); identifier table 226,721 transcripts, all 42,034 BED
transcripts present.

## Build 2: `isomiga` — built from the isoMiGA v38 GTF

```bash
cd /sc/arion/projects/pd-omics/alex/pipelines/apa-pipeline
GTF=/sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/isomiga_gencode.gtf
OUT=/sc/arion/projects/pd-omics/Mega_MyND/references/apa/isomiga/qapa
mkdir -p "$OUT"

# Poly(A) annotations for qapa build.
curl -sSL -o references/wgEncodeGencodePolyaV38.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodePolyaV38.txt.gz
curl -sSL -o references/atlas.clusters.2.0.GRCh38.96.bed.gz \
  https://www.polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz

# Identifier table.
bash scripts/qapa_identifiers_from_gtf.sh "$GTF" > "$OUT/ensembl_identifiers.tsv"

# 3' UTR library.
THREADS=4 bash scripts/build_qapa_utr_library.sh "$GTF" "$OUT/ensembl_identifiers.tsv" \
  references/wgEncodeGencodePolyaV38.txt.gz \
  references/atlas.clusters.2.0.GRCh38.96.bed.gz \
  "$OUT/qapa_3utrs.bed"

# Sequences and Salmon index (same script as build 1).
THREADS=4 bash scripts/prepare_qapa_reference.sh "$OUT/qapa_3utrs.bed" \
  "$OUT/ensembl_identifiers.tsv" \
  /sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/GRCh38.primary_assembly.genome.fa \
  "$OUT"
```

`build_qapa_utr_library.sh` runs, in order:

```bash
# 1. GTF -> extended genePred (15 columns; qapa build requires -genePredExt)
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons isomiga_gencode.gtf raw.genePred
#    name2 (column 12) is gene_id by default; replaced with the gene symbol from
#    the identifier table so event IDs look like the precompiled library
#    (ENST..._SAMD11_hsa_chr1_...). Underscores in symbols become hyphens
#    because qapa quant splits event IDs on underscores. *_PAR_Y rows dropped.

# 2. UCSC genePred dump -> BED6 (chrom, start, end, class, score, strand);
#    qapa build keeps rows whose name is polyA_site
gzip -dc wgEncodeGencodePolyaV38.txt.gz | awk -F'\t' -v OFS='\t' '{print $3,$5,$6,$13,$12,$4}' > gencode_polya.bed
gzip -dc atlas.clusters.2.0.GRCh38.96.bed.gz > polyasite.bed

# 3. qapa build (defaults: PolyASite min score 3, merge 3' ends within 24 nt,
#    skip contigs not matching ^(chr)*[0-9XYM]+$)
qapa build --db ensembl_identifiers.tsv -g gencode_polya.bed -p polyasite.bed \
  -c 4 annotation.genePred > qapa_3utrs.bed
```

Result: 41,829 BED entries (16,790 genes, 8,931 with more than one 3′ UTR
isoform). 4,386 entries include a PacBio (`MSTRG`) transcript; 526 are 3′ UTR entries at 3′
ends supported only by PacBio transcripts. 41,676 sequences after `qapa fasta`
skipped 153. Identifier table 272,181 transcripts (86,709 GENCODE + 29,996
PacBio-on-known-gene + 407 PacBio-on-novel-gene pass QAPA's protein-coding
filter).

### Decisions specific to the isoMiGA build

- **PacBio transcripts have no `gene_type`, `transcript_type`, or
  `gene_name` attributes.** `qapa build` keeps a transcript only if both types
  are `protein_coding`. `qapa_identifiers_from_gtf.sh` assigns the gene type
  and name from the parent gene record (33,956 transcripts sit on ENSG genes)
  and transcript type `protein_coding` when the GTF annotates a CDS for the
  transcript, otherwise `processed_transcript`. `novelGene_*` genes (1,305
  transcripts, no gene record type) are `protein_coding` when any of their
  transcripts has a CDS; 407 qualify.
- **Version stripping only for `ENS*` IDs.** `qapa build` strips versions
  from Ensembl IDs before joining; `MSTRG.25.1` is joined verbatim, so its
  trailing `.1` is not a version and is kept.
- **Gene symbol as the BED gene field, underscores replaced.** `qapa quant`
  parses `<transcript>_<gene>_hsa_<chr>_...` and aborts on underscores inside
  those fields (e.g. `novelGene_1`, `5_8S_rRNA`). The `Gene` and `Gene_Name`
  columns of the event table come from the identifier table, unmodified.
- **PAR_Y transcripts dropped**: duplicates of chrX models, and their IDs
  contain underscores.
- **Non-standard contigs** (`GL000009.2`, `KI...`) are dropped by
  `qapa build` itself.
- **Novel-gene transcripts do not survive `qapa build`.** The 407
  `novelGene_*` transcripts pass the protein-coding filter and 194 of them
  yield an intron-free 3′ UTR on a standard chromosome, but QAPA's annotation
  step keeps only UTRs overlapping a GENCODE poly(A) site or PolyASite
  cluster, and just 1 of the 194 has a site within 24 nt of its 3′ end. These
  microglia-only intergenic transcripts are absent from the PAS databases, so
  the library contains no `novelGene_*` entries. PacBio transcripts on known
  genes are unaffected (526 library entries are PacBio-only 3′ ends with PAS support). Keeping
  unsupported ends would need `qapa build -N` (drops PAS filtering for every
  gene) or a custom site file (`-o`) that adds long-read 3′ ends as sites.
- **No Salmon decoys.** `qapa fasta --decoys` can make a genome-augmented
  index for selective alignment; the precompiled-library workflow and the QAPA
  paper use the plain 3′ UTR index, so both builds do the same. Revisit if
  off-target mapping to UTRs becomes a concern.

## Smoke test (2026-09-10)

Two samples (13-066-CC, 13-080-GFM) run through `qapa_salmon_quant` and
`qapa_events` with `qapa_input: fastq`, `reference_name: isomiga`,
`run_apalyzer: false` (`Mega_MyND/data/microglia/apa_smoke/config.yaml`).
Salmon: 21.6M and 19.1M fragments, 32.3% and 28.2% mapped to the 3′ UTR
index, ~2.3 min per sample on 16 cores. `qapa quant`: 41,676 events, 8,904
genes with more than one event, `Gene`/`Gene_Name` populated on every row,
524 events on PacBio-only 3′ ends, no ID parsing errors.

## Verification

For both builds every FASTA header maps to a BED entry, every sequence length
equals its BED interval, headers are unique, and every BED transcript ID is in
the identifier table. UTR length median 1,249 nt (hg38) and 1,401 nt
(isomiga), minimum 101 nt.

## Output layout

```text
/sc/arion/projects/pd-omics/Mega_MyND/references/apa/<reference_name>/qapa/
├── qapa_3utrs.bed             3' UTR library
├── qapa_3utrs.fa              qapa fasta output
├── ensembl_identifiers.tsv    config: references.<name>.qapa.identifiers
├── salmon_index/              config: references.<name>.qapa.salmon_index
└── salmon_index.log
```

---

# APAlyzer reference build log

Two APAlyzer references were built on 2026-09-11 under
`/sc/arion/projects/pd-omics/Mega_MyND/references/apa/apalyzer/<reference_name>/prepared_reference.RData`
by the pipeline rule `prepare_apalyzer_reference`
(`scripts/prepare_apalyzer_reference.R`, R 4.4.1, APAlyzer 1.20.0).

## What APAlyzer needs

Following https://github.com/RJWANGbioinfo/APAlyzer: three tables —
`refUTRraw` (per gene: proximal and distal 3′ UTR PAS plus CDS end, from
which aUTR/cUTR regions are made), `dfIPA` (intronic PAS with flanking splice
sites), and `dfLE` (3′-most exon per gene). `REF4PAS()` turns the raw tables
into the `UTRdbraw`, `dfIPA`, `dfLE` objects `PASEXP_3UTR()` and
`PASEXP_IPA()` count against. The raw tables come from either:

- an official reference RData file (poly(A) sites from PolyA_DB v3), or
- `PAS2GEF(GTFfile, AnnoMethod = "V2")`, which takes every annotated
  transcript 3′ end in a GTF as a PAS.

`prepare_apalyzer_reference.R` accepts both: a `.RData` source is loaded, a
`.gtf`/`.gtf.gz` source goes through `PAS2GEF()`. The config selects one with
`references.<name>.apalyzer.rdata` or `.gtf`.

## Inputs

| File | Source | md5 |
|---|---|---|
| `hg38_REF.RData` | https://github.com/RJWANGbioinfo/PAS_reference_RData/raw/master/hg38_REF.RData — kept at `Mega_MyND/references/apa/hg38/apalyzer/` | `2d7e9bc7d775ce54be7b78837fc7eed5` |
| `isomiga_gencode.pas2gef.gtf.gz` | `scripts/prepare_apalyzer_gtf.sh` on `isomiga_gencode.gtf` (see below) — kept at `Mega_MyND/references/apa/isomiga/apalyzer/` | `b70f9b6f238f412954b5123fad5f30bd` |

## Build 1: `hg38` — official PolyA_DB v3 reference

```bash
cd /sc/arion/projects/pd-omics/alex/pipelines/apa-pipeline
A=/sc/arion/projects/pd-omics/Mega_MyND/references/apa/hg38/apalyzer
curl -sSL -o $A/hg38_REF.RData https://github.com/RJWANGbioinfo/PAS_reference_RData/raw/master/hg38_REF.RData
ml R/4.4.1
Rscript scripts/prepare_apalyzer_reference.R $A/hg38_REF.RData hg38 \
  /sc/arion/projects/pd-omics/Mega_MyND/references/apa/apalyzer/hg38/prepared_reference.RData
```

Result: 11,550 genes with an aUTR/cUTR pair; 32,170 IPA sites in 9,308 genes;
9,228 last exons. Median aUTR 1,079 nt, cUTR 312 nt.

## Build 2: `isomiga` — PAS2GEF on the isoMiGA v38 GTF

```bash
cd /sc/arion/projects/pd-omics/alex/pipelines/apa-pipeline
A=/sc/arion/projects/pd-omics/Mega_MyND/references/apa/isomiga/apalyzer
bash scripts/prepare_apalyzer_gtf.sh \
  /sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/references/isomiga_gencode.gtf \
  $A/isomiga_gencode.pas2gef.gtf.gz

# Then, with references.isomiga.apalyzer.gtf pointing at that file, the
# pipeline rule (LSF, 4 cores, 64 GB requested; used 4.9 GB, 6 min):
snakemake -s Snakefile --configfile config.yaml --cluster-config cluster.yaml --cluster '...' \
  --jobs 1 /sc/arion/projects/pd-omics/Mega_MyND/references/apa/apalyzer/isomiga/prepared_reference.RData
# which runs
Rscript scripts/prepare_apalyzer_reference.R $A/isomiga_gencode.pas2gef.gtf.gz isomiga \
  /sc/arion/projects/pd-omics/Mega_MyND/references/apa/apalyzer/isomiga/prepared_reference.RData
```

`prepare_apalyzer_gtf.sh` exists because `PAS2GEF()` was written for Ensembl
GTFs. It:

- keeps chr1–22, X, Y, M and strips the `chr` prefix — `PAS2GEF()` prepends
  `chr` to every seqname, so a GENCODE-style GTF would yield `chrchr1` and
  count nothing against the BAMs;
- drops `*_PAR_Y` transcripts;
- adds `gene_name` to the isoMiGA PacBio (`MSTRG`) transcript/exon/CDS records
  from their gene record (`gene_id` for `novelGene_*`) — `PAS2GEF()` uses
  `gene_name` as a merge key between its PAS, intron, and last-exon tables,
  and rows with `NA` are silently lost.

`prepare_apalyzer_reference.R` attaches `GenomicRanges` before calling
`PAS2GEF()`: APAlyzer 1.20.0's `.annotatePASRegion()` calls `queryHits()`
without importing it, which fails with "could not find function" otherwise.

Result: 13,677 genes with an aUTR/cUTR pair (8,611 shared with the official
reference, 5,066 additional); 69,006 IPA sites in 16,460 genes; 16,473 last
exons. Median aUTR 484 nt, cUTR 241 nt. Annotated transcript ends are denser
and closer to the CDS than PolyA_DB's 3′-seq sites, so more genes qualify but
aUTRs are shorter; IPA sites here are intron-terminating annotated isoforms
rather than experimentally mapped intronic PAS.

## Reference comparison and decision (2026-09-11)

`PASEXP_3UTR()` on the same two BAMs with each reference:

| Reference | Genes | Both regions zero | aUTR zero, cUTR covered | Finite RE | Finite and cUTR ≥ 20 reads |
|---|---|---|---|---|---|
| `hg38` (PolyA_DB v3) | 11,550 | 12.5% | 13.5% | 69% | 6,120 / 6,114 |
| `isomiga` (PAS2GEF) | 13,722 | 22% | 26% | 47% | 4,652 / 4,628 |

Across the full cohort with `isomiga`, 53% of 3′ UTR RE values were
non-finite (18.5% gene not expressed, 26.5% aUTR empty with cUTR covered) and
4,916 of 13,677 genes had a finite RE in ≥ 80% of samples. `PAS2GEF()` takes
every annotated transcript end as a PAS, so for many genes the distal site is
a rarely used isoform end and the aUTR has no coverage; PolyA_DB sites are
experimentally supported at both ends.

Decision: `apalyzer_reference_name: hg38` (official reference) is the
default. APAlyzer output is keyed by gene symbol, so the isoform-level
consistency that motivated the QAPA `isomiga` library does not apply. The
`isomiga` APAlyzer cohort tables are kept as
`events/miga_isomiga_apalyzer_{3utr,ipa}_events.tsv`.

## Smoke test (2026-09-11)

Same two samples, `apalyzer_strand_type: NONE` (libraries are unstranded:
isoMiGA methods record `IU`, Salmon auto-detects `IU` with compatible
fragment ratio 1.0), `apalyzer_seq_type: ThreeMostPairEnd`, 16 threads.
~5 min and <3 GB per sample. Per sample: 13,686 3′ UTR genes (6.8M / 5.9M
reads in aUTR+cUTR), 69,260 / 69,266 IPA rows. `collect_apalyzer_events`
stacked them into the cohort tables with a `sample` column.
