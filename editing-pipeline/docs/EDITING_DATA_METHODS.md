# RNA editing: data-generation methods

A>I RNA editing calls for the in-study MiGA microglia samples (`in_study == yes` in
`MiGA_isomiga_v38_annotation/sample_manifest.tsv`), from the isoMiGA re-aligned BAMs
(STAR 2.7.3a against `isomiga_gencode.gtf`, coordinate-sorted, indexed, duplicates
flagged not removed, unstranded paired-end). No phenotype or QC filtering happens
here; sample QC and the analysis cohort are defined downstream in
`Mega_MyND/analyses/microglia/editing/editing_qc.Rmd`.

Ported from [RajLabMSSM/editing-pipeline](https://github.com/RajLabMSSM/editing-pipeline)
(Jack Humphrey, Winston Cuddleston, Flora Hyomin Seo, Raj Lab), upstream commit
`91d3898`, and repaired for Rocky 9, R 4.4.1 and cohorts of several hundred samples.
Thresholds are the values in the pipeline's `config.yaml`; run-specific counts
(sites, samples, substitution spectrum, true zeros, AEI) are in this release's `qc/`.

## Arms

| Arm | What it produces | Status |
|---|---|---|
| **JACUSA** | Per-site editing rate and coverage matrices across the cohort, ANNOVAR-annotated | Primary deliverable. On by default. |
| **AEI** | One global Alu Editing Index per sample | Cheap QC / covariate. On by default. |
| **REDItools** | Per-sample annotated call set | Exploratory. **Off by default** — see [REDItools caveats](#reditools-caveats). |

Toggle with `tools.run_*` in `config.yaml`. SAILOR was dropped: its singularity image
lives under `/sc/hydra/`, a filesystem that no longer exists.

## How the JACUSA arm works

Two passes, which is the part of the upstream design most worth keeping:

1. **Discovery.** Genome-wide de novo calling (`jacusa_denovo`) plus a supervised pass
   at REDIportal's 15.7M known sites (`jacusa_known`). The known pass uses a lower depth
   threshold, which is what recovers real editing that de novo filtering would miss.
   The union is filtered per-sample, merged per-chromosome, filtered across the cohort,
   and annotated with ANNOVAR.
2. **Back-fill.** `jacusa_pileup` re-queries *every* sample at the surviving site list.
   Without this the matrices would be presence-only: a site would be missing in a
   sample simply because that sample never nominated it.

Final outputs, in the release's `matrices/`:

| File | Contents |
|---|---|
| `isomiga_pileup_editing.tsv.gz` | sites x samples editing rate |
| `isomiga_pileup_coverage.tsv.gz` | sites x samples total read depth |
| `isomiga_pileup_dtu.tsv.gz` | two rows per site (ref / alt counts) for DRIMSeq-style testing |
| `isomiga_pileup_annotation.tsv.gz` | per-site gene, region, repeat, REDIportal and conservation annotation |
| `isomiga_aei.tsv` | per-sample Alu Editing Index |

Site IDs are `chr:pos:ref:alt` **oriented to the gene's strand**, so an A>I edit reads
as `A:G` regardless of which genomic strand the gene sits on.

## Cohort filter

All in-study samples form one group. A site enters annotation if it is detected in
at least `min_coverage` (40) samples with a mean editing rate of at least
`min_edit_rate` (0.01) across them. This is the relaxed filter of the 2026-09-13
production run: requiring 40 detections within a single diagnosis group (the
2026-08-20 setting) kept 2.4x fewer sites and made discovery depend on phenotype.
Per-group detection rules belong to the downstream QC.

## ANNOVAR is the runtime bottleneck, not your cohort size

`hg38_gnomad30_genome.txt` is **43 GB** and ANNOVAR scans it linearly. In a 374-site
test run the `annovar` rule was still going after 35 minutes with zero bytes written,
while every other rule had finished. The cost is in the database scan, so it is
essentially the same whether you query 400 sites or 400,000.

The `annovar` rule allows 3 h. If you do not need the population-frequency filter,
drop that protocol — see the `annovar_protocols` block in `config.yaml`:

```yaml
annovar_protocols:  refGene,dbsnp153CommonSNV,phastConsElements30way,rmsk,rediportal_012920
annovar_operations: g,f,r,r,f
annovar_args:       "'--neargene 10000',,,,"
```

`filter_annovar.R` detects the missing `AF` column and skips the gnomAD filter, leaving
common-dbSNP as the only SNP filter. That is looser than upstream, which is why the
full protocol list is still the default.

## Divergences from upstream

Upstream is at commit `91d3898` (Sept 2025, "Winston updated version of Jacusa
pipeline used in manuscript"). Everything below is a deliberate change.

### The one that changes your numbers

**`parse_pileup.R` no longer turns true zeros into missing data.**

Commit `91d3898` added `alt_cov > 0` to the pileup filter. The version before it was
explicit about why that must not happen:

```r
# rows where alt_cov=0 are kept here but considered by missingness filter later on
df_long <- filter(df_long, ref != alt)
```

With the filter in place, a sample with 200 clean reads and zero edited reads emits no
row for the locus, so the downstream `left_join` yields `NA` — indistinguishable from a
sample with no coverage. Editing rates are biased upward, and the site/sample
missingness filters discard well-measured unedited samples.

Restoring the old filter alone is not enough. Per-sample allele resolution is itself
wrong for a back-fill: the site list's alternative allele is already decided by the
cohort, so if sample X's top alt at a locus is G and sample Y's is C, upstream builds
different site IDs and the join drops both. With all-zero alt counts the "top" allele
is whichever sorts first — arbitrary.

So allele selection moved out of the per-sample parse entirely. `parse_pileup.R` emits
the full A/C/G/T count vector per locus; `merge_pileup.R` reads off the count for the
allele the annotation specifies. Every sample is asked the same question at the same
locus, and `NA` now means exactly one thing: fewer than `pileup_min_cov` reads.

The de novo discovery parse still resolves alleles — that is where allele discovery
belongs — and `multiallelic: drop|top` in config selects between the `91d3898` and
pre-`91d3898` behaviours there.

`jacusa_pileup` also passes `-c 1`. Upstream left JACUSA's default `-c 5` in place, so
a locus covered by 1–4 reads produced no line at all and the R side could not tell
"shallow" from "absent".

### Correctness

| | |
|---|---|
| `merge_pileup.R` | Upstream assigned flipped site IDs **positionally** (`coverage_DF$ESid <- anno_DF$ESid2`), correct only while both frames happened to share row order. Now a `match()`. |
| `filter_cohort.R` | Upstream round-tripped site IDs through `make.names` (`data.frame(t(...))` turns `chr1:12345:A:G` into `chr1.12345.A.G`, then `gsub("\\.", ":")` back) — lossy for any contig containing a dot, and two full matrix transposes. Now filtered in place; IDs are never used as R names. |
| `filter_cohort.R` | Grouping semantics made explicit. Upstream evaluated the coverage and editing-rate criteria in *separate* passes and intersected, so a site could qualify on coverage in group A and on editing rate in group B. Both must now hold within the same group; `group_logic: any\|all` decides how groups combine. |
| `filter_annovar.R` | Multi-gene ANNOVAR annotations are resolved rather than dropped. `Gene.refGene` is `"GENE1;GENE2"` for sites spanning two genes and `"NONE(dist=NONE),LINC01128(dist=1000)"` for intergenic sites with `--neargene`; upstream joined that raw string to GENCODE `gene_name`, so all of them failed to match and were discarded. |
| `concatenate_cohort.R` | Upstream's `rat_files` held the coverage paths and `cov_files` held the ratio paths, then swapped them back at the read — correct but misleading. The chromosome list also lived here *and* in the Snakefile; it now comes from config. |
| Snakefile | `rm -f` on JACUSA's `.filtered` sidecar (upstream's bare `rm` failed the rule whenever it was absent). `_pileup_dtu.tsv.gz` added to `rule all` — upstream produced it but never made it a build target. |

### Blockers that stopped it running at all

| | |
|---|---|
| `humandb` path | Upstream's config pointed at `.../references/editing/humandb/`, which does not exist. The real directory is `annovar_humandb/`. |
| R version | Upstream pinned `R/4.0.3`, which cannot see this account's user library. Now `R/4.4.1`. |
| `merge_AEI.R` | Called `write_tsv(path = ...)`; `path` was removed in readr 2.0. Also **no rule ever invoked it**, so the merged AEI table was never produced. A `merge_aei` rule now does. |
| REDItools scripts | The `.smk` referenced `scripts/filter_reditools.awk` and `scripts/filter_DP.py`. Those were renamed to `reditools_filter_strand.awk` / `reditools_filter_depth.py` and the `.smk` was never updated — the rule could never have run. There was also a stray `}` in the bedtools call. |
| REDItools config | Input/output folders, sample list, and thresholds were hardcoded Python literals, with `OUTFOLDER` assigned twice. Now config- and manifest-driven. |

### Cluster and performance

| | |
|---|---|
| Threads | Upstream hardcoded `jacusa_threads = 40` while its `cluster.yaml` reserved 4 cores. Each rule now declares `threads:` once and passes `{threads}` to the tool, so the reservation and the tool agree. |
| `-F` flag mask | Upstream used `-F 1024` (duplicates only). The MiGA BAMs are STAR with `--outSAMunmapped Within` and secondary alignments retained — the first record of `13-066-CC.bam` is FLAG 355. Default is now `3844` (unmapped + QC-fail + duplicate + secondary + supplementary). Set `flag_filter: 1024` to reproduce upstream. |
| Merge scaling | `reduce(cbind)` over N per-sample data frames reallocates the whole matrix once per sample — O(N²), and the reason those rules needed 100–150 GB of himem. Now a single indexed fill into a preallocated matrix. |
| Sample discovery | Upstream globbed the output directory, picking up stale files and making column order filesystem-dependent. The sample list and its order now come from the sample manifest (written to `qc/samples.tsv`), asserted at every merge. |
| Parse scripts | `parse_jacusa.R` and `parse_knownSites.R` were byte-identical apart from one default. Merged into one script; the threshold comes from config. Reimplemented in data.table — the tidyverse `pivot_longer` over a genome-wide JACUSA output dominated the discovery pass. |

## REDItools caveats

Off by default, and worth understanding before turning it on:

- `ml reditools` loads **python/2.7.16**, whose `hashlib` is broken against Rocky 9's
  OpenSSL. It prints a wall of `unsupported hash type` errors to stderr and still
  exits 0. The rules redirect that noise; a genuine failure is still caught by `pipefail`.
- OpenMPI warns about a missing `libfabric.so.1`; the rule forces `--mca btl tcp,self`.
- Upstream created `temp_cov/` and `temp_results/` in the submission CWD, so every
  concurrent sample collided. They now go to node-local scratch keyed on `$LSB_JOBID`
  and sample.
- **There is no cohort merge.** This arm produces per-sample annotated call sets and
  stops. It is a second opinion on JACUSA, not a parallel deliverable.

## Environment

Verified on Minerva 2026-08-04.

| | |
|---|---|
| `jacusa2/2.0.1` | `/hpc/packages/minerva-centos7/jacusa2/2.0.1/JACUSA_v2.0.1.jar`, pulls `java/1.8.0_211`. Upstream JACUSA is on 2.0.4; only 2.0.1 is packaged here, and it is what the manuscript used. |
| `rnaeditingindexer/2019-11-19` | Resources tree under `ad-omics/data/software/RNAEditingIndexer/Resources`. |
| `annovar/2019Oct24` | humandb has all six requested protocols. |
| `R/4.4.1` | optparse, data.table 1.18, tidyverse present. |
| Snakemake | 7, via the shared `pipelines/profiles/lsf` profile; resources per rule in the Snakefile. |
| MiGA BAMs | isoMiGA re-alignment: STAR 2.7.3a, `isomiga_gencode.gtf`, coordinate-sorted, MarkDuplicates-flagged (not removed), chr-prefixed GRCh38, unstranded paired-end. |

## Interpretation caveats

- **MiGA is unstranded**, so strand comes entirely from the gene-strand flip in
  `merge_pileup.R`. An A>G call in a gene with antisense overlap is ambiguous — the
  same caveat already recorded for APAlyzer in `apa-pipeline/config.yaml`.
- Sites with no resolvable GENCODE gene are dropped, because strand is what makes an
  A>I call interpretable and an unstranded library offers no other source for it.
- `max_gnomad_af: 0` requires absence from gnomAD entirely (upstream's behaviour).
  That is conservative: it also discards genuine editing sites carrying a singleton
  population call.


## Release checks

`validate_release` fails the release unless the matrix columns are the in-study samples
in manifest order (minus any removed by `sample_missingness`, which are listed), at
least `min_ag_fraction` (0.8) of final sites are gene-oriented A>G, the editing matrix
contains true zeros, and every in-study sample has an AEI. Its counts go to
`qc/release_summary.tsv` and `qc/substitution_spectrum.tsv`. Two checks remain manual:
the strongest sites should be in REDIportal, and the per-sample mean editing rate at
Alu sites should track the AEI.
