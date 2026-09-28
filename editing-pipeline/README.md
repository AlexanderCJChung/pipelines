# RNA editing data-generation pipeline

A>I RNA editing for the in-study MiGA microglia samples: JACUSA2 site-level editing
matrices and the Alu Editing Index, from the isoMiGA re-aligned BAMs. No phenotype or
QC filtering; that happens downstream in `Mega_MyND/analyses/microglia/editing/`.

```
sample_sheet ─┬─> jacusa_denovo ─> parse_denovo ─┐
              ├─> jacusa_known ──> parse_known ──┴> cat_all_sites ×N ─> merge_jacusa ×chrom
              │     ─> filter_cohort ×chrom ─> concatenate_cohort ─> annovar ─> filter_annovar
              │     ─> jacusa_pileup ×N ─> parse_pileup ×N ─> merge_pileup ─┐
              └─> aei_subset_bam ×N ─> aei ×N ─> merge_aei ─────────────────┴> validate_release ─> publish
REDItools (off): reditools ─> reditools_filter ─> reditools_annovar, per sample
```

(`snakemake --rulegraph | dot -Tpdf` regenerates this from the Snakefile.)

## Inputs and outputs

- Samples: `/sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/sample_manifest.tsv`
  (`in_study == yes`; `bam` column). Shared with the SUPPA and APA pipelines.
- References: paths in `config.yaml` (GRCh38 primary assembly, ENCODE blacklist,
  REDIportal, ANNOVAR humandb, GENCODE v38 gene strand, RNAEditingIndexer resources).
- Release: `release_root` in `config.yaml`, under `MiGA_isomiga_v38_annotation/editing/`.
  Holds `matrices/ qc/ cohort/ samples/ logs/`, `METHODS.md` and `SHA256SUMS.txt`
  (checksums cover everything except `samples/` and `logs/`). Per-sample JACUSA calls in
  `samples/` are ~70 GB and are kept so a later cohort-filter change can reuse them.
  Start a new dated `release_root` for a new release rather than rerunning into an old one.

`docs/EDITING_DATA_METHODS.md` is the full methods text, including every divergence
from upstream, and is copied into each release.

## Running

```bash
./submit.sh -n                        # dry run
./submit.sh --until filter_annovar    # discovery only, before the back-fill
./submit.sh                           # full run
./submit.sh --configfile tests/smoke_config.yaml   # smoke test, ~1.5 h (ANNOVAR)
```

Jobs go to LSF through the shared profile `../profiles/lsf`; resources are declared per
rule in the Snakefile and each job's log lands in `<release>/logs/`. A full run is about
1,800 JACUSA jobs (8 cores, up to ~35 min each) plus a ~1 h ANNOVAR step.

The smoke test runs both arms on 4 BAMs cut to chr21:33-37 Mb
(`/sc/arion/projects/pd-omics/alex/editing_smoke_tmp`) into scratch. JACUSA cannot run
on a login node (the JVM exhausts `ulimit -u`), so always submit.
