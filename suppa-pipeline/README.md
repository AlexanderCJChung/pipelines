# SUPPA data-generation pipeline

SUPPA v2.4 local-event PSI and isoform PSI matrices from the isoMiGA Salmon
quantifications of the in-study MiGA microglia samples. No phenotype filtering
or differential analysis.

```
build_inputs ─┬─> generate_events ──> psi_per_event ×7 ─┐
 (TPM + GTF)  ├─> generate_ioi ─────────────────────────┼─> validate_release ─> publish
              └─> psi_per_isoform ──────────────────────┘
```

(`snakemake --rulegraph | dot -Tpdf` regenerates this from the Snakefile.)

## Inputs and outputs

- Samples: `/sc/arion/projects/pd-omics/MiGA_isomiga_v38_annotation/sample_manifest.tsv`
  (`in_study == yes`; `salmon_quant` column). Shared with the editing and APA pipelines.
- Annotation: `MiGA_isomiga_v38_annotation/references/isomiga_gencode.gtf` + `tx2gene.isomiga.tsv`.
- Release: `release_root` in `config.yaml`, under `MiGA_isomiga_v38_annotation/suppa/`.
  Holds `expression/ reference/ events/ psi/ qc/ logs/`, `METHODS.md` and `SHA256SUMS.txt`.
  Start a new dated `release_root` for a new release rather than rerunning into an old one.

`docs/SUPPA_DATA_METHODS.md` is the full methods text and is copied into each release.

## Running

```bash
./submit.sh -n                       # dry run
./submit.sh --until psi_per_event    # stop after a rule
./submit.sh                          # full run
```

Jobs go to LSF through the shared profile `../profiles/lsf`; resources are
declared per rule in the Snakefile and each job's log lands in `<release>/logs/`.
The whole run takes a few minutes of compute; PSI jobs peak at ~18 GB.
`tests/run_smoke.sh` runs SUPPA on a toy GTF locally.

## Vendored SUPPA

`vendor/SUPPA` (gitignored) is GitHub release tag `v2.4`, commit `cda335d`
(`suppa.py` SHA-256 `9eb78e52ba0f25b5271c12d1126559bf8be343ccdff8d505e85d71b1f430db79`).
Upstream did not update the embedded version string at that tag, so
`suppa.py --version` reports 2.3.
