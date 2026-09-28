# SUPPA data-generation methods

## Scope

This release generates SUPPA v2.4 local-event PSI and isoform PSI matrices for
the in-study isoMiGA microglia RNA-seq samples. It stops before phenotype
assembly, filtering, case-control modeling, contrasts, multiple-test correction,
or biological interpretation. Run-specific counts (samples, transcripts, events,
missingness) are in `qc/`, not in this document.

## Software and reproducibility

- SUPPA: official GitHub release tag `v2.4`, commit `cda335d`; `suppa.py`
  SHA-256 `9eb78e52ba0f25b5271c12d1126559bf8be343ccdff8d505e85d71b1f430db79`.
  Upstream left the embedded version string at 2.3, so `suppa.py --version`
  reports 2.3 despite the v2.4 release archive.
- Python 3.10.17 in `/sc/arion/projects/pd-omics/alex/conda/suppa2.4`; exact
  Conda package URLs in `suppa-pipeline/environment.explicit.txt`.
- Workflow: Snakemake 7 on Minerva LSF via the shared
  `pipelines/profiles/lsf` profile. Per-rule resources are declared in the
  Snakefile; each job's LSF log is under `logs/`.

## Step 1: reconcile Salmon expression and the annotation (`build_inputs`)

Inputs:

- `MiGA_isomiga_v38_annotation/sample_manifest.tsv`, restricted to
  `in_study == yes`, in manifest order. Its `salmon_quant` column gives one
  `quant.sf` per sample (the Gibbs-resampled Salmon 1.9.0 run under
  `salmon_gibbs/`; SUPPA uses only the point-estimate TPM).
- `aux_info/meta_info.json` beside each `quant.sf`, used to record mapping
  statistics and to require a single Salmon index across samples.
- `references/isomiga_gencode.gtf` and `references/tx2gene.isomiga.tsv`.

Processing:

- All `quant.sf` files are read in lockstep and must have identical, uniquely
  named transcript rows with numeric, non-negative TPM; per-sample TPM sums are
  recorded.
- Only terminal numeric versions on ENST and ENSG identifiers are removed.
  `MSTRG.<locus>.<isoform>` identifiers are kept unchanged.
- SUPPA reads exon records only, so an exon-only GTF is written, restricted to
  the Salmon transcript namespace. Annotation transcripts that Salmon collapsed
  or dropped when building its index are listed as `gtf_not_quantified` in
  `qc/annotation_transcript_audit.tsv`. GTF and tx2gene gene assignments must
  agree for every retained transcript.

Outputs: `expression/isomiga_transcript_tpm.tsv` (header = sample IDs; rows =
transcript ID then TPM), `reference/isomiga_suppa.gtf`,
`reference/transcript_id_map.tsv`, and `qc/quant_input_audit.tsv`,
`qc/annotation_transcript_audit.tsv`, `qc/input_build_summary.tsv`.

## Step 2: event and isoform definitions (`generate_events`, `generate_ioi`)

```text
suppa.py generateEvents -i isomiga_suppa.gtf -o isomiga -f ioe -e SE SS MX RI FL -b S -p
suppa.py generateEvents -i isomiga_suppa.gtf -o isomiga -f ioi -p
```

`-p` pools events sharing coordinates; `-b S` uses the strict boundary
definition. `SS` expands to A3 and A5 and `FL` to AF and AL, giving seven
classes: SE (skipped exon), A3/A5 (alternative 3'/5' splice site), MX (mutually
exclusive exons), RI (retained intron), AF/AL (alternative first/last exon).
Outputs: `events/isomiga_<TYPE>_strict.ioe` and matching `.gtf` tracks, and
`events/isomiga.ioi` mapping each transcript to its gene-level denominator.

## Step 3: local-event PSI (`psi_per_event`)

```text
suppa.py psiPerEvent -i isomiga_<TYPE>_strict.ioe -e isomiga_transcript_tpm.tsv -o isomiga_<TYPE>
```

PSI is the TPM of transcripts supporting the alternative form divided by the
TPM of all transcripts in the event. Outputs `psi/local/isomiga_<TYPE>.psi`
have the same sample header as the TPM matrix and the SUPPA event ID as first
column. A zero denominator gives a missing value, kept rather than imputed.

## Step 4: isoform PSI (`psi_per_isoform`)

```text
suppa.py psiPerIsoform -g isomiga_suppa.gtf -e isomiga_transcript_tpm.tsv -o isomiga
```

Transcript TPM divided by the summed TPM of its gene; one row per quantified
transcript, row ID `<gene>;<transcript>`, in `psi/isoform/isomiga_isoform.psi`.

## Step 5: validation and publication (`validate_release`, `publish`)

Validation requires:

- the in-study manifest sample order as the header of every matrix;
- unique transcript IDs in the TPM matrix, and every IOE/IOI transcript in it;
- exact row-set agreement between each IOE and its local PSI (SUPPA's row order
  is its own), and between the TPM, IOI and isoform PSI transcript sets;
- every observed PSI in [0, 1], missing values allowed.

Results: `qc/event_counts.tsv`, `qc/local_psi_summary.tsv`,
`qc/isoform_psi_summary.tsv`, and `qc/validation_passed.txt`. `publish` then
copies this document to `METHODS.md` and writes `SHA256SUMS.txt` over every
release file except `logs/`; check it with `sha256sum -c SHA256SUMS.txt` from
the release directory.
