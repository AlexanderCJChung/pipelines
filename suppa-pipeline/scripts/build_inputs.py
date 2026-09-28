#!/usr/bin/env python3
"""Build SUPPA expression and annotation inputs from Salmon quantifications.

The expression matrix is streamed in transcript-major order across every
quant.sf, so the 271,438 x 587 matrix is never materialized in memory. The
Samples and quant.sf paths come from the shared sample manifest. The
SUPPA GTF contains exon rows only (the only feature SUPPA reads) and is
restricted to transcripts present in the Salmon output namespace.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from contextlib import ExitStack
from pathlib import Path


ENST_VERSION = re.compile(r"(ENST\d+)\.\d+")
ENSG_VERSION = re.compile(r"(ENSG\d+)\.\d+")
TRANSCRIPT_ATTR = re.compile(r'(?:^|;\s*)transcript_id "([^"]+)"')
GENE_ATTR = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')


def normalize_transcript(value: str) -> str:
    return ENST_VERSION.sub(r"\1", value)


def normalize_gene(value: str) -> str:
    return ENSG_VERSION.sub(r"\1", value)


def normalize_attributes(value: str) -> str:
    return normalize_gene(normalize_transcript(value))


def read_manifest(path: Path) -> tuple[list[str], list[Path]]:
    """In-study samples and their quant.sf paths, in manifest order."""
    with path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["in_study"] == "yes"]
    if not rows:
        raise ValueError(f"No in-study samples in {path}")
    samples = [row["sample_id"] for row in rows]
    if len(set(samples)) != len(samples):
        raise ValueError(f"Duplicated sample IDs in {path}")
    return samples, [Path(row["salmon_quant"]) for row in rows]


def ensure_parents(paths: list[Path]) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def read_meta_info(quant_path: Path) -> dict[str, object]:
    meta_path = quant_path.parent / "aux_info" / "meta_info.json"
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing Salmon metadata: {meta_path}")
    with meta_path.open() as handle:
        return json.load(handle)


def build_expression(
    samples: list[str],
    quant_paths: list[Path],
    output_tpm: Path,
    sample_audit: Path,
    id_map: Path,
) -> tuple[list[str], list[str], list[float]]:
    missing = [str(path) for path in quant_paths if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing/empty quant.sf files:\n" + "\n".join(missing[:20]))

    audit_rows: list[dict[str, object]] = []
    for sample, quant_path in zip(samples, quant_paths):
        meta = read_meta_info(quant_path)
        audit_rows.append(
            {
                "sample": sample,
                "quant_sf": str(quant_path),
                "quant_bytes": quant_path.stat().st_size,
                "num_bootstraps": meta.get("num_bootstraps", ""),
                "num_processed": meta.get("num_processed", ""),
                "num_mapped": meta.get("num_mapped", ""),
                "percent_mapped": meta.get("percent_mapped", ""),
                "index_seq_hash": meta.get("index_seq_hash", ""),
                "index_name_hash": meta.get("index_name_hash", ""),
            }
        )

    index_hashes = {row["index_seq_hash"] for row in audit_rows}
    if len(index_hashes) != 1:
        raise ValueError(f"quant.sf files come from {len(index_hashes)} Salmon indexes")

    raw_ids: list[str] = []
    normalized_ids: list[str] = []
    tpm_sums = [0.0] * len(samples)

    with ExitStack() as stack:
        handles = [stack.enter_context(path.open()) for path in quant_paths]
        headers = [handle.readline().rstrip("\n").split("\t") for handle in handles]
        expected = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
        for sample, header in zip(samples, headers):
            if header != expected:
                raise ValueError(f"{sample}: unexpected quant.sf header {header}")

        with output_tpm.open("w", buffering=1024 * 1024) as output:
            output.write("\t".join(samples) + "\n")
            row_number = 0
            while True:
                lines = [handle.readline() for handle in handles]
                at_eof = [line == "" for line in lines]
                if all(at_eof):
                    break
                if any(at_eof):
                    ended = [samples[i] for i, value in enumerate(at_eof) if value]
                    raise ValueError(
                        f"quant.sf files have unequal row counts at row {row_number + 1}; "
                        f"early EOF: {', '.join(ended[:10])}"
                    )

                fields = [line.rstrip("\n").split("\t") for line in lines]
                if any(len(row) != 5 for row in fields):
                    raise ValueError(f"Malformed quant.sf row {row_number + 1}")
                raw_id = fields[0][0]
                mismatches = [samples[i] for i, row in enumerate(fields) if row[0] != raw_id]
                if mismatches:
                    raise ValueError(
                        f"Transcript order mismatch at row {row_number + 1} ({raw_id}); "
                        f"samples: {', '.join(mismatches[:10])}"
                    )

                normalized_id = normalize_transcript(raw_id)
                values = [row[3] for row in fields]
                try:
                    numeric = [float(value) for value in values]
                except ValueError as error:
                    raise ValueError(f"Non-numeric TPM for {raw_id}") from error
                if any(value < 0 for value in numeric):
                    raise ValueError(f"Negative TPM for {raw_id}")

                raw_ids.append(raw_id)
                normalized_ids.append(normalized_id)
                for i, value in enumerate(numeric):
                    tpm_sums[i] += value
                output.write(normalized_id + "\t" + "\t".join(values) + "\n")
                row_number += 1

    if len(set(raw_ids)) != len(raw_ids):
        raise ValueError("Duplicated raw transcript IDs in quant.sf")
    if len(set(normalized_ids)) != len(normalized_ids):
        duplicates = sorted(
            transcript for transcript in set(normalized_ids) if normalized_ids.count(transcript) > 1
        )
        raise ValueError(
            "Transcript version removal created duplicates: " + ", ".join(duplicates[:20])
        )

    for row, total in zip(audit_rows, tpm_sums):
        row["transcript_rows"] = len(normalized_ids)
        row["tpm_sum"] = f"{total:.6f}"

    with sample_audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]))
        writer.writeheader()
        writer.writerows(audit_rows)

    with id_map.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["salmon_transcript_id", "suppa_transcript_id", "origin"])
        for raw_id, normalized_id in zip(raw_ids, normalized_ids):
            origin = "GENCODE" if raw_id.startswith("ENST") else "PacBio_or_custom"
            writer.writerow([raw_id, normalized_id, origin])

    return raw_ids, normalized_ids, tpm_sums


def build_gtf(
    source_gtf: Path,
    quantified_ids: set[str],
    output_gtf: Path,
    annotation_audit: Path,
) -> tuple[int, int, int, dict[str, str]]:
    source_transcripts: set[str] = set()
    retained_transcripts: set[str] = set()
    transcript_gene: dict[str, str] = {}
    exon_rows = 0

    with source_gtf.open() as source, output_gtf.open("w", buffering=1024 * 1024) as output:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            transcript_match = TRANSCRIPT_ATTR.search(fields[8])
            gene_match = GENE_ATTR.search(fields[8])
            if transcript_match is None or gene_match is None:
                raise ValueError(f"Missing gene_id/transcript_id on exon line {line_number}")
            transcript_id = normalize_transcript(transcript_match.group(1))
            gene_id = normalize_gene(gene_match.group(1))
            source_transcripts.add(transcript_id)
            previous_gene = transcript_gene.get(transcript_id)
            if previous_gene is not None and previous_gene != gene_id:
                raise ValueError(
                    f"Transcript {transcript_id} maps to multiple genes: {previous_gene}, {gene_id}"
                )
            transcript_gene[transcript_id] = gene_id

            if transcript_id not in quantified_ids:
                continue
            fields[8] = normalize_attributes(fields[8])
            output.write("\t".join(fields) + "\n")
            retained_transcripts.add(transcript_id)
            exon_rows += 1

    missing_from_gtf = sorted(quantified_ids - retained_transcripts)
    if missing_from_gtf:
        raise ValueError(
            f"{len(missing_from_gtf)} quantified transcripts have no retained GTF exon; "
            f"examples: {', '.join(missing_from_gtf[:20])}"
        )

    gtf_only = sorted(source_transcripts - quantified_ids)
    with annotation_audit.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["suppa_transcript_id", "gene_id", "status"])
        for transcript_id in sorted(quantified_ids):
            writer.writerow([transcript_id, transcript_gene[transcript_id], "quantified_retained"])
        for transcript_id in gtf_only:
            writer.writerow([transcript_id, transcript_gene[transcript_id], "gtf_not_quantified"])

    return len(source_transcripts), len(retained_transcripts), exon_rows, transcript_gene


def check_tx2gene(
    path: Path,
    quantified_ids: set[str],
    gtf_map: dict[str, str],
) -> tuple[int, int]:
    seen: dict[str, str] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed tx2gene line {line_number}")
            transcript_id = normalize_transcript(fields[0])
            gene_id = normalize_gene(fields[1])
            if transcript_id not in quantified_ids:
                continue
            previous = seen.get(transcript_id)
            if previous is not None and previous != gene_id:
                raise ValueError(f"tx2gene maps {transcript_id} to multiple genes")
            seen[transcript_id] = gene_id

    missing = quantified_ids - set(seen)
    if missing:
        raise ValueError(
            f"{len(missing)} quantified transcripts absent from tx2gene; "
            f"examples: {', '.join(sorted(missing)[:20])}"
        )
    disagreements = [
        transcript_id
        for transcript_id, gene_id in seen.items()
        if gtf_map.get(transcript_id) != gene_id
    ]
    if disagreements:
        raise ValueError(
            f"{len(disagreements)} GTF/tx2gene disagreements; "
            f"examples: {', '.join(sorted(disagreements)[:20])}"
        )
    return len(seen), len(disagreements)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--tx2gene", type=Path, required=True)
    parser.add_argument("--output-tpm", type=Path, required=True)
    parser.add_argument("--output-gtf", type=Path, required=True)
    parser.add_argument("--sample-audit", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True)
    parser.add_argument("--annotation-audit", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    outputs = [
        args.output_tpm,
        args.output_gtf,
        args.sample_audit,
        args.id_map,
        args.annotation_audit,
        args.summary,
    ]
    ensure_parents(outputs)

    samples, quant_paths = read_manifest(args.manifest)
    raw_ids, normalized_ids, tpm_sums = build_expression(
        samples, quant_paths, args.output_tpm, args.sample_audit, args.id_map
    )
    quantified_ids = set(normalized_ids)
    source_n, retained_n, exon_rows, gtf_map = build_gtf(
        args.gtf, quantified_ids, args.output_gtf, args.annotation_audit
    )
    tx2gene_n, disagreements = check_tx2gene(
        args.tx2gene, quantified_ids, gtf_map
    )

    summary = {
        "samples": len(samples),
        "quantified_transcripts": len(normalized_ids),
        "source_gtf_transcripts": source_n,
        "gtf_transcripts_retained": retained_n,
        "gtf_transcripts_omitted": source_n - retained_n,
        "retained_exon_rows": exon_rows,
        "tx2gene_transcripts_checked": tx2gene_n,
        "tx2gene_disagreements": disagreements,
        "minimum_tpm_sum": f"{min(tpm_sums):.6f}",
        "maximum_tpm_sum": f"{max(tpm_sums):.6f}",
    }
    with args.summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary.items())

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise

