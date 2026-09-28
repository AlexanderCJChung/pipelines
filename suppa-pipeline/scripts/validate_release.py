#!/usr/bin/env python3
"""Validate a SUPPA data release: sample headers, row membership, PSI range.

TPM values themselves were checked while build_inputs.py wrote the matrix; here
only its header and transcript IDs are read, as the reference for everything else.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Callable

MISSING = {"NA", "NaN", "nan", "NAN", ""}
IOE_HEADER = ["seqname", "gene_id", "event_id", "alternative_transcripts", "total_transcripts"]
IOI_HEADER = ["seqname", "gene_id", "isoform_id", "inclusion_transcripts", "total_transcripts"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--tpm", type=Path, required=True)
    parser.add_argument("--events", nargs="+", required=True)
    parser.add_argument("--ioe", type=Path, nargs="+", required=True, help="same order as --events")
    parser.add_argument("--local-psi", type=Path, nargs="+", required=True, help="same order as --events")
    parser.add_argument("--ioi", type=Path, required=True)
    parser.add_argument("--isoform-psi", type=Path, required=True)
    parser.add_argument("--event-counts", type=Path, required=True)
    parser.add_argument("--local-summary", type=Path, required=True)
    parser.add_argument("--isoform-summary", type=Path, required=True)
    parser.add_argument("--passed", type=Path, required=True)
    args = parser.parse_args()
    if not len(args.events) == len(args.ioe) == len(args.local_psi):
        parser.error("--events, --ioe and --local-psi need the same length")
    return args


def read_samples(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t") if row["in_study"] == "yes"]


def read_tpm_ids(path: Path, samples: list[str]) -> set[str]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader) != samples:
            raise ValueError(f"Sample header mismatch in {path}")
        transcripts = [row[0] for row in reader]
    if not transcripts or len(set(transcripts)) != len(transcripts):
        raise ValueError(f"Empty or duplicated transcript IDs in {path}")
    return set(transcripts)


def read_definitions(path: Path, header: list[str], id_col: int, transcripts: set[str]) -> list[str]:
    """IDs of an IOE/IOI file; every transcript it references must be quantified."""
    ids: list[str] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader) != header:
            raise ValueError(f"Unexpected header in {path}")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Wrong field count in {path}:{line_number}")
            referenced = {t for field in row[3:5] for t in field.split(",") if t}
            if referenced - transcripts:
                raise ValueError(f"{path}:{line_number} references unquantified {sorted(referenced - transcripts)[0]}")
            ids.append(row[id_col])
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicated IDs in {path}")
    return ids


def validate_psi(
    path: Path, samples: list[str], expected: set[str], key: Callable[[str], str] = lambda x: x
) -> dict[str, object]:
    """Header, exact row membership (SUPPA's row order is its own) and PSI in [0, 1]."""
    seen: set[str] = set()
    observed = missing = 0
    lo, hi = math.inf, -math.inf
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader) != samples:
            raise ValueError(f"Sample header mismatch in {path}")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(samples) + 1:
                raise ValueError(f"Wrong field count in {path}:{line_number}")
            identifier = key(row[0])
            if identifier in seen:
                raise ValueError(f"Duplicated row {identifier} in {path}")
            seen.add(identifier)
            for token in row[1:]:
                if token in MISSING:
                    missing += 1
                    continue
                value = float(token)
                if not math.isfinite(value) or value < -1e-12 or value > 1 + 1e-12:
                    raise ValueError(f"Invalid PSI {token} in {path}:{line_number}")
                observed += 1
                lo, hi = min(lo, value), max(hi, value)
    if seen != expected:
        raise ValueError(f"{path}: rows differ from definitions (observed {len(seen)}, expected {len(expected)})")
    return {
        "rows": len(seen),
        "samples": len(samples),
        "observed_values": observed,
        "missing_values": missing,
        "missing_fraction": missing / (observed + missing) if observed + missing else 0.0,
        "minimum_observed_psi": lo if observed else "NA",
        "maximum_observed_psi": hi if observed else "NA",
    }


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    samples = read_samples(args.manifest)
    transcripts = read_tpm_ids(args.tpm, samples)

    event_rows, local_rows = [], []
    for event_type, ioe, psi in zip(args.events, args.ioe, args.local_psi):
        events = read_definitions(ioe, IOE_HEADER, 2, transcripts)
        event_rows.append({"event_type": event_type, "events": len(events)})
        local_rows.append({"event_type": event_type, **validate_psi(psi, samples, set(events))})

    isoforms = read_definitions(args.ioi, IOI_HEADER, 3, transcripts)
    if set(isoforms) != transcripts:
        raise ValueError("IOI transcript set does not match the TPM transcript set")
    # psiPerIsoform row IDs are "<gene>;<transcript>"
    isoform_row = validate_psi(args.isoform_psi, samples, transcripts, key=lambda x: x.rsplit(";", 1)[-1])

    write_tsv(args.event_counts, event_rows)
    write_tsv(args.local_summary, local_rows)
    write_tsv(args.isoform_summary, [isoform_row])
    args.passed.write_text(
        "PASS\n"
        f"samples\t{len(samples)}\n"
        f"transcripts\t{len(transcripts)}\n"
        f"local_events\t{sum(row['events'] for row in event_rows)}\n"
    )


if __name__ == "__main__":
    main()
