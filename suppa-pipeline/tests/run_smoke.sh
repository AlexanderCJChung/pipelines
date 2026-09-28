#!/usr/bin/env bash
set -euo pipefail

ROOT=/sc/arion/projects/pd-omics/alex/pipelines/suppa-pipeline
PYTHON=/sc/arion/projects/pd-omics/alex/conda/suppa2.4/bin/python
SUPPA="${ROOT}/vendor/SUPPA/suppa.py"

tmpdir=$(mktemp -d /tmp/suppa_smoke.XXXXXX)
trap 'rm -rf "${tmpdir}"' EXIT

"${PYTHON}" "${SUPPA}" generateEvents \
  -i "${ROOT}/tests/fixtures/smoke.gtf" \
  -o "${tmpdir}/smoke" -f ioe -e SE SS MX RI FL -b S -p

test -s "${tmpdir}/smoke_SE_strict.ioe"

"${PYTHON}" "${SUPPA}" psiPerEvent \
  -i "${tmpdir}/smoke_SE_strict.ioe" \
  -e "${ROOT}/tests/fixtures/smoke.tpm" \
  -o "${tmpdir}/smoke_SE"

"${PYTHON}" "${SUPPA}" psiPerIsoform \
  -g "${ROOT}/tests/fixtures/smoke.gtf" \
  -e "${ROOT}/tests/fixtures/smoke.tpm" \
  -o "${tmpdir}/smoke"

"${PYTHON}" - "${tmpdir}/smoke_SE.psi" "${tmpdir}/smoke_isoform.psi" <<'PY'
import math
import sys

for path in sys.argv[1:]:
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        assert header == ["sample_a", "sample_b"], (path, header)
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    assert rows, path
    for row in rows:
        for value in row[1:]:
            x = float(value)
            assert math.isnan(x) or 0.0 <= x <= 1.0, (path, row)
print("SUPPA smoke test passed")
PY

