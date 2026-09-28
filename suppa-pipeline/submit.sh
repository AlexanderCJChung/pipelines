#!/usr/bin/env bash
# Run the SUPPA pipeline on LSF. Extra arguments go to snakemake, e.g.
#   ./submit.sh -n                     dry run
#   ./submit.sh --until psi_per_event  stop after a rule
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
export XDG_CACHE_HOME="${PWD}/.snakemake/cache"
exec /sc/arion/projects/pd-omics/alex/conda/snakemake7/bin/snakemake \
  --profile ../profiles/lsf "$@"
