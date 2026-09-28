#!/usr/bin/env bash
# Run the RNA editing pipeline on LSF. Extra arguments go to snakemake, e.g.
#   ./submit.sh -n                     dry run
#   ./submit.sh --until filter_annovar discovery only, before the back-fill
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
export XDG_CACHE_HOME="${PWD}/.snakemake/cache"
exec /sc/arion/projects/pd-omics/alex/conda/snakemake7/bin/snakemake \
  --profile ../profiles/lsf "$@"
