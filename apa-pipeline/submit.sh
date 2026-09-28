#!/usr/bin/env bash
#
# Submit the APA workflow to Minerva.
#
# CONFIG         configuration file (default: ./config.yaml)
# SNAKEMAKE_ENV  Snakemake 7 environment
# JOBS           maximum active jobs (default: 50)
# PHASES         phases to run (default: "1 2")
# DRYRUN         set to 1 to plan without submitting jobs

set -euo pipefail

pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${pipeline_dir}"

config_file="${CONFIG:-${pipeline_dir}/config.yaml}"
snakemake_env="${SNAKEMAKE_ENV:-/sc/arion/projects/pd-omics/alex/conda/snakemake7}"
jobs="${JOBS:-50}"
phases="${PHASES:-1 2}"
dry_run="${DRYRUN:-0}"

[ -f "${config_file}" ] || {
  echo "Config not found: ${config_file}" >&2
  exit 1
}

snakemake="${snakemake_env}/bin/snakemake"
[ -x "${snakemake}" ] || {
  echo "Snakemake not found: ${snakemake}" >&2
  exit 1
}

unset PYTHONPATH
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

snakemake_version="$("${snakemake}" --version)"
snakemake_major="${snakemake_version%%.*}"
if [ "${snakemake_major}" -ge 8 ]; then
  echo "This driver requires Snakemake 7.x. Found ${snakemake_version}." >&2
  exit 1
fi

read_config() {
  "${snakemake_env}/bin/python" - "$1" "$2" <<'PY'
import sys
import yaml

with open(sys.argv[1]) as handle:
    value = yaml.safe_load(handle)
for key in sys.argv[2].split("."):
    value = value.get(key, {}) if isinstance(value, dict) else {}
if isinstance(value, bool):
    print(str(value).lower())
elif isinstance(value, dict):
    print("")
else:
    print(value)
PY
}

qapa_input="$(read_config "${config_file}" qapa_input)"
run_qapa="$(read_config "${config_file}" tools.run_qapa)"
run_apalyzer="$(read_config "${config_file}" tools.run_apalyzer)"
output_dir="$(read_config "${config_file}" output_dir)"
reference_dir="$(read_config "${config_file}" reference_dir)"

[[ "${qapa_input}" == "bam" || "${qapa_input}" == "fastq" ]] || {
  echo "qapa_input must be bam or fastq" >&2
  exit 1
}
for setting in "${run_qapa}" "${run_apalyzer}"; do
  [[ "${setting}" == "true" || "${setting}" == "false" ]] || {
    echo "Tool settings must be true or false" >&2
    exit 1
  }
done
if [ "${run_qapa}" = "false" ] && [ "${run_apalyzer}" = "false" ]; then
  echo "No tools are enabled" >&2
  exit 1
fi
[ -n "${output_dir}" ] && [ -n "${reference_dir}" ] || {
  echo "output_dir and reference_dir are required" >&2
  exit 1
}

if [ "${dry_run}" != "1" ]; then
  mkdir -p logs "${output_dir}" "${reference_dir}"
fi

cluster_command='bsub -P acc_pd-omics -q {cluster.queue} -n {cluster.cores} -W {cluster.time} -R "rusage[mem={cluster.mem}] span[hosts=1]" -o {cluster.output} -e {cluster.error}'

run_phase() {
  local label="$1"
  local phase_jobs="$2"
  shift 2

  echo
  echo "${label}"
  if [ "${dry_run}" = "1" ]; then
    "${snakemake}" -s Snakefile --configfile "${config_file}" -n "$@"
  else
    "${snakemake}" -s Snakefile \
      --configfile "${config_file}" \
      --cluster-config cluster.yaml \
      --cluster "${cluster_command}" \
      --jobs "${phase_jobs}" \
      --latency-wait 60 \
      --rerun-incomplete \
      "$@"
  fi
}

has_phase() {
  [[ " ${phases} " == *" $1 "* ]]
}

add_target() {
  if [ "$1" = "true" ]; then
    targets+=("$2")
  fi
}

if has_phase 1; then
  targets=()
  add_target "${run_qapa}" qapa_salmon_quant
  add_target "${run_apalyzer}" apalyzer_count
  run_phase "Phase 1: sample quantification" "${jobs}" --until "${targets[@]}"
fi

if has_phase 2; then
  targets=()
  add_target "${run_qapa}" qapa_events
  add_target "${run_apalyzer}" collect_apalyzer_events
  run_phase "Phase 2: event tables" 10 --forcerun "${targets[@]}"
fi

echo
echo "Complete. Event tables: ${output_dir}/events/"
