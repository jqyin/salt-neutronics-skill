#!/usr/bin/env bash
#
# run_on_hpc.sh -- portable launcher for the salt-neutronics TBR skill.
#
# Mimics a Shift Monte Carlo neutronics run on an HPC system by interpolating a
# precomputed parameter study. Everything after `--` (or any extra args) is
# forwarded verbatim to the CLI, so this script works for every subcommand.
#
# Usage:
#   bash scripts/run_on_hpc.sh tbr --bef2 33.3 --li6 0.075 --json
#   bash scripts/run_on_hpc.sh sweep --bef2-range 30 47 60 --plot
#
# To submit as a Slurm batch job, prepend SBATCH directives, e.g.:
#   #SBATCH -A <account>
#   #SBATCH -t 00:05:00
#   #SBATCH -N 1
# This is a lightweight, CPU-only, single-node post-processing task -- it does
# NOT need GPUs or multiple nodes.

set -euo pipefail

# --- Locate the skill root (parent of this script's directory) ----------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${SKILL_ROOT}"

ENV_NAME="${ENV_NAME:-salt-neutronics}"
ENV_YML="${ENV_YML:-environment.yml}"

# --- Optionally load site modules (OLCF Frontier / Andes etc.) ----------------
# Guarded so the script also runs on a laptop or any cluster without these.
if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null || true
    # Provide miniforge/conda; ignore failures on sites that name modules differently.
    module load miniforge3 2>/dev/null || module load python 2>/dev/null || true
fi

# --- Activate (or create) the conda environment, if conda is available --------
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
        echo "[run_on_hpc] Creating conda env '${ENV_NAME}' from ${ENV_YML} ..." >&2
        conda env create -f "${ENV_YML}"
    fi
    conda activate "${ENV_NAME}"
else
    echo "[run_on_hpc] conda not found; using system Python (numpy/scipy/h5py/matplotlib must be installed)." >&2
fi

# --- Run the requested query --------------------------------------------------
# Make the package importable without installation.
export PYTHONPATH="${SKILL_ROOT}:${PYTHONPATH:-}"
# Headless plotting backend (compute nodes have no display).
export MPLBACKEND="${MPLBACKEND:-Agg}"

# Strip a leading literal "--" if the caller used it as a separator.
if [[ "${1:-}" == "--" ]]; then shift; fi

if [[ "$#" -eq 0 ]]; then
    echo "[run_on_hpc] No arguments given; running the nominal-FLiBe example query." >&2
    set -- tbr --bef2 33.3 --li6 0.075
fi

# Prefer `python`, fall back to `python3` (some systems ship only the latter).
if command -v python >/dev/null 2>&1; then
    PYTHON=python
elif command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
else
    echo "[run_on_hpc] No python interpreter found on PATH." >&2
    exit 1
fi

echo "[run_on_hpc] ${PYTHON} -m salt_neutronics.cli $*" >&2
exec "${PYTHON}" -u -m salt_neutronics.cli "$@"
