#!/usr/bin/env bash
# =============================================================================
# run_NuMI_nue.sh
#
# Convenience wrapper for submitting the NuMI nue selection to the grid using
# the medulla batch system.  All grid nodes will check out the branch specified
# by --tag (default: medulla_claude).
#
# Usage:
#   ./run_NuMI_nue.sh --project-dir=<pnfs_path>  [options]
#
# Modes (pick one; status check is always performed when the project exists):
#   --create              Create the project database (requires --batch-size).
#   --launch [N]          Launch all pending jobs, or N jobs if given.
#   --test                Launch a single test job.
#
# Options:
#   --project-dir=PATH    (required) dCache directory for this project.
#   --toml=FILE           Selection TOML file (default: NuMI_nue.toml).
#                         May be a filename relative to selection/toml/ or an
#                         absolute path.  Available NuMI nue configs:
#                           NuMI_nue.toml           CV extended sample only
#                           NuMI_nue_all.toml       All datasets (CV + data + det-sys)
#                           NuMI_nue_data.toml      On-beam + off-beam data only
#                           NuMI_nue_Null_CV.toml   Null / CV det-sys reference
#                           NuMI_nue_gain.toml      Gain ± variations
#                           NuMI_nue_coh.toml       Coherent noise ± variations
#                           NuMI_nue_int.toml       Interaction rate ± variations
#                           NuMI_nue_liftime.toml   Electron lifetime ± variations
#                           NuMI_nue_YZ_scint.toml  YZ uniformity + scintillation
#                           NuMI_nue_Wire_Recomb.toml  Wire gap + recombination
#                           NuMI_nue_Cath_Bend.toml Cathode bend variation
#   --batch-size=N        Files per grid job (default: 1).
#   --tag=BRANCH          medulla git branch/tag for grid nodes (default: medulla_claude).
#   --experiment=EXP      Experiment token for jobsub [icarus|sbnd] (default: icarus).
#   --memory=MB           Memory per job in MB (default: 1800).
#   --disk=GB             Disk per job in GB (default: 25).
#   --lifetime=DUR        Expected job wall time, e.g. 1h, 30m (default: 1h).
#   -h, --help            Print this message and exit.
#
# Examples:
#   # Create project then launch all jobs (CV sample)
#   ./run_NuMI_nue.sh --project-dir=/pnfs/icarus/scratch/users/$USER/numi_nue_run1 \
#       --create --batch-size=5
#
#   ./run_NuMI_nue.sh --project-dir=/pnfs/icarus/scratch/users/$USER/numi_nue_run1 \
#       --launch
#
#   # Run over all datasets
#   ./run_NuMI_nue.sh --project-dir=/pnfs/icarus/scratch/users/$USER/numi_nue_all \
#       --toml=NuMI_nue_all.toml --create --batch-size=5
#
#   # Single test job on a different branch with a specific TOML
#   ./run_NuMI_nue.sh --project-dir=/pnfs/icarus/scratch/users/$USER/numi_nue_test \
#       --toml=NuMI_nue_data.toml --create --batch-size=1 --tag=develop
#   ./run_NuMI_nue.sh --project-dir=/pnfs/icarus/scratch/users/$USER/numi_nue_test \
#       --test
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve the repo root so the script works from any working directory
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_DIR="${SCRIPT_DIR}/batch"
TOML_DIR="${SCRIPT_DIR}/selection/toml"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
PROJECT_DIR=""
TOML_ARG="NuMI_nue.toml"   # filename or absolute path; resolved below after parsing
BATCH_SIZE=1
TAG="medulla_claude"
EXPERIMENT="icarus"
MEMORY=1800
DISK=25
LIFETIME="1h"
MODE=""        # create | launch | test | (empty = status only)
NJOBS=""       # only used when MODE=launch

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
    sed -n '/^# Usage:/,/^# =\{5\}/p' "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project-dir=*)  PROJECT_DIR="${1#*=}"; shift ;;
        --project-dir)    PROJECT_DIR="$2"; shift 2 ;;
        --toml=*)         TOML_ARG="${1#*=}"; shift ;;
        --toml)           TOML_ARG="$2"; shift 2 ;;
        --batch-size=*)   BATCH_SIZE="${1#*=}"; shift ;;
        --batch-size)     BATCH_SIZE="$2"; shift 2 ;;
        --tag=*)          TAG="${1#*=}"; shift ;;
        --tag)            TAG="$2"; shift 2 ;;
        --experiment=*)   EXPERIMENT="${1#*=}"; shift ;;
        --experiment)     EXPERIMENT="$2"; shift 2 ;;
        --memory=*)       MEMORY="${1#*=}"; shift ;;
        --memory)         MEMORY="$2"; shift 2 ;;
        --disk=*)         DISK="${1#*=}"; shift ;;
        --disk)           DISK="$2"; shift 2 ;;
        --lifetime=*)     LIFETIME="${1#*=}"; shift ;;
        --lifetime)       LIFETIME="$2"; shift 2 ;;
        --create)         MODE="create"; shift ;;
        --launch)
            MODE="launch"
            # Optional inline job count: --launch 10  or  --launch=10
            if [[ -n "${2:-}" && "$2" =~ ^[0-9]+$ ]]; then
                NJOBS="$2"; shift 2
            else
                shift
            fi
            ;;
        --launch=*)       MODE="launch"; NJOBS="${1#*=}"; shift ;;
        --test)           MODE="test"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Resolve TOML path: absolute paths pass through; bare filenames are looked
# up inside selection/toml/.
# ---------------------------------------------------------------------------
if [[ "$TOML_ARG" = /* ]]; then
    TOML="$TOML_ARG"
else
    TOML="${TOML_DIR}/${TOML_ARG}"
fi

if [[ ! -f "$TOML" ]]; then
    echo "ERROR: TOML file not found: ${TOML}" >&2
    echo "       Available configs in ${TOML_DIR}/:" >&2
    ls "${TOML_DIR}"/NuMI_nue*.toml 2>/dev/null | xargs -n1 basename >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
if [[ -z "$PROJECT_DIR" ]]; then
    echo "ERROR: --project-dir is required." >&2
    usage >&2
    exit 1
fi

if [[ "$EXPERIMENT" != "icarus" && "$EXPERIMENT" != "sbnd" ]]; then
    echo "ERROR: --experiment must be 'icarus' or 'sbnd'." >&2
    exit 1
fi

if [[ "$MODE" == "create" && -z "$BATCH_SIZE" ]]; then
    echo "ERROR: --batch-size is required with --create." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build the base medulla.py command
# ---------------------------------------------------------------------------
PYTHON_CMD=(
    python3 "${BATCH_DIR}/medulla.py"
    --project-dir "${PROJECT_DIR}"
    --experiment  "${EXPERIMENT}"
    --tag         "${TAG}"
    --memory      "${MEMORY}"
    --disk        "${DISK}"
    --lifetime    "${LIFETIME}"
)

# ---------------------------------------------------------------------------
# Execute the requested mode
# ---------------------------------------------------------------------------
echo "[INFO] Repository : ${SCRIPT_DIR}"
echo "[INFO] TOML       : ${TOML}"
echo "[INFO] Project    : ${PROJECT_DIR}"
echo "[INFO] Branch/tag : ${TAG}"
echo "[INFO] Experiment : ${EXPERIMENT}"
echo ""

case "$MODE" in
    create)
        echo "[INFO] Creating project with batch size ${BATCH_SIZE}..."
        "${PYTHON_CMD[@]}" \
            --create-project \
            --toml       "${TOML}" \
            --batch-size "${BATCH_SIZE}"
        echo ""
        echo "[INFO] Project created.  Run with --launch to submit jobs."
        ;;

    launch)
        if [[ -n "$NJOBS" ]]; then
            echo "[INFO] Launching ${NJOBS} job(s)..."
            "${PYTHON_CMD[@]}" --launch-jobs "${NJOBS}"
        else
            echo "[INFO] Launching all pending jobs..."
            "${PYTHON_CMD[@]}" --launch-jobs
        fi
        ;;

    test)
        echo "[INFO] Launching single test job..."
        "${PYTHON_CMD[@]}" --test-job
        ;;

    *)
        # No mode flag — just report project status (requires the project to exist)
        echo "[INFO] Checking project status (no mode flag given)..."
        "${PYTHON_CMD[@]}"
        ;;
esac
