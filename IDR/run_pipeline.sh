#!/usr/bin/env bash
# IDR Analysis Pipeline – Run Script
#
# Activates the venv created by setup_idr_env.sh and runs idr_analysis.py.
# By default the run resumes from any existing ta_idr_results.csv so that
# interrupted runs (e.g. due to a session timeout) pick up where they left off.
#
# Usage:
#   bash IDR/run_pipeline.sh                    # resume TA-only run if output exists
#   bash IDR/run_pipeline.sh --no-resume        # start a clean TA-only run from scratch
#   bash IDR/run_pipeline.sh --no-predictions   # DisProt annotation only, no AI tools
#   bash IDR/run_pipeline.sh --all-proteins     # whole dataset → IDR/all_idr_results.csv
#   bash IDR/run_pipeline.sh --output /path/to/out.csv
#
# Options are passed straight through to idr_analysis.py; see
#   python IDR/idr_analysis.py --help
# for the full list.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"
AIUPRED_CLONE="$SCRIPT_DIR/tools/aiupred"
MOBIDB_BINDIR="$SCRIPT_DIR/tools/mobidb-lite/src/mobidb_lite/bin"

# ── Sanity checks ─────────────────────────────────────────────────────────────

if [[ ! -d "$VENV_DIR" ]]; then
    echo "✗  Virtual environment not found at $VENV_DIR"
    echo "   Run setup first: bash IDR/setup_idr_env.sh"
    exit 1
fi

if [[ ! -d "$AIUPRED_CLONE" ]]; then
    echo "⚠  AIUPred clone not found at $AIUPRED_CLONE – disorder predictions may be skipped."
fi

# ── Activate venv and export tool paths ───────────────────────────────────────

# shellcheck source=/dev/null
source "$VENV_DIR/bin/activate"

export AIUPRED_DIR="$AIUPRED_CLONE"
export MOBIDB_BINDIR="$MOBIDB_BINDIR"

# ── Run ───────────────────────────────────────────────────────────────────────

echo ""
echo "🧬  IDR Analysis Pipeline"
echo "========================="
echo "  AIUPRED_DIR   = $AIUPRED_DIR"
echo "  MOBIDB_BINDIR = $MOBIDB_BINDIR"
echo ""

# Pass all arguments straight through to idr_analysis.py (--no-resume,
# --no-predictions, --output, etc.)
"$VENV_DIR/bin/python" "$SCRIPT_DIR/idr_analysis.py" "$@"
