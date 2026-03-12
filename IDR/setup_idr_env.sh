#!/usr/bin/env bash
# IDR pipeline – environment setup and launcher
#
# First run (or re-run after clean): creates a virtual environment, clones
# AIUPred and MobiDB-lite, installs all Python dependencies, and then runs
# the IDR analysis pipeline.
#
# Re-running is safe: existing venv/clones are reused.
#
# Usage:
#   bash IDR/setup_idr_env.sh          # setup + run pipeline
#   bash IDR/setup_idr_env.sh --setup-only  # setup without running pipeline
#
# After the first run you can also activate the environment manually:
#   source IDR/.venv/bin/activate
#   export AIUPRED_DIR="$PWD/IDR/tools/aiupred"
#   export MOBIDB_BINDIR="$PWD/IDR/tools/mobidb-lite/src/mobidb_lite/bin"
#   python IDR/idr_analysis.py

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"
TOOLS_DIR="$SCRIPT_DIR/tools"
AIUPRED_CLONE="$TOOLS_DIR/aiupred"
MOBIDB_CLONE="$TOOLS_DIR/mobidb-lite"
MOBIDB_BINDIR="$MOBIDB_CLONE/src/mobidb_lite/bin"
REQUIREMENTS="$SCRIPT_DIR/requirements.txt"

# Pinned to latest stable releases (update these when upgrading).
AIUPRED_TAG="2.1.2"
MOBIDB_TAG="v4.0"

SETUP_ONLY=true       # default: setup only – use run_pipeline.sh to run
for arg in "$@"; do
    [[ "$arg" == "--run" ]] && SETUP_ONLY=false
done

# ── helpers ───────────────────────────────────────────────────────────────────

print_header() { echo ""; echo "══ $* ══"; }
ok()  { echo "  ✓  $*"; }
info(){ echo "  →  $*"; }
warn(){ echo "  ⚠  $*"; }
fail(){ echo "  ✗  $*"; exit 1; }

# ── 1. Python version check ───────────────────────────────────────────────────

echo ""
echo "🧬  IDR Analysis Pipeline – Environment Setup"
echo "=============================================="

print_header "Python"
PYTHON=python3
if ! command -v "$PYTHON" &>/dev/null; then
    fail "python3 not found – please install Python 3.10 or later."
fi

PYVER=$("$PYTHON" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYMAJOR=$("$PYTHON" -c "import sys; print(sys.version_info.major)")
PYMINOR=$("$PYTHON" -c "import sys; print(sys.version_info.minor)")
ok "Found Python $PYVER"

if [[ "$PYMAJOR" -lt 3 || ( "$PYMAJOR" -eq 3 && "$PYMINOR" -lt 10 ) ]]; then
    fail "Python 3.10+ required (found $PYVER)."
fi

# ── 2. Virtual environment ────────────────────────────────────────────────────

print_header "Virtual environment"
if [[ -d "$VENV_DIR" ]]; then
    ok "Existing venv found at $VENV_DIR"
else
    info "Creating venv at $VENV_DIR …"
    "$PYTHON" -m venv "$VENV_DIR"
    ok "Venv created"
fi

# Activate
# shellcheck source=/dev/null
source "$VENV_DIR/bin/activate"
ok "Venv activated"

PIP="$VENV_DIR/bin/pip"
"$PIP" install -q --upgrade pip setuptools wheel
ok "pip/setuptools up to date"

# ── 3. Clone third-party tools ────────────────────────────────────────────────

print_header "Third-party tool repositories"
mkdir -p "$TOOLS_DIR"

# AIUPred
if [[ -d "$AIUPRED_CLONE/.git" ]]; then
    ok "AIUPred already cloned at $AIUPRED_CLONE"
else
    info "Cloning AIUPred $AIUPRED_TAG …"
    git clone --depth 1 --branch "$AIUPRED_TAG" \
        https://github.com/doszilab/AIUPred.git "$AIUPRED_CLONE"
    ok "AIUPred $AIUPRED_TAG cloned"
fi

# MobiDB-lite
if [[ -d "$MOBIDB_CLONE/.git" ]]; then
    ok "MobiDB-lite already cloned at $MOBIDB_CLONE"
else
    info "Cloning MobiDB-lite $MOBIDB_TAG …"
    git clone --depth 1 --branch "$MOBIDB_TAG" \
        https://github.com/BioComputingUP/MobiDB-lite.git "$MOBIDB_CLONE"
    ok "MobiDB-lite $MOBIDB_TAG cloned"
fi

# Make MobiDB-lite compiled sub-predictor binaries executable.
# We target only regular files within the immediate per-predictor subdirectories
# (DisEMBL, ESpritz, IUPred, SEG, TISEAN) to avoid making unrelated files executable.
if [[ -d "$MOBIDB_BINDIR" ]]; then
    find "$MOBIDB_BINDIR" -maxdepth 2 -type f ! -name "*.py" -exec chmod +x {} \;
    ok "MobiDB-lite binaries marked executable"
else
    warn "MobiDB-lite bin/ not found at $MOBIDB_BINDIR – check clone."
fi

# ── 4. Install Python dependencies ───────────────────────────────────────────

print_header "Python dependencies"

# Install MobiDB-lite from the local clone so bin/ is included in site-packages.
info "Installing MobiDB-lite from local clone (includes bin/ directory) …"
"$PIP" install -q "$MOBIDB_CLONE"
ok "MobiDB-lite installed"

# Install remaining deps (requirements.txt – the git+ MobiDB-lite line is
# commented out so we don't double-install from the remote URL).
info "Installing remaining dependencies from requirements.txt …"
"$PIP" install -q -r "$REQUIREMENTS"
ok "All dependencies installed"

# ── 5. Summary ────────────────────────────────────────────────────────────────

print_header "Environment ready"
echo ""
echo "  AIUPRED_DIR  = $AIUPRED_CLONE"
echo "  MOBIDB_BINDIR = $MOBIDB_BINDIR"
echo ""
echo "  To activate this environment in a new shell:"
echo "    source $VENV_DIR/bin/activate"
echo "    export AIUPRED_DIR=\"$AIUPRED_CLONE\""
echo "    export MOBIDB_BINDIR=\"$MOBIDB_BINDIR\""
echo ""

# ── 6. Run pipeline (unless --setup-only) ─────────────────────────────────────

if [[ "$SETUP_ONLY" == true ]]; then
    ok "Setup complete."
    echo ""
    echo "  Run the pipeline with:"
    echo "    bash IDR/run_pipeline.sh"
    echo "  or resume an interrupted run with:"
    echo "    bash IDR/run_pipeline.sh --resume"
    exit 0
fi

print_header "Running IDR pipeline"
export AIUPRED_DIR="$AIUPRED_CLONE"
export MOBIDB_BINDIR="$MOBIDB_BINDIR"

"$VENV_DIR/bin/python" "$SCRIPT_DIR/idr_analysis.py"
