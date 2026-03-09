#!/usr/bin/env python3
"""
IDR Analysis Pipeline for Tail-Anchored (TA) Proteins

Filters the main membrane-protein dataset to genuine TA proteins, cross-references
them against human DisProt entries, and runs intrinsic-disorder predictions using
the AIUPred and MobiDB-lite standalone tools installed locally.

Tool setup
----------
AIUPred (https://github.com/doszilab/AIUPred):
    git clone https://github.com/doszilab/AIUPred.git /path/to/aiupred
    pip install torch>=2.6.0 numpy scipy
    export AIUPRED_DIR=/path/to/aiupred

MobiDB-lite (https://github.com/BioComputingUP/MobiDB-lite):
    pip install git+https://github.com/BioComputingUP/MobiDB-lite.git
    # Obtain compiled sub-predictor binaries and set:
    export MOBIDB_BINDIR=/path/to/mobidb-lite/bin

Outputs a DataFrame with all TA-protein columns plus appended disorder annotations.
"""

import io
import os
import subprocess
import sys
import tempfile
import time

import pandas as pd
import requests

# ── Tunable parameters ────────────────────────────────────────────────────────
# A TA protein has exactly one transmembrane domain (TMD), located towards
# the C-terminus with minimal sequence extension after it.

EXACT_TMD_COUNT: int = 1       # Total membrane domains (TMD + intramembrane)
MAX_N_TERM_DOMAINS: int = 0    # N-terminal membrane domains (0 = C-terminal only)
MAX_CTERM_EXTENSION: int = 30  # Max residues between TMD end and protein C-terminus

AIUPRED_THRESHOLD: float = 0.5   # Per-residue score ≥ this → predicted disordered
API_DELAY: float = 0.3           # Seconds to pause between UniProt sequence fetches

HUMAN_TAXON_ID: str = "9606"

# ── Local tool paths (override via environment variables) ─────────────────────
# Directory of the cloned doszilab/AIUPred repository.
AIUPRED_DIR: str = os.environ.get("AIUPRED_DIR", "")
# Directory containing the compiled MobiDB-lite binary sub-predictors.
MOBIDB_BINDIR: str = os.environ.get("MOBIDB_BINDIR", "")

# Default file paths (relative to this script's location)
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_PATH: str = os.path.join(_HERE, "..", "raw_data",
                               "membrane_protein_analysis_with_reduced_cc.csv")
DISPROT_PATH: str = os.path.join(_HERE, "..", "DisProt_release_2025_12.tsv")


# ── TA filtering ──────────────────────────────────────────────────────────────

def filter_ta_proteins(
    df: pd.DataFrame,
    tmd_count: int = EXACT_TMD_COUNT,
    max_n_term: int = MAX_N_TERM_DOMAINS,
    max_cterm: int = MAX_CTERM_EXTENSION,
) -> pd.DataFrame:
    """Return proteins that match tail-anchored (TA) criteria.

    A TA protein has:
    - Exactly *tmd_count* total membrane domains.
    - At most *max_n_term* domains situated N-terminally (i.e. the single TMD
      is towards the C-terminus).
    - A C-terminal tail of at most *max_cterm* residues after the TMD end.

    Parameters
    ----------
    df:
        The main membrane-protein DataFrame.
    tmd_count:
        Required total membrane domain count.
    max_n_term:
        Maximum allowed N-terminal membrane domains.
    max_cterm:
        Maximum allowed residues between the TMD end and the protein C-terminus.

    Returns
    -------
    pd.DataFrame
        Filtered copy containing only TA proteins.
    """
    working = df.copy()
    for col in ("membrane_domain_count", "N_term_md", "cterm_distance"):
        working[col] = pd.to_numeric(working[col], errors="coerce")

    mask = (
        (working["membrane_domain_count"] == tmd_count)
        & (working["N_term_md"] <= max_n_term)
        & (working["cterm_distance"] <= max_cterm)
    )
    return working[mask].copy()


# ── DisProt lookup ────────────────────────────────────────────────────────────

def filter_disprot_human(disprot_df: pd.DataFrame) -> pd.DataFrame:
    """Return only human entries from a DisProt DataFrame.

    Filters on 'NCBI Taxon ID' == HUMAN_TAXON_ID (9606).
    """
    return disprot_df[
        disprot_df["NCBI Taxon ID"].astype(str).str.strip() == HUMAN_TAXON_ID
    ].copy()


def match_disprot(ta_df: pd.DataFrame, disprot_human: pd.DataFrame) -> pd.DataFrame:
    """Annotate *ta_df* with matching DisProt entries.

    Joins on UniProt accession (``Entry`` in *ta_df*, ``UniProt ACC`` in
    *disprot_human*).  Multiple DisProt regions for the same protein are
    aggregated: the first DisProt ID / disorder-content value is kept and all
    Region IDs are joined with ``';'``.

    Adds columns
    ------------
    disprot_id : str or NaN
    disprot_disorder_content : float or NaN
    disprot_disorder_regions : str or NaN
        Semicolon-separated list of Region IDs.
    """
    agg = (
        disprot_human.groupby("UniProt ACC")
        .agg(
            disprot_id=("DisProt ID", "first"),
            disprot_disorder_content=("Protein Disorder Content", "first"),
            disprot_disorder_regions=("Region ID", lambda s: ";".join(s.dropna().unique())),
        )
        .reset_index()
        .rename(columns={"UniProt ACC": "_disprot_acc"})
    )

    merged = ta_df.merge(agg, left_on="Entry", right_on="_disprot_acc", how="left")
    merged.drop(columns=["_disprot_acc"], inplace=True)
    return merged


# ── UniProt sequence fetch ────────────────────────────────────────────────────

def fetch_sequence(uniprot_acc: str, delay: float = API_DELAY) -> str | None:
    """Download the protein FASTA sequence from the UniProt REST API.

    Returns the sequence string (no header) or ``None`` on failure.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        return "".join(line for line in lines if not line.startswith(">"))
    except Exception:
        return None
    finally:
        time.sleep(delay)


# ── AIUPred (local) ───────────────────────────────────────────────────────────

def run_aiupred_local(sequence: str, aiupred_dir: str = AIUPRED_DIR) -> dict:
    """Run AIUPred disorder prediction using the locally cloned standalone tool.

    Clone the repository and set ``AIUPRED_DIR`` (or pass *aiupred_dir*):

    .. code-block:: bash

        git clone https://github.com/doszilab/AIUPred.git /path/to/aiupred
        pip install torch>=2.6.0 numpy scipy
        export AIUPRED_DIR=/path/to/aiupred

    The function writes a temporary FASTA file, calls ``aiupred.py`` as a
    subprocess with ``--force-cpu``, and parses the tab-delimited output.

    Returns
    -------
    dict
        ``{"scores": [<float>, ...]}`` with one score per residue, or ``{}``
        when *aiupred_dir* is unset, the script is not found, or the call fails.
    """
    if not aiupred_dir or not os.path.isdir(aiupred_dir):
        return {}

    script = os.path.join(aiupred_dir, "aiupred.py")
    if not os.path.isfile(script):
        return {}

    fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(f">query\n{sequence}\n")

        result = subprocess.run(
            [sys.executable, script, "-i", fasta_path, "--force-cpu"],
            capture_output=True,
            text=True,
            timeout=300,
            cwd=aiupred_dir,
        )
        if result.returncode != 0:
            return {}

        # Output: lines starting with '#' are comments; data lines are
        # tab-separated "Position\tResidue\tDisorder"
        scores = []
        for line in result.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    scores.append(float(parts[2]))
                except ValueError:
                    continue
        return {"scores": scores}
    except (subprocess.TimeoutExpired, OSError):
        return {}
    finally:
        os.unlink(fasta_path)


def aiupred_summary(scores: list, threshold: float = AIUPRED_THRESHOLD) -> dict:
    """Summarise per-residue AIUPred *scores* into scalar metrics.

    Returns
    -------
    dict with keys:
    ``aiupred_mean_score``         – mean disorder score across all residues.
    ``aiupred_disordered_fraction``– fraction of residues ≥ *threshold*.
    Values are ``None`` when *scores* is empty.
    """
    if not scores:
        return {"aiupred_mean_score": None, "aiupred_disordered_fraction": None}
    mean_score = sum(scores) / len(scores)
    disordered_fraction = sum(s >= threshold for s in scores) / len(scores)
    return {
        "aiupred_mean_score": round(mean_score, 4),
        "aiupred_disordered_fraction": round(disordered_fraction, 4),
    }


# ── MobiDB-lite (local) ───────────────────────────────────────────────────────

def run_mobidb_lite_local(sequence: str, bindir: str = MOBIDB_BINDIR) -> dict:
    """Run MobiDB-lite disorder prediction using the locally installed package.

    Install the package and provide the path to its compiled binary tools:

    .. code-block:: bash

        pip install git+https://github.com/BioComputingUP/MobiDB-lite.git
        # place compiled sub-predictor binaries in a directory then:
        export MOBIDB_BINDIR=/path/to/mobidb-lite/bin

    The function uses the ``mobidb_lite.consensus.run`` Python API directly,
    passing a StringIO FASTA handle and the *bindir* path.

    Returns
    -------
    dict
        ``{"prediction-disorder-mobidb_lite": {"regions": [[start, end], ...]}}``
        compatible with :func:`mobidb_lite_summary`, or ``{}`` when *bindir* is
        unset, the package is not installed, or the prediction fails.
    """
    if not bindir or not os.path.isdir(bindir):
        return {}

    try:
        from mobidb_lite.consensus import run as _mobi_run  # type: ignore[import]
    except ImportError:
        return {}

    fasta_handle = io.StringIO(f">query\n{sequence}\n")
    try:
        for _seq_id, regions, _scores in _mobi_run(fasta_handle, bindir, "", threads=1):
            if regions is None:
                return {}
            mobi_regions = regions.get("mobidblite", [])
            return {"prediction-disorder-mobidb_lite": {"regions": mobi_regions}}
    except Exception:
        return {}

    return {}


def mobidb_lite_summary(data: dict, protein_length: int) -> dict:
    """Parse a MobiDB-lite result dict into a disordered-fraction metric.

    Returns
    -------
    dict with key:
    ``mobidb_lite_disordered_fraction``– fraction of residues in disordered
      regions, or ``None`` when the prediction is absent or length is zero.
    """
    if not data or not protein_length:
        return {"mobidb_lite_disordered_fraction": None}

    pred = data.get("prediction-disorder-mobidb_lite", {})
    regions = pred.get("regions", []) if isinstance(pred, dict) else []

    if not regions:
        return {"mobidb_lite_disordered_fraction": None}

    disordered_residues = sum(
        int(end) - int(start) + 1 for start, end in regions
    )
    fraction = min(disordered_residues / protein_length, 1.0)
    return {"mobidb_lite_disordered_fraction": round(fraction, 4)}


# ── Full prediction loop ──────────────────────────────────────────────────────

def run_disorder_predictions(ta_df: pd.DataFrame) -> pd.DataFrame:
    """Add AIUPred and MobiDB-lite predictions to each row of *ta_df*.

    Fetches the amino-acid sequence from UniProt for each protein, then runs
    both local predictors.  Failed calls produce ``None`` values rather than
    raising exceptions.

    Returns a new DataFrame with three additional columns appended:
    ``aiupred_mean_score``, ``aiupred_disordered_fraction``,
    ``mobidb_lite_disordered_fraction``.
    """
    records = []
    total = len(ta_df)
    for i, (_, row) in enumerate(ta_df.iterrows(), start=1):
        acc = row["Entry"]
        length = int(row["Length"]) if pd.notna(row.get("Length")) else 0
        print(f"  [{i}/{total}] {acc}", end=" ", flush=True)

        # Fetch sequence once; used by both local tools
        sequence = fetch_sequence(acc)
        if sequence:
            ai_raw = run_aiupred_local(sequence)
            ai_metrics = aiupred_summary(ai_raw.get("scores", []))

            mobi_raw = run_mobidb_lite_local(sequence)
            mobi_metrics = mobidb_lite_summary(mobi_raw, length)
        else:
            ai_metrics = {"aiupred_mean_score": None, "aiupred_disordered_fraction": None}
            mobi_metrics = {"mobidb_lite_disordered_fraction": None}

        records.append({**ai_metrics, **mobi_metrics})
        print("✓")

    pred_df = pd.DataFrame(records, index=ta_df.index)
    return pd.concat([ta_df, pred_df], axis=1)


# ── Main entry point ──────────────────────────────────────────────────────────

def main(
    data_path: str = DATA_PATH,
    disprot_path: str = DISPROT_PATH,
    run_predictions: bool = True,
    tmd_count: int = EXACT_TMD_COUNT,
    max_n_term: int = MAX_N_TERM_DOMAINS,
    max_cterm: int = MAX_CTERM_EXTENSION,
) -> pd.DataFrame:
    """Run the full IDR analysis pipeline.

    1. Load the membrane-protein DataFrame and the DisProt TSV.
    2. Filter to TA proteins using the supplied (or default) thresholds.
    3. Cross-reference against human DisProt entries.
    4. Optionally run AIUPred and MobiDB-lite predictions using the local
       standalone tools (set ``AIUPRED_DIR`` and ``MOBIDB_BINDIR``).

    Parameters
    ----------
    data_path:
        Path to ``membrane_protein_analysis_with_reduced_cc.csv``.
    disprot_path:
        Path to the DisProt TSV file.
    run_predictions:
        When ``True``, run the AIUPred and MobiDB-lite local tools for each
        protein.  Set to ``False`` to skip predictions (DisProt annotation only).
    tmd_count, max_n_term, max_cterm:
        Override the module-level TA filtering constants.

    Returns
    -------
    pd.DataFrame
        TA proteins with DisProt annotations and (if requested) disorder
        predictions appended as new columns.
    """
    print("Loading data …")
    df = pd.read_csv(data_path)
    disprot_df = pd.read_csv(disprot_path, sep="\t", dtype=str)

    print(f"  {len(df):,} proteins in main dataset")

    print("Filtering TA proteins …")
    ta_df = filter_ta_proteins(df, tmd_count=tmd_count,
                               max_n_term=max_n_term, max_cterm=max_cterm)
    print(f"  {len(ta_df):,} TA proteins identified "
          f"(TMD={tmd_count}, N-term≤{max_n_term}, C-ext≤{max_cterm})")

    print("Matching against human DisProt entries …")
    disprot_human = filter_disprot_human(disprot_df)
    print(f"  {len(disprot_human):,} human DisProt annotations")
    ta_df = match_disprot(ta_df, disprot_human)
    hits = ta_df["disprot_id"].notna().sum()
    print(f"  {hits} TA proteins found in DisProt")
    if hits:
        print(ta_df.loc[ta_df["disprot_id"].notna(),
                        ["Entry", "Entry.Name", "disprot_id",
                         "disprot_disorder_content"]].to_string(index=False))

    if run_predictions:
        print(f"\nRunning disorder predictions for {len(ta_df)} proteins …")
        ta_df = run_disorder_predictions(ta_df)
        print("  Predictions complete.")

    return ta_df


if __name__ == "__main__":
    result_df = main()
    output_path = os.path.join(_HERE, "ta_idr_results.csv")
    result_df.to_csv(output_path, index=False)
    print(f"\nSaved {len(result_df):,} rows → {output_path}")
    print(f"Columns: {list(result_df.columns)}")
