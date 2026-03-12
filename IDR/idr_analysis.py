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

    ``aiupred_lib.py`` is imported directly (no subprocess) so the models are
    loaded from ``<AIUPRED_DIR>/data/`` once and reused for all sequences.

MobiDB-lite (https://github.com/BioComputingUP/MobiDB-lite):
    git clone https://github.com/BioComputingUP/MobiDB-lite.git /path/to/mobidb-lite
    pip install /path/to/mobidb-lite          # installs Python package + bin/ directory
    export MOBIDB_BINDIR=/path/to/mobidb-lite/src/mobidb_lite/bin

    The ``bin/`` directory contains compiled sub-predictor binaries (IUPred,
    DisEMBL, ESpritz, TISEAN, SEG).  It is NOT included when you run
    ``pip install git+https://...``; you must clone the repo and either
    install locally (``pip install .``) or set ``MOBIDB_BINDIR`` manually.
    The pipeline also tries to auto-detect ``bin/`` from the installed package.

Outputs a DataFrame with all TA-protein columns plus appended disorder annotations.
"""

import contextlib
import gc
import importlib.util
import io
import logging
import multiprocessing as _mp
import os
import queue as _queue
import time

import pandas as pd
import requests

logger = logging.getLogger(__name__)

# ── Tunable parameters ────────────────────────────────────────────────────────
# A TA protein has exactly one transmembrane domain (TMD), located towards
# the C-terminus with minimal sequence extension after it.

EXACT_TMD_COUNT: int = 1       # Total membrane domains (TMD + intramembrane)
MAX_N_TERM_DOMAINS: int = 0    # N-terminal membrane domains (0 = C-terminal only)
MAX_CTERM_EXTENSION: int = 10  # Max residues between TMD end and protein C-terminus

AIUPRED_THRESHOLD: float = 0.5   # Per-residue score ≥ this → predicted disordered
API_DELAY: float = 0.3           # Seconds to pause between UniProt sequence fetches

HUMAN_TAXON_ID: str = "9606"

# ── Local tool paths (override via environment variables) ─────────────────────
# Directory of the cloned doszilab/AIUPred repository.
AIUPRED_DIR: str = os.environ.get("AIUPRED_DIR", "")

# Directory containing the compiled MobiDB-lite binary sub-predictors.
# If not set explicitly, the pipeline tries to find bin/ inside the installed
# mobidb_lite package (works when pip-installed from a local clone of the repo).
def _detect_mobidb_bindir() -> str:
    explicit = os.environ.get("MOBIDB_BINDIR", "")
    if explicit:
        # Always resolve to an absolute path so sub-predictor subprocesses can
        # find their binaries regardless of the working directory at call time.
        return os.path.abspath(explicit)
    try:
        import mobidb_lite as _mobi  # type: ignore[import]
        candidate = os.path.join(os.path.dirname(_mobi.__file__), "bin")
        if os.path.isdir(candidate):
            return os.path.abspath(candidate)
    except ImportError:
        pass
    return ""

MOBIDB_BINDIR: str = _detect_mobidb_bindir()

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

# Module-level cache: aiupred_dir → loaded aiupred_lib module (or None on failure)
_aiupred_cache: dict = {}


def _load_aiupred_lib(aiupred_dir: str):
    """Import ``aiupred_lib`` directly from the cloned AIUPred repository.

    Uses :mod:`importlib.util` so the caller's ``sys.path`` is not polluted.
    The loaded module is cached keyed by *aiupred_dir* so that the expensive
    PyTorch model weights are only read from disk once per interpreter session.

    Returns the module, or ``None`` when the library cannot be loaded.
    """
    if aiupred_dir in _aiupred_cache:
        return _aiupred_cache[aiupred_dir]

    lib_path = os.path.join(aiupred_dir, "aiupred_lib.py")
    if not os.path.isfile(lib_path):
        logger.warning(
            "AIUPred: aiupred_lib.py not found in %s – AIUPred predictions will be skipped. "
            "Clone https://github.com/doszilab/AIUPred.git and set AIUPRED_DIR.",
            aiupred_dir,
        )
        _aiupred_cache[aiupred_dir] = None
        return None

    try:
        spec = importlib.util.spec_from_file_location("_aiupred_lib", lib_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)  # executes module code; __file__=lib_path so data/ is found
        # Pre-load the disorder models once; stored in the cached module object
        mod._disorder_models = mod.init_models("disorder", force_cpu=True)
        _aiupred_cache[aiupred_dir] = mod
        return mod
    except Exception as exc:
        logger.warning("AIUPred: failed to load aiupred_lib from %s: %s", aiupred_dir, exc)
        _aiupred_cache[aiupred_dir] = None
        return None


def run_aiupred_local(sequence: str, aiupred_dir: str = AIUPRED_DIR) -> dict:
    """Run AIUPred disorder prediction using the locally cloned standalone tool.

    Clone the repository and set ``AIUPRED_DIR`` (or pass *aiupred_dir*):

    .. code-block:: bash

        git clone https://github.com/doszilab/AIUPred.git /path/to/aiupred
        pip install torch>=2.6.0 numpy scipy
        export AIUPRED_DIR=/path/to/aiupred

    The function imports ``aiupred_lib.py`` directly via :mod:`importlib.util`
    and calls ``predict_disorder`` with pre-loaded models (loaded once per
    session and cached).  No subprocess is spawned and no text output is parsed.

    Returns
    -------
    dict
        ``{"scores": [<float>, ...]}`` with one score per residue, or ``{}``
        when *aiupred_dir* is unset, the library is not found, or prediction fails.
    """
    if not aiupred_dir or not os.path.isdir(aiupred_dir):
        if aiupred_dir:
            logger.warning("AIUPred: AIUPRED_DIR '%s' does not exist.", aiupred_dir)
        return {}

    mod = _load_aiupred_lib(aiupred_dir)
    if mod is None:
        return {}

    try:
        embedding_model, reg_model, device = mod._disorder_models
        # torch.no_grad() prevents gradient graph accumulation during inference,
        # which would otherwise grow memory linearly with protein count.
        try:
            import torch
            _no_grad: contextlib.AbstractContextManager = torch.no_grad()
        except ImportError:
            _no_grad = contextlib.nullcontext()
        with _no_grad:
            scores = mod.predict_disorder(
                sequence, embedding_model, reg_model, device, smoothing=True
            )
        return {"scores": list(scores)}
    except Exception as exc:
        logger.warning("AIUPred: prediction failed for sequence of length %d: %s", len(sequence), exc)
        return {}


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

    Install the package from a **local clone** so that the compiled
    sub-predictor binaries are available alongside the Python package:

    .. code-block:: bash

        git clone https://github.com/BioComputingUP/MobiDB-lite.git /path/to/mobidb-lite
        pip install /path/to/mobidb-lite

    The ``pip install git+https://...`` form does **not** ship the ``bin/``
    directory; install from a local clone instead, or set ``MOBIDB_BINDIR``
    explicitly to ``/path/to/mobidb-lite/src/mobidb_lite/bin``.

    The function uses ``mobidb_lite.consensus.run`` directly (no subprocess),
    passing a :class:`io.StringIO` FASTA handle and the *bindir* path.

    Returns
    -------
    dict
        ``{"prediction-disorder-mobidb_lite": {"regions": [[start, end], ...]}}``
        compatible with :func:`mobidb_lite_summary`, or ``{}`` when *bindir* is
        unset, the package is not installed, or the prediction fails.
    """
    if not bindir or not os.path.isdir(bindir):
        logger.warning(
            "MobiDB-lite: bin/ directory not found (bindir=%r). "
            "Clone https://github.com/BioComputingUP/MobiDB-lite.git, "
            "run 'pip install .' from the clone, and set MOBIDB_BINDIR if needed.",
            bindir,
        )
        return {}

    try:
        from mobidb_lite.consensus import run as _mobi_run  # type: ignore[import]
    except ImportError:
        logger.warning("MobiDB-lite: mobidb_lite package is not installed.")
        return {}

    fasta_handle = io.StringIO(f">query\n{sequence}\n")
    try:
        for _seq_id, regions, _scores in _mobi_run(fasta_handle, bindir, "", threads=1):
            if regions is None:
                logger.warning(
                    "MobiDB-lite: prediction failed for sequence of length %d "
                    "(some sub-predictors may be missing from %s).",
                    len(sequence),
                    bindir,
                )
                return {}
            mobi_regions = regions.get("mobidblite", [])
            return {"prediction-disorder-mobidb_lite": {"regions": mobi_regions}}
    except Exception as exc:
        if "TimeoutExpired" in type(exc).__name__:
            logger.warning(
                "MobiDB-lite: sub-predictor timed out for sequence of length %d", len(sequence)
            )
        else:
            logger.warning("MobiDB-lite: prediction raised an exception: %s", exc)
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
    if not isinstance(pred, dict):
        # Unexpected shape — treat as tool failure
        return {"mobidb_lite_disordered_fraction": None}
    regions = pred.get("regions", [])

    # Empty regions list means the tool ran but found no disordered residues
    # (fully ordered) → 0.0, not None.  None is reserved for true failures
    # where run_mobidb_lite_local returned {} (bindir missing, exception, etc.).
    if not regions:
        return {"mobidb_lite_disordered_fraction": 0.0}

    disordered_residues = sum(
        int(end) - int(start) + 1 for start, end in regions
    )
    fraction = min(disordered_residues / protein_length, 1.0)
    return {"mobidb_lite_disordered_fraction": round(fraction, 4)}


# ── Full prediction loop ──────────────────────────────────────────────────────

_PRED_COLS = ("aiupred_mean_score", "aiupred_disordered_fraction",
              "mobidb_lite_disordered_fraction")

# Maximum wall-clock seconds to wait for one protein's predictions.
# If the worker process hangs longer than this it is killed and restarted.
# Each MobiDB-lite sub-predictor has a 30s timeout (8 predictors max = 240s),
# plus AIUPred time and communication overhead.  Set worker budget generously.
_PER_PROTEIN_TIMEOUT: int = 360  # 6 minutes

_NULL_PRED: dict = {
    "aiupred_mean_score": None,
    "aiupred_disordered_fraction": None,
    "mobidb_lite_disordered_fraction": None,
}


# ── Persistent worker ─────────────────────────────────────────────────────────

def _worker_fn(job_q: "_mp.Queue", result_q: "_mp.Queue",
               aiupred_dir: str, mobidb_bindir: str) -> None:
    """Long-running child process: loads AIUPred models once, then loops.

    Receives ``(acc, sequence, length)`` tuples from *job_q*, puts prediction
    dicts into *result_q*.  A sentinel ``None`` job causes the worker to exit.
    """
    mod = _load_aiupred_lib(aiupred_dir)

    while True:
        job = job_q.get()
        if job is None:
            return

        acc, sequence, length = job
        try:
            if sequence:
                # AIUPred ── no_grad prevents gradient graph accumulation
                if mod is not None:
                    try:
                        import torch
                        _no_grad: contextlib.AbstractContextManager = torch.no_grad()
                    except ImportError:
                        _no_grad = contextlib.nullcontext()
                    with _no_grad:
                        emb, reg, dev = mod._disorder_models
                        scores = mod.predict_disorder(
                            sequence, emb, reg, dev, smoothing=True
                        )
                    ai_metrics = aiupred_summary(list(scores))
                else:
                    ai_metrics = {"aiupred_mean_score": None,
                                  "aiupred_disordered_fraction": None}

                mobi_raw = run_mobidb_lite_local(sequence, mobidb_bindir)
                mobi_metrics = mobidb_lite_summary(mobi_raw, length)
            else:
                ai_metrics = {"aiupred_mean_score": None,
                              "aiupred_disordered_fraction": None}
                mobi_metrics = {"mobidb_lite_disordered_fraction": None}

            result_q.put({**ai_metrics, **mobi_metrics})

        except Exception as exc:  # noqa: BLE001
            logger.warning("Worker: unhandled error for %s: %s", acc, exc)
            result_q.put(dict(_NULL_PRED))

        finally:
            gc.collect()


class _PredictionWorker:
    """Manages a persistent child process that runs disorder predictions.

    The worker is forked at creation (inheriting already-loaded module state)
    and restarted automatically whenever it hangs or crashes.
    """

    def __init__(self, timeout: int = _PER_PROTEIN_TIMEOUT) -> None:
        self.timeout = timeout
        self._job_q: "_mp.Queue" = _mp.Queue()
        self._result_q: "_mp.Queue" = _mp.Queue()
        self._proc: "_mp.Process | None" = None
        self._start()

    def _start(self) -> None:
        self._job_q = _mp.Queue()
        self._result_q = _mp.Queue()
        self._proc = _mp.Process(
            target=_worker_fn,
            args=(self._job_q, self._result_q, AIUPRED_DIR, MOBIDB_BINDIR),
            daemon=True,
        )
        self._proc.start()

    def _stop(self) -> None:
        if self._proc is None:
            return
        self._proc.terminate()
        self._proc.join(5)
        if self._proc.is_alive():
            self._proc.kill()
            self._proc.join()
        self._proc = None

    def predict(self, acc: str, sequence: "str | None", length: int) -> dict:
        """Return prediction dict for *sequence*.

        Sends the job to the worker and waits up to ``self.timeout`` seconds.
        If the worker does not respond in time it is killed and restarted;
        ``_NULL_PRED`` is returned for the timed-out protein.
        """
        self._job_q.put((acc, sequence, length))
        t0 = time.monotonic()
        try:
            return self._result_q.get(timeout=self.timeout)
        except _queue.Empty:
            elapsed = time.monotonic() - t0
            logger.warning(
                "Worker timeout after %.0fs (budget=%ds) for %s (len=%d) – killing and restarting worker.",
                elapsed, self.timeout, acc, length,
            )
            self._stop()
            self._start()
            return dict(_NULL_PRED)

    def shutdown(self) -> None:
        """Send the sentinel and wait for the worker to exit cleanly."""
        try:
            self._job_q.put(None)
            if self._proc is not None:
                self._proc.join(10)
        finally:
            self._stop()


def run_disorder_predictions(
    ta_df: pd.DataFrame,
    checkpoint_path: str | None = None,
) -> pd.DataFrame:
    """Add AIUPred and MobiDB-lite predictions to each row of *ta_df*.

    Fetches the amino-acid sequence from UniProt for each protein, then runs
    both predictors in an isolated worker process with a hard per-protein
    timeout.  Failed or timed-out proteins produce ``None`` values.

    Checkpoints are written to a separate ``.ckpt`` sidecar file after every
    protein so an interrupted run can resume cheaply.
    """
    # Sidecar checkpoint – never shares a name with the output file
    ckpt_file = (checkpoint_path + ".ckpt") if checkpoint_path else None

    # ── resume ────────────────────────────────────────────────────────────────
    completed: dict[str, dict] = {}
    if ckpt_file and os.path.isfile(ckpt_file):
        try:
            ckpt = pd.read_csv(ckpt_file)
            for _, r in ckpt.iterrows():
                completed[r["Entry"]] = {
                    "aiupred_mean_score": r["aiupred_mean_score"] if pd.notna(r["aiupred_mean_score"]) else None,
                    "aiupred_disordered_fraction": r["aiupred_disordered_fraction"] if pd.notna(r["aiupred_disordered_fraction"]) else None,
                    "mobidb_lite_disordered_fraction": r["mobidb_lite_disordered_fraction"] if pd.notna(r["mobidb_lite_disordered_fraction"]) else None,
                }
            if completed:
                print(f"  Resuming: {len(completed)} proteins already in checkpoint, skipping them.")
        except Exception as exc:  # noqa: BLE001
            logger.warning("Checkpoint load failed (%s) – starting from scratch.", exc)
            completed = {}

    worker = _PredictionWorker(timeout=_PER_PROTEIN_TIMEOUT)
    records = []
    total = len(ta_df)
    try:
        for i, (idx, row) in enumerate(ta_df.iterrows(), start=1):
            acc = row["Entry"]
            length = int(row["Length"]) if pd.notna(row.get("Length")) else 0

            if acc in completed:
                print(f"  [{i}/{total}] {acc} (cached ✓)")
                records.append((idx, completed[acc]))
                continue

            print(f"  [{i}/{total}] {acc} (len={length})", end=" ", flush=True)

            # Sequence fetch stays in the main process (controls API rate)
            sequence = fetch_sequence(acc)

            t_pred = time.monotonic()
            pred = worker.predict(acc, sequence, length)
            elapsed_pred = time.monotonic() - t_pred

            timed_out = all(v is None for v in pred.values())
            status = f"TIMEOUT – skipped ({elapsed_pred:.0f}s)" if timed_out else f"✓ ({elapsed_pred:.0f}s)"
            print(status)

            records.append((idx, pred))

            if ckpt_file:
                ckpt_row = pd.DataFrame([{"Entry": acc, **pred}])
                write_header = not os.path.isfile(ckpt_file)
                ckpt_row.to_csv(ckpt_file, mode="a", index=False, header=write_header)

    finally:
        worker.shutdown()

    indices, dicts = zip(*records) if records else ([], [])
    pred_df = pd.DataFrame(list(dicts), index=list(indices))
    return pd.concat([ta_df, pred_df], axis=1)


# ── Main entry point ──────────────────────────────────────────────────────────

def main(
    data_path: str = DATA_PATH,
    disprot_path: str = DISPROT_PATH,
    run_predictions: bool = True,
    filter_ta: bool = True,
    tmd_count: int = EXACT_TMD_COUNT,
    max_n_term: int = MAX_N_TERM_DOMAINS,
    max_cterm: int = MAX_CTERM_EXTENSION,
    output_path: str | None = None,
    resume: bool = True,
) -> pd.DataFrame:
    """Run the full IDR analysis pipeline.

    1. Load the membrane-protein DataFrame and the DisProt TSV.
    2. Optionally filter to TA proteins (default) or process the whole dataset.
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
    filter_ta:
        When ``True`` (default), restrict to TA proteins before running
        predictions.  Set to ``False`` to process every protein in *data_path*
        — useful for building a whole-population baseline for enrichment
        comparisons.  When ``False`` a boolean ``is_ta`` column is added so
        that TA proteins can still be identified downstream.
    tmd_count, max_n_term, max_cterm:
        Override the module-level TA filtering constants.
    output_path:
        Where to write the final CSV.  Defaults to ``ta_idr_results.csv`` (TA
        mode) or ``all_idr_results.csv`` (whole-population mode) next to this
        script.  Also used as the per-protein checkpoint file so that
        interrupted runs can be resumed automatically.
    resume:
        When ``True`` (default) and *output_path* already exists, skip proteins
        that are already present in that file and append new predictions to it.
        Set to ``False`` to force a clean run from scratch.

    Returns
    -------
    pd.DataFrame
        Proteins with DisProt annotations and (if requested) disorder
        predictions appended as new columns.
    """
    if output_path is None:
        default_name = "ta_idr_results.csv" if filter_ta else "all_idr_results.csv"
        output_path = os.path.join(_HERE, default_name)

    print("Loading data …")
    df = pd.read_csv(data_path)
    disprot_df = pd.read_csv(disprot_path, sep="\t", dtype=str)

    print(f"  {len(df):,} proteins in main dataset")

    if filter_ta:
        print("Filtering TA proteins …")
        working_df = filter_ta_proteins(df, tmd_count=tmd_count,
                                        max_n_term=max_n_term, max_cterm=max_cterm)
        print(f"  {len(working_df):,} TA proteins identified "
              f"(TMD={tmd_count}, N-term≤{max_n_term}, C-ext≤{max_cterm})")
    else:
        print("Processing whole dataset (no TA filter) …")
        working_df = df.copy()
        # Label which proteins meet TA criteria so downstream analysis can filter
        for col in ("membrane_domain_count", "N_term_md", "cterm_distance"):
            working_df[col] = pd.to_numeric(working_df[col], errors="coerce")
        ta_mask = (
            (working_df["membrane_domain_count"] == tmd_count)
            & (working_df["N_term_md"] <= max_n_term)
            & (working_df["cterm_distance"] <= max_cterm)
        )
        working_df["is_ta"] = ta_mask
        ta_count = ta_mask.sum()
        print(f"  {ta_count:,} proteins labelled is_ta=True "
              f"(TMD={tmd_count}, N-term≤{max_n_term}, C-ext≤{max_cterm})")

    print("Matching against human DisProt entries …")
    disprot_human = filter_disprot_human(disprot_df)
    print(f"  {len(disprot_human):,} human DisProt annotations")
    working_df = match_disprot(working_df, disprot_human)
    hits = working_df["disprot_id"].notna().sum()
    print(f"  {hits} proteins found in DisProt")
    if hits and filter_ta:
        print(working_df.loc[working_df["disprot_id"].notna(),
                             ["Entry", "Entry.Name", "disprot_id",
                              "disprot_disorder_content"]].to_string(index=False))

    if run_predictions:
        checkpoint = output_path if resume else None
        print(f"\nRunning disorder predictions for {len(working_df):,} proteins …")
        ckpt_sidecar = (output_path + ".ckpt") if checkpoint else None
        if ckpt_sidecar and os.path.isfile(ckpt_sidecar):
            print(f"  Checkpoint file found: {ckpt_sidecar}")
        working_df = run_disorder_predictions(working_df, checkpoint_path=checkpoint)
        print("  Predictions complete.")

    return working_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="IDR analysis pipeline for tail-anchored proteins."
    )
    parser.add_argument(
        "--output", default=os.path.join(_HERE, "ta_idr_results.csv"),
        help="Path for the output CSV (default: IDR/ta_idr_results.csv).",
    )
    parser.add_argument(
        "--no-resume", action="store_true",
        help="Start a clean run, ignoring any existing output file.",
    )
    parser.add_argument(
        "--no-predictions", action="store_true",
        help="Skip AIUPred/MobiDB-lite – DisProt annotation only.",
    )
    parser.add_argument(
        "--all-proteins", action="store_true",
        help=(
            "Process every protein in the dataset instead of filtering to TA proteins. "
            "Adds an is_ta column so TA proteins can be identified downstream. "
            "Default output file: IDR/all_idr_results.csv."
        ),
    )
    args = parser.parse_args()

    # When --all-proteins is set and no explicit --output was given, use all_idr_results.csv
    output = args.output
    if args.all_proteins and output == os.path.join(_HERE, "ta_idr_results.csv"):
        output = os.path.join(_HERE, "all_idr_results.csv")

    result_df = main(
        output_path=output,
        resume=not args.no_resume,
        run_predictions=not args.no_predictions,
        filter_ta=not args.all_proteins,
    )
    result_df.to_csv(output, index=False)
    print(f"\nSaved {len(result_df):,} rows → {output}")
    print(f"Columns: {list(result_df.columns)}")
