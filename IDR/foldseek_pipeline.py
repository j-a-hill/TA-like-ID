#!/usr/bin/env python3
"""
Foldseek structural similarity pipeline for TA protein N-terminal (preTMD) regions.

Workflow
--------
1. ``setup_foldseek``                – locate or download foldseek binary
2. ``download_alphafold_structures`` – fetch AlphaFold PDB files via REST API
3. ``trim_structures_to_nterm``      – crop each PDB to residues before the TMD
4. ``run_foldseek_allvall``          – all-vs-all structural comparison
5. ``parse_foldseek_results``        – parse tab-separated output into DataFrame
6. ``build_tmscore_matrix``          – pivot to a symmetric TM-score matrix
7. ``plot_tmscore_heatmap``          – annotated heatmap ordered by clustering
8. ``plot_tmscore_dendrogram``       – UPGMA dendrogram of structural distances
9. ``run_pipeline``                  – convenience wrapper for the full pipeline

AlphaFold structures are fetched once and cached locally.  Foldseek is
downloaded automatically if it is not already on PATH.

Dependencies
------------
    requests, numpy, pandas, matplotlib  (already in IDR/requirements.txt)
    scipy  (for linkage / dendrogram)
    foldseek binary  (installed via setup_foldseek or pre-installed)

Usage example
-------------
>>> import pandas as pd
>>> from IDR.foldseek_pipeline import run_pipeline, plot_tmscore_heatmap
>>> ta_df = pd.read_csv("IDR/ta_idr_results.csv")
>>> seq_df = pd.read_csv("IDR/ta_sequences.csv")
>>> sequences = dict(zip(seq_df["Entry"], seq_df["sequence"]))
>>> results = run_pipeline(ta_df, sequences, work_dir="IDR/foldseek_work")
>>> plot_tmscore_heatmap(results["tmscore_matrix"], results["labels"])
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

# AlphaFold REST API endpoint — returns JSON with pdbUrl for the latest version
_AF_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# Foldseek GitHub release binary URLs (try AVX2 first, fall back to SSE2)
_FOLDSEEK_RELEASES_BASE = (
    "https://github.com/steineggerlab/foldseek/releases/latest/download"
)
_FOLDSEEK_AVX2 = f"{_FOLDSEEK_RELEASES_BASE}/foldseek-linux-avx2.tar.gz"
_FOLDSEEK_SSE2 = f"{_FOLDSEEK_RELEASES_BASE}/foldseek-linux-sse2.tar.gz"

_DEFAULT_BIN_DIR = Path.home() / ".local" / "bin"
_DEFAULT_RETRY_DELAY = 0.3   # seconds between AlphaFold API calls

# Output columns for foldseek convertalis
_FOLDSEEK_OUTPUT_FIELDS = (
    "query,target,alntmscore,lddt,rmsd,prob,evalue,bits,alnlen,qlen,tlen"
)


# ── Foldseek binary management ─────────────────────────────────────────────────


def find_foldseek() -> Optional[str]:
    """Return the path to the foldseek binary if it is on PATH, else ``None``."""
    return shutil.which("foldseek")


def install_foldseek(bin_dir: Optional[Path] = None) -> str:
    """Download and install the foldseek binary to *bin_dir*.

    Tries the AVX2 build first; falls back to SSE2 if AVX2 is not supported.

    Parameters
    ----------
    bin_dir:
        Directory to install the binary.  Defaults to ``~/.local/bin``.

    Returns
    -------
    str
        Absolute path to the installed ``foldseek`` executable.

    Raises
    ------
    RuntimeError
        If the binary cannot be downloaded or does not run.
    """
    bin_dir = Path(bin_dir) if bin_dir else _DEFAULT_BIN_DIR
    bin_dir.mkdir(parents=True, exist_ok=True)
    dest = bin_dir / "foldseek"

    for url in [_FOLDSEEK_AVX2, _FOLDSEEK_SSE2]:
        logger.info("Downloading foldseek from %s …", url)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tgz = Path(tmp) / "foldseek.tar.gz"
                resp = requests.get(url, stream=True, timeout=300)
                resp.raise_for_status()
                with open(tgz, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=65536):
                        fh.write(chunk)
                with tarfile.open(tgz) as tf:
                    # Binary lives at foldseek/bin/foldseek; match on path and
                    # require it to be a real file (not the root directory entry).
                    member = next(
                        m
                        for m in tf.getmembers()
                        if m.isfile() and (
                            m.name.endswith("/bin/foldseek")
                            or m.name == "foldseek"
                        )
                    )
                    # Use extractfile to avoid path-collision with the
                    # foldseek/ directory that also lives in the archive.
                    fobj = tf.extractfile(member)
                    if fobj is None:
                        raise RuntimeError("Could not read binary from archive")
                    with open(dest, "wb") as out:
                        shutil.copyfileobj(fobj, out)
            dest.chmod(0o755)
            # Quick smoke-test
            result = subprocess.run(
                [str(dest), "version"],
                capture_output=True, text=True, timeout=15,
            )
            if result.returncode == 0:
                logger.info("Installed foldseek → %s", dest)
                return str(dest)
            logger.warning("foldseek binary ran but returned non-zero: %s", result.stderr)
        except StopIteration:
            logger.warning("Binary not found inside archive from %s", url)
        except Exception as exc:
            logger.warning("Failed to install from %s: %s", url, exc)

    raise RuntimeError(
        "Could not install foldseek automatically.\n"
        "Please install manually: https://github.com/steineggerlab/foldseek"
    )


def setup_foldseek(bin_dir: Optional[Path] = None) -> str:
    """Return path to foldseek, downloading it if it is not already on PATH.

    Parameters
    ----------
    bin_dir:
        Where to install if not already on PATH.

    Returns
    -------
    str
        Absolute path to the foldseek executable.
    """
    path = find_foldseek()
    if path:
        logger.info("Found foldseek at %s", path)
        return path
    logger.info("foldseek not found on PATH — downloading …")
    return install_foldseek(bin_dir)


# ── AlphaFold structure download ───────────────────────────────────────────────


def download_alphafold_structure(
    accession: str,
    out_dir: Path,
    retry_delay: float = _DEFAULT_RETRY_DELAY,
) -> Optional[Path]:
    """Download the latest AlphaFold PDB structure for *accession* to *out_dir*.

    Queries the AlphaFold REST API to obtain the ``pdbUrl`` for the current
    model version (v4, v5, v6, …), so this stays correct as AlphaFold releases
    new model versions.

    Parameters
    ----------
    accession:
        UniProt accession (e.g. ``"P61981"``).
    out_dir:
        Directory to write the ``.pdb`` file.
    retry_delay:
        Seconds to sleep before each request (polite rate-limiting).

    Returns
    -------
    Path to the downloaded file, or ``None`` if download failed (e.g. 404).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{accession}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest   # already cached

    try:
        time.sleep(retry_delay)
        # Step 1: look up the current pdbUrl via the AlphaFold API
        api_resp = requests.get(_AF_API_URL.format(acc=accession), timeout=30)
        if api_resp.status_code == 404:
            logger.debug("No AlphaFold entry for %s (HTTP 404)", accession)
            return None
        api_resp.raise_for_status()
        entries = api_resp.json()
        if not entries:
            logger.debug("Empty AlphaFold API response for %s", accession)
            return None
        pdb_url = entries[0].get("pdbUrl")
        if not pdb_url:
            logger.warning("No pdbUrl in AlphaFold API response for %s", accession)
            return None

        # Step 2: download the PDB file
        time.sleep(retry_delay)
        resp = requests.get(pdb_url, timeout=60)
        if resp.status_code == 404:
            logger.debug("No AlphaFold structure for %s (HTTP 404 on pdbUrl)", accession)
            return None
        resp.raise_for_status()
        dest.write_bytes(resp.content)
        return dest
    except requests.RequestException as exc:
        logger.warning("Download failed for %s: %s", accession, exc)
        return None


def download_alphafold_structures(
    accessions: list[str],
    out_dir: Path,
    retry_delay: float = _DEFAULT_RETRY_DELAY,
) -> dict[str, Path]:
    """Download AlphaFold PDB structures for all accessions.

    Parameters
    ----------
    accessions:
        List of UniProt accessions.
    out_dir:
        Directory to write PDB files.
    retry_delay:
        Seconds between requests.

    Returns
    -------
    dict mapping accession → Path for successful downloads.
    """
    out_dir = Path(out_dir)
    results: dict[str, Path] = {}
    n = len(accessions)
    for i, acc in enumerate(accessions, 1):
        path = download_alphafold_structure(acc, out_dir, retry_delay)
        if path:
            results[acc] = path
        if i % 25 == 0 or i == n:
            logger.info("Structures: %d/%d downloaded (%d available)", i, n, len(results))
    return results


# ── PDB trimming ───────────────────────────────────────────────────────────────


def _pdb_residue_number(line: str) -> int:
    """Parse residue sequence number from columns 23–26 of a PDB ATOM line."""
    return int(line[22:26].strip())


def trim_pdb_to_nterm(
    pdb_path: Path,
    tmd_start: int,
    out_path: Path,
) -> bool:
    """Write a trimmed PDB containing only residues **before** the TMD.

    Keeps ``ATOM`` / ``HETATM`` records whose residue sequence number is
    strictly less than *tmd_start*.  Header records (MODEL, REMARK, TITLE,
    HEADER) are preserved; an ``END`` record is appended.

    Parameters
    ----------
    pdb_path:
        Input full-length PDB file.
    tmd_start:
        1-based position of the first TMD residue.
    out_path:
        Output path for the trimmed structure.

    Returns
    -------
    bool
        ``True`` if at least one residue was kept; ``False`` if the preTMD
        region is empty (e.g. TMD starts at residue 1).
    """
    pdb_path = Path(pdb_path)
    out_path = Path(out_path)
    kept: list[str] = []
    n_residues: set[int] = set()

    with open(pdb_path) as fh:
        for line in fh:
            record = line[:6].strip()
            if record in ("ATOM", "HETATM"):
                try:
                    resnum = _pdb_residue_number(line)
                except ValueError:
                    continue
                if resnum < tmd_start:
                    kept.append(line)
                    n_residues.add(resnum)
            elif record in ("MODEL", "ENDMDL", "HEADER", "TITLE", "REMARK"):
                kept.append(line)

    if not n_residues:
        return False

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return True


def trim_structures_to_nterm(
    structure_paths: dict[str, Path],
    tmd_starts: dict[str, int],
    out_dir: Path,
) -> dict[str, Path]:
    """Trim a collection of PDB structures to the preTMD region.

    Parameters
    ----------
    structure_paths:
        Mapping accession → full PDB path.
    tmd_starts:
        Mapping accession → 1-based TMD start position.
    out_dir:
        Directory to write trimmed PDB files.

    Returns
    -------
    dict mapping accession → trimmed PDB path.  Proteins with no TMD position
    or an empty preTMD region are omitted.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    trimmed: dict[str, Path] = {}

    for acc, pdb_path in structure_paths.items():
        tmd_start = tmd_starts.get(acc)
        if tmd_start is None or tmd_start <= 1:
            logger.debug("Skipping %s: no TMD start, or TMD at position 1", acc)
            continue
        out_path = out_dir / f"{acc}_nterm.pdb"
        if out_path.exists() and out_path.stat().st_size > 0:
            trimmed[acc] = out_path
            continue
        ok = trim_pdb_to_nterm(pdb_path, tmd_start, out_path)
        if ok:
            trimmed[acc] = out_path
        else:
            logger.warning(
                "Empty preTMD region for %s (tmd_start=%d) — skipping", acc, tmd_start
            )

    logger.info(
        "Trimmed %d / %d structures to preTMD region",
        len(trimmed), len(structure_paths),
    )
    return trimmed


# ── Foldseek all-vs-all search ─────────────────────────────────────────────────


def run_foldseek_allvall(
    structure_dir: Path,
    work_dir: Path,
    foldseek_bin: Optional[str] = None,
    threads: int = 4,
) -> Path:
    """Run all-vs-all Foldseek structural comparison on PDB files in *structure_dir*.

    Uses TM-align (``--alignment-type 1``) with exhaustive prefilter so that
    every pair of structures is scored regardless of sequence similarity.

    Parameters
    ----------
    structure_dir:
        Directory containing the preTMD PDB files.
    work_dir:
        Working directory for Foldseek databases and temporary files.
    foldseek_bin:
        Path to foldseek binary.  Auto-detected (or installed) if ``None``.
    threads:
        CPU threads for Foldseek.

    Returns
    -------
    Path
        Tab-separated results file (``results.m8`` inside *work_dir*).
    """
    structure_dir = Path(structure_dir)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if foldseek_bin is None:
        foldseek_bin = setup_foldseek()

    db_path     = work_dir / "structureDB"
    result_db   = work_dir / "resultDB"
    result_m8   = work_dir / "results.m8"
    tmp_path    = work_dir / "tmp"

    # ── Build structure database ──────────────────────────────────────────────
    logger.info("Building Foldseek database from %s …", structure_dir)
    _run_cmd([foldseek_bin, "createdb", str(structure_dir), str(db_path)])

    # ── All-vs-all search ─────────────────────────────────────────────────────
    logger.info("Running Foldseek all-vs-all search …")
    _run_cmd([
        foldseek_bin, "search",
        str(db_path), str(db_path),
        str(result_db),
        str(tmp_path),
        "--exhaustive-search", "1",   # disable prefilter, all pairs scored
        "--alignment-type",  "1",     # TM-align
        "-e",                "inf",   # no e-value cutoff
        "--tmscore-threshold", "0",   # report all TM scores (including 0)
        "-s",                "9.5",   # maximum sensitivity
        "--threads",         str(threads),
    ])

    # ── Convert to readable m8 format ─────────────────────────────────────────
    logger.info("Converting results to m8 format …")
    _run_cmd([
        foldseek_bin, "convertalis",
        str(db_path), str(db_path),
        str(result_db),
        str(result_m8),
        "--format-output", _FOLDSEEK_OUTPUT_FIELDS,
    ])

    logger.info("Foldseek complete → %s", result_m8)
    return result_m8


def _run_cmd(cmd: list[str]) -> None:
    """Run *cmd* as a subprocess, raising ``RuntimeError`` on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {result.returncode}):\n"
            f"  {' '.join(cmd)}\n"
            f"stdout: {result.stdout[-3000:]}\n"
            f"stderr: {result.stderr[-3000:]}"
        )


# ── Results parsing ────────────────────────────────────────────────────────────


def parse_foldseek_results(results_path: Path) -> pd.DataFrame:
    """Parse a Foldseek m8-format results file into a DataFrame.

    Parameters
    ----------
    results_path:
        Path to the tab-separated file written by ``foldseek convertalis``.

    Returns
    -------
    pd.DataFrame with columns from ``_FOLDSEEK_OUTPUT_FIELDS``.  The
    ``query`` and ``target`` columns have ``.pdb`` suffixes and ``_nterm``
    suffixes stripped so they match plain UniProt accessions.
    """
    cols = _FOLDSEEK_OUTPUT_FIELDS.split(",")
    df = pd.read_csv(results_path, sep="\t", header=None, names=cols)
    for col in ("query", "target"):
        df[col] = (
            df[col]
            .str.replace(r"\.pdb$", "", regex=True)
            .str.replace(r"_nterm$", "", regex=True)
        )
    return df


# ── TM-score matrix ────────────────────────────────────────────────────────────


def build_tmscore_matrix(
    foldseek_df: pd.DataFrame,
    accessions: list[str],
) -> np.ndarray:
    """Build a symmetric all-vs-all TM-score matrix from Foldseek results.

    Missing pairs are set to 0.  The diagonal (self vs self) is set to 1.
    Asymmetric scores are symmetrised by taking the maximum of the two
    directions (query→target, target→query).

    Parameters
    ----------
    foldseek_df:
        DataFrame from :func:`parse_foldseek_results`.
    accessions:
        Ordered list of accessions defining matrix rows/columns.

    Returns
    -------
    np.ndarray of shape ``(n, n)`` with values in [0, 1].
    """
    n = len(accessions)
    idx = {acc: i for i, acc in enumerate(accessions)}
    mat = np.zeros((n, n), dtype=float)
    np.fill_diagonal(mat, 1.0)

    for _, row in foldseek_df.iterrows():
        q_acc = row["query"]
        t_acc = row["target"]
        if q_acc == t_acc:
            continue
        score = float(row["alntmscore"])
        qi = idx.get(q_acc)
        ti = idx.get(t_acc)
        if qi is not None and ti is not None:
            # Symmetrise by keeping the higher of the two directional scores
            mat[qi, ti] = max(mat[qi, ti], score)
            mat[ti, qi] = max(mat[ti, qi], score)

    return mat


# ── Visualisation ──────────────────────────────────────────────────────────────


def plot_tmscore_heatmap(
    tmscore_matrix: np.ndarray,
    labels: list[str],
    ax=None,
    title: Optional[str] = None,
    cmap: str = "RdYlGn",
) -> "matplotlib.axes.Axes":  # type: ignore[name-defined]
    """Plot a TM-score heatmap, rows/columns reordered by UPGMA clustering.

    Parameters
    ----------
    tmscore_matrix:
        Symmetric (n, n) TM-score matrix (values 0–1).
    labels:
        Tick labels (accessions or gene names).
    ax:
        Pre-existing :class:`matplotlib.axes.Axes`; created if ``None``.
    title:
        Optional axes title.
    cmap:
        Colormap name (default ``"RdYlGn"``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform

    n = len(labels)
    dist = 1.0 - tmscore_matrix
    np.fill_diagonal(dist, 0.0)
    Z = linkage(squareform(dist, checks=False), method="average")
    order = leaves_list(Z)

    reordered = tmscore_matrix[np.ix_(order, order)]
    reordered_labels = [labels[i] for i in order]

    if ax is None:
        size = max(5, n * 0.22)
        _, ax = plt.subplots(figsize=(size, size))

    im = ax.imshow(reordered, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    plt.colorbar(im, ax=ax, label="TM-score", shrink=0.8)

    if n <= 80:
        ax.set_xticks(range(n))
        ax.set_xticklabels(reordered_labels, rotation=90, fontsize=6)
        ax.set_yticks(range(n))
        ax.set_yticklabels(reordered_labels, fontsize=6)
    if title:
        ax.set_title(title)
    return ax


def plot_tmscore_dendrogram(
    tmscore_matrix: np.ndarray,
    labels: list[str],
    ax=None,
    title: Optional[str] = None,
    color_threshold: float = 0.5,
) -> "matplotlib.axes.Axes":  # type: ignore[name-defined]
    """Plot a UPGMA dendrogram based on TM-score distance (1 − TM-score).

    Parameters
    ----------
    tmscore_matrix:
        Symmetric (n, n) TM-score matrix.
    labels:
        Leaf labels.
    ax:
        Pre-existing axes (landscape orientation recommended).
    title:
        Optional axes title.
    color_threshold:
        Distance threshold for cluster colouring.  Default 0.5 corresponds
        to proteins with average TM-score ≥ 0.5 within a cluster.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

    n = len(labels)
    dist = 1.0 - tmscore_matrix
    np.fill_diagonal(dist, 0.0)
    Z = linkage(squareform(dist, checks=False), method="average")

    if ax is None:
        _, ax = plt.subplots(figsize=(max(10, n * 0.22), 5))

    dendrogram(
        Z,
        labels=labels,
        ax=ax,
        color_threshold=color_threshold,
        leaf_rotation=90,
        leaf_font_size=6,
    )
    ax.axhline(color_threshold, linestyle="--", color="grey", linewidth=0.8,
               label=f"threshold = {color_threshold:.2f}")
    ax.set_ylabel("Distance (1 − TM-score)")
    ax.legend(fontsize=8)
    if title:
        ax.set_title(title)
    return ax


def tmscore_distribution(
    tmscore_matrix: np.ndarray,
    ax=None,
    title: Optional[str] = None,
) -> "matplotlib.axes.Axes":  # type: ignore[name-defined]
    """Histogram of off-diagonal TM-scores from the all-vs-all matrix.

    Parameters
    ----------
    tmscore_matrix:
        Symmetric (n, n) matrix.
    ax:
        Pre-existing axes; created if ``None``.
    title:
        Optional axes title.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    n = tmscore_matrix.shape[0]
    mask = ~np.eye(n, dtype=bool)
    scores = tmscore_matrix[mask]

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    ax.hist(scores, bins=40, edgecolor="k", color="steelblue", alpha=0.8)
    ax.axvline(np.median(scores), color="tomato", linestyle="--",
               label=f"median = {np.median(scores):.3f}")
    ax.axvline(0.5, color="grey", linestyle=":", linewidth=0.8,
               label="TM-score = 0.5")
    ax.set_xlabel("TM-score")
    ax.set_ylabel("Pair count")
    ax.legend(fontsize=9)
    if title:
        ax.set_title(title)
    return ax


# ── Top-hit summary table ──────────────────────────────────────────────────────


def top_hits(
    foldseek_df: pd.DataFrame,
    top_n: int = 5,
    min_tmscore: float = 0.5,
) -> pd.DataFrame:
    """Return the top structural hits per query protein.

    Parameters
    ----------
    foldseek_df:
        DataFrame from :func:`parse_foldseek_results`.
    top_n:
        Maximum hits per query.
    min_tmscore:
        Minimum TM-score to include.

    Returns
    -------
    pd.DataFrame with columns ``query``, ``target``, ``alntmscore``,
    ``lddt``, ``rmsd``, sorted by query then TM-score descending.
    """
    df = foldseek_df[
        (foldseek_df["query"] != foldseek_df["target"])
        & (foldseek_df["alntmscore"] >= min_tmscore)
    ].copy()
    df = df.sort_values(["query", "alntmscore"], ascending=[True, False])
    return (
        df.groupby("query", sort=False)
          .head(top_n)
          .reset_index(drop=True)
        [["query", "target", "alntmscore", "lddt", "rmsd"]]
    )


# ── Convenience pipeline wrapper ───────────────────────────────────────────────


def run_pipeline(
    ta_df: pd.DataFrame,
    sequences: dict[str, str],
    work_dir: "str | Path",
    foldseek_bin: Optional[str] = None,
    force_rerun: bool = False,
    threads: int = 4,
) -> dict:
    """End-to-end Foldseek pipeline for TA N-terminal structural similarity.

    Steps
    -----
    1. Parse TMD start positions from ``ta_df["Transmembrane"]``
    2. Download AlphaFold structures for all proteins with a TMD position
    3. Trim each structure to the preTMD region (residues < tmd_start)
    4. Run Foldseek all-vs-all on the trimmed structures
    5. Parse results and build a symmetric TM-score matrix

    Parameters
    ----------
    ta_df:
        DataFrame with columns ``Entry`` and ``Transmembrane``.
    sequences:
        Mapping accession → full protein sequence.  Used only to determine
        which proteins exist (structures are fetched from AlphaFold).
    work_dir:
        Root directory for all pipeline outputs.  Sub-directories are created
        automatically:

        * ``alphafold_structures/``  – downloaded full-length PDB files
        * ``trimmed_structures/``    – preTMD-trimmed PDB files
        * ``foldseek/``              – Foldseek databases and results

    foldseek_bin:
        Path to the foldseek binary.  Auto-detected or installed if ``None``.
    force_rerun:
        Re-run Foldseek even if ``foldseek/results.m8`` already exists.
    threads:
        CPU threads for Foldseek.

    Returns
    -------
    dict with keys:

    ``tmd_starts``       – ``dict[str, int]`` of TMD start positions
    ``structure_paths``  – ``dict[str, Path]`` of downloaded full-length PDBs
    ``trimmed_paths``    – ``dict[str, Path]`` of trimmed preTMD PDBs
    ``foldseek_df``      – raw Foldseek results DataFrame
    ``tmscore_matrix``   – symmetric TM-score matrix (numpy array)
    ``accessions``       – ordered list of accessions in the matrix
    ``labels``           – display labels (gene first word, fallback accession)
    """
    # Lazy import to avoid circular dependency when used as a sub-module
    try:
        from nterm_analysis import parse_tmd_positions  # type: ignore
    except ImportError:
        from IDR.nterm_analysis import parse_tmd_positions  # type: ignore

    work_dir = Path(work_dir)
    raw_dir      = work_dir / "alphafold_structures"
    trimmed_dir  = work_dir / "trimmed_structures"
    foldseek_dir = work_dir / "foldseek"
    results_file = foldseek_dir / "results.m8"

    # ── 1. Parse TMD positions ────────────────────────────────────────────────
    logger.info("Parsing TMD positions for %d proteins …", len(ta_df))
    tmd_starts: dict[str, int] = {}
    for _, row in ta_df.iterrows():
        acc = row["Entry"]
        positions = parse_tmd_positions(row.get("Transmembrane", ""))
        if positions:
            # Use the last (most C-terminal) TMD — the defining TA anchor
            tmd_starts[acc] = max(positions, key=lambda p: p[0])[0]

    accessions = [acc for acc in ta_df["Entry"] if acc in tmd_starts]
    logger.info("%d / %d proteins have a parseable TMD start", len(accessions), len(ta_df))

    # ── 2. Download AlphaFold structures ──────────────────────────────────────
    logger.info("Downloading AlphaFold structures to %s …", raw_dir)
    structure_paths = download_alphafold_structures(accessions, raw_dir)
    logger.info("%d / %d structures available from AlphaFold", len(structure_paths), len(accessions))

    # ── 3. Trim to N-terminal region ──────────────────────────────────────────
    logger.info("Trimming structures to preTMD region …")
    trimmed_paths = trim_structures_to_nterm(structure_paths, tmd_starts, trimmed_dir)

    if not trimmed_paths:
        raise RuntimeError(
            "No trimmed structures were produced — cannot run Foldseek."
        )

    # ── 4. Run Foldseek ───────────────────────────────────────────────────────
    if results_file.exists() and not force_rerun:
        logger.info("Using cached Foldseek results: %s", results_file)
    else:
        if foldseek_bin is None:
            foldseek_bin = setup_foldseek()
        results_file = run_foldseek_allvall(
            trimmed_dir, foldseek_dir, foldseek_bin, threads=threads,
        )

    # ── 5. Parse results ──────────────────────────────────────────────────────
    foldseek_df = parse_foldseek_results(results_file)

    # Build matrix only for proteins that have trimmed structures
    matrix_accs = [acc for acc in accessions if acc in trimmed_paths]
    tmscore_matrix = build_tmscore_matrix(foldseek_df, matrix_accs)

    # ── Labels ────────────────────────────────────────────────────────────────
    acc_to_gene = {}
    if "Gene.Names" in ta_df.columns:
        acc_to_gene = dict(zip(ta_df["Entry"], ta_df["Gene.Names"]))
    labels = [
        str(acc_to_gene.get(acc, acc)).split()[0]
        if pd.notna(acc_to_gene.get(acc))
        else acc
        for acc in matrix_accs
    ]

    logger.info(
        "Pipeline complete: %d structures in a %d×%d TM-score matrix",
        len(trimmed_paths), *tmscore_matrix.shape,
    )
    return {
        "tmd_starts":      tmd_starts,
        "structure_paths": structure_paths,
        "trimmed_paths":   trimmed_paths,
        "foldseek_df":     foldseek_df,
        "tmscore_matrix":  tmscore_matrix,
        "accessions":      matrix_accs,
        "labels":          labels,
    }
