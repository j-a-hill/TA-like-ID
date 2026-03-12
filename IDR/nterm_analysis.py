#!/usr/bin/env python3
"""
N-Terminal Similarity Analysis for Tail-Anchored (TA) Proteins

Provides focused, composable functions for characterising the N-terminal
region of TA proteins:

  - Window extraction relative to the TMD boundary
  - Right-aligned sequence matrices for position-wise comparison
  - Sequence logos via ``logomaker`` (information-content representation)
  - All-vs-all pairwise identity heatmaps
  - CIDER-style IDR metrics: FCR, NCPR, κ (kappa), Ω (omega)
  - FASTA formatting helper
  - EBI MEME REST API client for de-novo motif discovery

None of the functions here are "monster functions" — each does one thing and
can be used independently or chained in a notebook or pipeline.

Dependencies
------------
  numpy, pandas, matplotlib  (already in IDR/requirements.txt)
  logomaker                  (add logomaker>=0.8 to IDR/requirements.txt)
  requests                   (already in IDR/requirements.txt)

Usage example
-------------
>>> import pandas as pd
>>> from IDR.nterm_analysis import (
...     parse_tmd_positions, extract_nterm_window, right_align_sequences,
...     plot_sequence_logo, pairwise_identity_matrix, plot_identity_heatmap,
...     cider_profile, build_cider_dataframe,
... )
"""

from __future__ import annotations

import logging
import re
import time
from typing import Optional

import numpy as np
import pandas as pd
import requests

logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────

# Standard 20 amino acids, ordered as localCIDER/IUPAC convention
_AA_ORDER: list[str] = list("ACDEFGHIKLMNPQRSTVWY")

# Charged residues
_POSITIVE_AA: frozenset[str] = frozenset("KR")
_NEGATIVE_AA: frozenset[str] = frozenset("DE")
_CHARGED_OR_PRO: frozenset[str] = frozenset("KRDEP")

# EBI JDispatcher base URL for the MEME web service
_EBI_MEME_URL: str = "https://www.ebi.ac.uk/Tools/services/rest/meme"

# ── TMD parsing ───────────────────────────────────────────────────────────────


def parse_tmd_positions(transmembrane_str: str) -> list[tuple[int, int]]:
    """Parse TMD start/end positions from a UniProt *Transmembrane* annotation.

    UniProt encodes one or more TMDs per cell in the format::

        TRANSMEM 426..446; /note="Helical"; ...
        TRANSMEM 3..23; /note="Helical; Name=1"; ...

    Parameters
    ----------
    transmembrane_str:
        Raw string from the ``Transmembrane`` column of the membrane-protein
        dataset.

    Returns
    -------
    list of (start, end) tuples (1-based, inclusive), one per TRANSMEM entry.
    Empty list if the string is missing/unparseable.
    """
    if not isinstance(transmembrane_str, str):
        return []
    matches = re.findall(r"TRANSMEM\s+(\d+)\.\.(\d+)", transmembrane_str)
    return [(int(s), int(e)) for s, e in matches]


def _last_tmd_start(transmembrane_str: str) -> Optional[int]:
    """Return the 1-based start of the last (C-terminal) TMD, or None."""
    positions = parse_tmd_positions(transmembrane_str)
    if not positions:
        return None
    return max(positions, key=lambda p: p[0])[0]


# ── Window extraction ─────────────────────────────────────────────────────────


def extract_nterm_window(
    sequence: str,
    tmd_start: int,
    window_size: int = 30,
) -> str:
    """Return the *window_size* residues immediately upstream of the TMD.

    The window is anchored at the TMD N-terminal boundary (``tmd_start``,
    1-based).  If the N-terminal domain is shorter than *window_size*, the
    entire available N-terminal sequence is returned.

    Parameters
    ----------
    sequence:
        Full protein sequence string.
    tmd_start:
        1-based position of the first TMD residue.
    window_size:
        Maximum number of N-terminal residues to include (default 30).

    Returns
    -------
    str
        Sequence fragment, at most *window_size* characters long.
    """
    if not sequence or tmd_start < 1:
        return ""
    nterm_end = tmd_start - 1           # index of last N-term residue (1-based)
    start_idx = max(0, nterm_end - window_size)
    return sequence[start_idx:nterm_end].upper()


def extract_windows_from_df(
    df: pd.DataFrame,
    sequences: dict[str, str],
    window_size: int = 30,
) -> dict[str, str]:
    """Extract N-terminal windows for each protein in *df*.

    Parameters
    ----------
    df:
        DataFrame containing at least ``Entry`` and ``Transmembrane`` columns.
    sequences:
        Mapping of UniProt accession → full protein sequence.
    window_size:
        Passed through to :func:`extract_nterm_window`.

    Returns
    -------
    dict mapping UniProt accession → window sequence.  Proteins with no
    parseable TMD or no sequence entry are omitted.
    """
    windows: dict[str, str] = {}
    for _, row in df.iterrows():
        acc = row["Entry"]
        seq = sequences.get(acc)
        if not seq:
            continue
        tmd_start = _last_tmd_start(row.get("Transmembrane", ""))
        if tmd_start is None:
            continue
        win = extract_nterm_window(seq, tmd_start, window_size=window_size)
        if win:
            windows[acc] = win
    return windows


# ── Alignment helpers ─────────────────────────────────────────────────────────


def right_align_sequences(sequences: list[str], width: Optional[int] = None) -> list[str]:
    """Right-align *sequences* by padding with ``'-'`` on the left.

    Useful for comparing N-terminal windows relative to the TMD boundary:
    the rightmost position of every aligned sequence corresponds to the
    residue immediately before the TMD in each protein.

    Parameters
    ----------
    sequences:
        List of (variable-length) sequence strings.
    width:
        Total alignment width.  Defaults to ``max(len(s) for s in sequences)``.

    Returns
    -------
    list of str, each of length *width*.
    """
    if not sequences:
        return []
    target = width if width is not None else max(len(s) for s in sequences)
    return [s.rjust(target, "-") for s in sequences]


# ── Sequence logo ─────────────────────────────────────────────────────────────


def build_logo_matrix(aligned_seqs: list[str]) -> pd.DataFrame:
    """Build a per-position amino-acid **counts** matrix from aligned sequences.

    Gap characters (``'-'``) are ignored when counting.

    Parameters
    ----------
    aligned_seqs:
        List of equal-length strings (output of :func:`right_align_sequences`).

    Returns
    -------
    pd.DataFrame with shape ``(n_positions, 20)``, columns are single-letter
    amino-acid codes, index is 0-based column position.
    """
    if not aligned_seqs:
        raise ValueError("aligned_seqs must not be empty")
    width = len(aligned_seqs[0])
    counts = np.zeros((width, len(_AA_ORDER)), dtype=float)
    aa_idx = {aa: i for i, aa in enumerate(_AA_ORDER)}

    for seq in aligned_seqs:
        if len(seq) != width:
            raise ValueError(
                f"All aligned sequences must be the same length; "
                f"expected {width}, got {len(seq)}"
            )
        for col, aa in enumerate(seq):
            if aa in aa_idx:
                counts[col, aa_idx[aa]] += 1

    return pd.DataFrame(counts, columns=_AA_ORDER)


def plot_sequence_logo(
    aligned_seqs: list[str],
    ax=None,
    title: Optional[str] = None,
) -> "matplotlib.axes.Axes":  # type: ignore[name-defined]
    """Plot an information-content sequence logo.

    Requires ``logomaker`` (``pip install logomaker``).  The logo is drawn on
    *ax* (created if not provided) and returned.

    Parameters
    ----------
    aligned_seqs:
        Equal-length strings (e.g. from :func:`right_align_sequences`).  Gap
        characters are treated as missing data.
    ax:
        Pre-existing :class:`matplotlib.axes.Axes` to draw on.
    title:
        Optional axes title.

    Returns
    -------
    matplotlib.axes.Axes
    """
    try:
        import logomaker  # noqa: PLC0415
    except ImportError as exc:
        raise ImportError("logomaker is required: pip install logomaker") from exc

    import matplotlib.pyplot as plt  # noqa: PLC0415

    counts_df = build_logo_matrix(aligned_seqs)
    info_df = logomaker.transform_matrix(counts_df, from_type="counts", to_type="information")

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6, len(aligned_seqs[0]) * 0.5), 2.5))

    logomaker.Logo(info_df, ax=ax, color_scheme="chemistry")
    ax.set_xlabel("Position relative to TMD (right-aligned)")
    ax.set_ylabel("Information (bits)")
    if title:
        ax.set_title(title)
    return ax


# ── Pairwise identity ─────────────────────────────────────────────────────────


def pairwise_identity_matrix(sequences: list[str]) -> np.ndarray:
    """Compute an all-vs-all pairwise percent-identity matrix.

    Sequences are first right-aligned to the same width (padded with ``'-'``).
    Identity is calculated over *non-gap* positions in *both* sequences.

    Parameters
    ----------
    sequences:
        List of (variable-length) amino-acid strings.

    Returns
    -------
    np.ndarray, shape ``(n, n)``, values in [0.0, 100.0].  The diagonal is 100.
    """
    aligned = right_align_sequences(sequences)
    n = len(aligned)
    mat = np.zeros((n, n), dtype=float)
    arrs = [np.frombuffer(s.encode(), dtype=np.uint8) for s in aligned]

    for i in range(n):
        mat[i, i] = 100.0
        for j in range(i + 1, n):
            a, b = arrs[i], arrs[j]
            gap = ord("-")
            # Compare positions where neither sequence has a gap
            valid = (a != gap) & (b != gap)
            n_valid = valid.sum()
            if n_valid == 0:
                pct = 0.0
            else:
                pct = float((a[valid] == b[valid]).sum()) / n_valid * 100.0
            mat[i, j] = mat[j, i] = pct

    return mat


def plot_identity_heatmap(
    pct_matrix: np.ndarray,
    labels: list[str],
    ax=None,
    title: Optional[str] = None,
) -> "matplotlib.axes.Axes":  # type: ignore[name-defined]
    """Plot a pairwise-identity heatmap.

    Parameters
    ----------
    pct_matrix:
        Square float matrix of percent identities (from
        :func:`pairwise_identity_matrix`).
    labels:
        Tick labels (e.g. protein accessions or gene names).
    ax:
        Pre-existing axes, or ``None`` to create one.
    title:
        Optional axes title.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt  # noqa: PLC0415

    n = len(labels)
    if ax is None:
        size = max(4, n * 0.25)
        _, ax = plt.subplots(figsize=(size, size))

    im = ax.imshow(pct_matrix, vmin=0, vmax=100, cmap="viridis", aspect="auto")
    plt.colorbar(im, ax=ax, label="% identity")

    if n <= 50:
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=6)
    if title:
        ax.set_title(title)
    return ax


# ── CIDER-style IDR metrics ───────────────────────────────────────────────────


def compute_fcr(sequence: str) -> float:
    """Fraction of Charged Residues (FCR).

    FCR = (K + R + D + E) / len(sequence).
    Returns 0.0 for an empty sequence.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    return sum(1 for aa in seq if aa in _POSITIVE_AA | _NEGATIVE_AA) / len(seq)


def compute_ncpr(sequence: str) -> float:
    """Net Charge Per Residue (NCPR).

    NCPR = (K + R − D − E) / len(sequence).  Positive → net cationic;
    negative → net anionic.  Returns 0.0 for an empty sequence.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    pos = sum(1 for aa in seq if aa in _POSITIVE_AA)
    neg = sum(1 for aa in seq if aa in _NEGATIVE_AA)
    return (pos - neg) / len(seq)


def compute_kappa(sequence: str, blob_size: int = 5) -> float:
    """Charge asymmetry κ (kappa).

    κ quantifies the degree to which positively and negatively charged
    residues are spatially segregated along the chain (Das & Pappu, 2013).

    κ = δ / δ_max

    where:
      • δ = mean[(σ_j − σ̄)²] over all blobs j of width *blob_size*, and
        σ_j = (f⁺_j − f⁻_j) is the local charge asymmetry;
      • δ_max is the δ for a maximally-segregated reference sequence in which
        all positively charged residues precede all negatively charged ones.

    κ → 0 : charges well-mixed; κ → 1 : maximally segregated.
    Returns 0.0 if the sequence contains no charged residues.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0

    pos_count = sum(1 for aa in seq if aa in _POSITIVE_AA)
    neg_count = sum(1 for aa in seq if aa in _NEGATIVE_AA)
    if pos_count == 0 and neg_count == 0:
        return 0.0

    def _blob_delta(s: str) -> float:
        """Mean squared deviation of per-blob σ from global σ."""
        blobs = [s[i : i + blob_size] for i in range(0, len(s), blob_size)]
        sigmas = []
        for blob in blobs:
            n = len(blob)
            fp = sum(1 for aa in blob if aa in _POSITIVE_AA) / n
            fm = sum(1 for aa in blob if aa in _NEGATIVE_AA) / n
            sigmas.append(fp - fm)
        mean_s = sum(sigmas) / len(sigmas)
        return sum((s_j - mean_s) ** 2 for s_j in sigmas) / len(sigmas)

    delta = _blob_delta(seq)

    neutral_count = len(seq) - pos_count - neg_count
    ref_seq = "K" * pos_count + "D" * neg_count + "A" * neutral_count
    delta_max = _blob_delta(ref_seq)

    if delta_max == 0.0:
        return 0.0
    return min(delta / delta_max, 1.0)


def compute_omega(sequence: str) -> float:
    """Charge-and-proline clustering Ω (omega).

    Ω measures the fraction of nearest-neighbour (i, i+1) pairs between
    residues in the set C = {K, R, D, E, P} that are C–C contacts, relative
    to all pairs that involve at least one C residue.

    Ω → 1 : C residues cluster together; Ω → 0 : C residues are dispersed.
    Returns NaN if no C residues are present, and 0.0 if no C-containing
    pairs are found.
    """
    seq = sequence.upper()
    if not seq:
        return float("nan")

    is_c = [aa in _CHARGED_OR_PRO for aa in seq]
    if not any(is_c):
        return float("nan")

    cc_pairs = 0
    cu_pairs = 0
    for i in range(len(seq) - 1):
        c_i = is_c[i]
        c_j = is_c[i + 1]
        if c_i and c_j:
            cc_pairs += 1
        elif c_i or c_j:
            cu_pairs += 1

    total = cc_pairs + cu_pairs
    return cc_pairs / total if total > 0 else 0.0


def cider_profile(sequence: str) -> dict:
    """Compute all four CIDER-style metrics for *sequence*.

    Returns
    -------
    dict with keys ``fcr``, ``ncpr``, ``kappa``, ``omega``.
    """
    return {
        "fcr": compute_fcr(sequence),
        "ncpr": compute_ncpr(sequence),
        "kappa": compute_kappa(sequence),
        "omega": compute_omega(sequence),
    }


def build_cider_dataframe(
    accessions: list[str],
    sequences: dict[str, str],
) -> pd.DataFrame:
    """Compute CIDER metrics for every sequence in *sequences*.

    Parameters
    ----------
    accessions:
        Ordered list of UniProt accessions to include (preserves row order).
    sequences:
        Mapping of accession → sequence string.

    Returns
    -------
    pd.DataFrame with columns ``Entry``, ``fcr``, ``ncpr``, ``kappa``,
    ``omega``.  Rows for accessions absent from *sequences* have NaN metrics.
    """
    rows = []
    for acc in accessions:
        seq = sequences.get(acc, "")
        if seq:
            profile = cider_profile(seq)
        else:
            profile = {"fcr": float("nan"), "ncpr": float("nan"),
                       "kappa": float("nan"), "omega": float("nan")}
        rows.append({"Entry": acc, **profile})
    return pd.DataFrame(rows)


# ── FASTA helper ──────────────────────────────────────────────────────────────


def format_fasta(sequences: dict[str, str], line_width: int = 60) -> str:
    """Format a dict of ``{accession: sequence}`` as a FASTA string.

    Parameters
    ----------
    sequences:
        Mapping of accession → bare amino-acid sequence.
    line_width:
        Number of sequence characters per line (default 60).

    Returns
    -------
    str
        Multi-FASTA formatted text, ready to write to a file or POST to an API.
    """
    lines = []
    for acc, seq in sequences.items():
        lines.append(f">{acc}")
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
    return "\n".join(lines) + "\n"


# ── EBI MEME REST client ──────────────────────────────────────────────────────


def submit_meme_ebi(
    fasta_text: str,
    email: str,
    *,
    nmotifs: int = 3,
    minw: int = 6,
    maxw: int = 15,
    url: str = _EBI_MEME_URL,
) -> str:
    """Submit a MEME motif-discovery job to the EBI JDispatcher REST API.

    Parameters
    ----------
    fasta_text:
        Multi-FASTA formatted sequences (use :func:`format_fasta`).
    email:
        Submitter e-mail address (required by the EBI API).
    nmotifs:
        Number of distinct motifs to find (default 3).
    minw:
        Minimum motif width (default 6).
    maxw:
        Maximum motif width (default 15).
    url:
        EBI MEME endpoint; override for testing.

    Returns
    -------
    str
        Job ID string (e.g. ``"meme-XXXX-XXXX-XXXX-XXXX-XXXX"``).

    Raises
    ------
    requests.HTTPError
        If the submission fails.
    """
    payload = {
        "email": email,
        "sequence": fasta_text,
        "nmotifs": str(nmotifs),
        "minw": str(minw),
        "maxw": str(maxw),
    }
    resp = requests.post(url + "/run", data=payload, timeout=30)
    resp.raise_for_status()
    return resp.text.strip()


def poll_meme_ebi(
    job_id: str,
    *,
    poll_interval: float = 30.0,
    max_wait: float = 600.0,
    url: str = _EBI_MEME_URL,
) -> str:
    """Poll the EBI JDispatcher REST API until the MEME job completes.

    Parameters
    ----------
    job_id:
        Job ID returned by :func:`submit_meme_ebi`.
    poll_interval:
        Seconds between status checks (default 30).
    max_wait:
        Maximum total wait time in seconds (default 600 = 10 min).
    url:
        EBI MEME endpoint; override for testing.

    Returns
    -------
    str
        The plain-text MEME result (the ``out`` result type).

    Raises
    ------
    TimeoutError
        If the job does not finish within *max_wait* seconds.
    requests.HTTPError
        If any HTTP request fails.
    RuntimeError
        If the job finishes with an error status.
    """
    deadline = time.monotonic() + max_wait
    while True:
        status_resp = requests.get(f"{url}/status/{job_id}", timeout=15)
        status_resp.raise_for_status()
        status = status_resp.text.strip()
        logger.debug("MEME job %s status: %s", job_id, status)

        if status in ("FINISHED", "finished"):
            result_resp = requests.get(
                f"{url}/result/{job_id}/out", timeout=30
            )
            result_resp.raise_for_status()
            return result_resp.text

        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"MEME job {job_id} failed with status: {status}")

        if time.monotonic() >= deadline:
            raise TimeoutError(
                f"MEME job {job_id} did not finish within {max_wait:.0f} s"
            )
        time.sleep(poll_interval)
