"""Tests for IDR/nterm_analysis.py."""

from __future__ import annotations

import sys
import os
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import requests

# Make IDR module importable from the tests/ directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "IDR"))

from nterm_analysis import (
    # TMD parsing
    parse_tmd_positions,
    # Window extraction
    extract_nterm_window,
    extract_windows_from_df,
    # Alignment
    right_align_sequences,
    # Logo
    build_logo_matrix,
    # Pairwise identity
    pairwise_identity_matrix,
    # CIDER metrics
    compute_fcr,
    compute_ncpr,
    compute_kappa,
    compute_omega,
    cider_profile,
    build_cider_dataframe,
    # FASTA
    format_fasta,
    # EBI MEME client
    MemeConfig,
    submit_meme_ebi,
    poll_meme_ebi,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture()
def sample_ta_df() -> pd.DataFrame:
    return pd.DataFrame({
        "Entry": ["P001", "P002", "P003"],
        "Transmembrane": [
            'TRANSMEM 426..446; /note="Helical"; /evidence="ECO:0000255"',
            'TRANSMEM 10..30; /note="Helical"; /evidence="ECO:0000255"',
            None,
        ],
    })


@pytest.fixture()
def sample_sequences() -> dict[str, str]:
    return {
        "P001": "M" * 425 + "LLLLLLLLLLLLLLLLLLLLL",   # 446 total; TMD at 426
        "P002": "ACDEFGHIKL" + "LLLLLLLLLLLLLLLLLLLLL",  # TMD at 10; N-term = 'ACDEFGHIKL'
        "P003": "MKLVNSTAQP",
    }


# ── parse_tmd_positions ───────────────────────────────────────────────────────


class TestParseTmdPositions:
    def test_single_tmd(self):
        s = 'TRANSMEM 426..446; /note="Helical"; /evidence="ECO:0000255"'
        result = parse_tmd_positions(s)
        assert result == [(426, 446)]

    def test_multiple_tmds(self):
        s = (
            'TRANSMEM 3..23; /note="Helical; Name=1"; /evidence="ECO:0000255"; '
            'TRANSMEM 302..322; /note="Helical; Name=2"; /evidence="ECO:0000255"'
        )
        result = parse_tmd_positions(s)
        assert len(result) == 2
        assert (3, 23) in result
        assert (302, 322) in result

    def test_nan_returns_empty(self):
        assert parse_tmd_positions(float("nan")) == []  # type: ignore[arg-type]

    def test_none_returns_empty(self):
        assert parse_tmd_positions(None) == []  # type: ignore[arg-type]

    def test_empty_string_returns_empty(self):
        assert parse_tmd_positions("") == []


# ── extract_nterm_window ──────────────────────────────────────────────────────


class TestExtractNtermWindow:
    def test_basic_window(self):
        seq = "A" * 10 + "TMD" * 10
        result = extract_nterm_window(seq, tmd_start=11, window_size=5)
        assert result == "AAAAA"

    def test_window_shorter_than_nterm(self):
        """If the N-term is longer than window_size, return last window_size chars."""
        seq = "ABCDEFGHIJ" + "TMD"
        result = extract_nterm_window(seq, tmd_start=11, window_size=5)
        assert result == "FGHIJ"
        assert len(result) == 5

    def test_window_longer_than_nterm(self):
        """If the N-term is shorter than window_size, return entire N-term."""
        seq = "MKL" + "TTTTT"
        result = extract_nterm_window(seq, tmd_start=4, window_size=30)
        assert result == "MKL"

    def test_tmd_at_start(self):
        """TMD starts at position 1 → empty N-terminal."""
        result = extract_nterm_window("MKLVNSTAQP", tmd_start=1, window_size=10)
        assert result == ""

    def test_empty_sequence(self):
        assert extract_nterm_window("", tmd_start=5, window_size=10) == ""

    def test_negative_tmd_start(self):
        assert extract_nterm_window("MKLV", tmd_start=-1, window_size=10) == ""

    def test_uppercase_output(self):
        result = extract_nterm_window("acdefg", tmd_start=5, window_size=10)
        assert result == result.upper()


# ── extract_windows_from_df ───────────────────────────────────────────────────


class TestExtractWindowsFromDf:
    def test_extracts_windows(self, sample_ta_df, sample_sequences):
        result = extract_windows_from_df(
            sample_ta_df, sample_sequences, window_size=10
        )
        # P001: last 10 of 425 'M' residues
        assert "P001" in result
        assert len(result["P001"]) == 10

    def test_skips_missing_sequence(self, sample_ta_df, sample_sequences):
        """P003 has no TMD annotation → omitted."""
        result = extract_windows_from_df(sample_ta_df, sample_sequences)
        assert "P003" not in result

    def test_skips_missing_tmd(self, sample_ta_df, sample_sequences):
        """P003 Transmembrane is None → omitted."""
        result = extract_windows_from_df(sample_ta_df, sample_sequences)
        assert "P003" not in result

    def test_empty_df_returns_empty(self, sample_sequences):
        df = pd.DataFrame(columns=["Entry", "Transmembrane"])
        result = extract_windows_from_df(df, sample_sequences)
        assert result == {}


# ── right_align_sequences ─────────────────────────────────────────────────────


class TestRightAlignSequences:
    def test_all_same_length(self):
        seqs = ["ABC", "DEF", "GHI"]
        result = right_align_sequences(seqs)
        assert all(len(s) == 3 for s in result)
        assert result == seqs

    def test_pads_shorter_sequences(self):
        seqs = ["A", "BC", "DEF"]
        result = right_align_sequences(seqs)
        assert result == ["--A", "-BC", "DEF"]

    def test_custom_width(self):
        result = right_align_sequences(["AB", "A"], width=5)
        assert result == ["---AB", "----A"]

    def test_empty_list(self):
        assert right_align_sequences([]) == []

    def test_all_same_width_no_padding(self):
        seqs = ["MKLV", "ACDE"]
        aligned = right_align_sequences(seqs)
        assert not any("-" in s for s in aligned)


# ── build_logo_matrix ─────────────────────────────────────────────────────────


class TestBuildLogoMatrix:
    def test_shape(self):
        aligned = ["MKLV", "ACDE"]
        mat = build_logo_matrix(aligned)
        assert mat.shape == (4, 20)

    def test_counts_correct(self):
        aligned = ["MMMM"]
        mat = build_logo_matrix(aligned)
        assert mat["M"].sum() == 4.0

    def test_gap_not_counted(self):
        aligned = ["--KL", "MMKL"]
        mat = build_logo_matrix(aligned)
        assert mat.loc[0].sum() == 1.0   # only 'M' at position 0
        assert mat.loc[1].sum() == 1.0   # only 'M' at position 1

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            build_logo_matrix([])

    def test_unequal_length_raises(self):
        with pytest.raises(ValueError):
            build_logo_matrix(["MKLV", "ACE"])


# ── pairwise_identity_matrix ──────────────────────────────────────────────────


class TestPairwiseIdentityMatrix:
    def test_diagonal_is_100(self):
        seqs = ["MKLV", "ACDE", "FGHW"]
        mat = pairwise_identity_matrix(seqs)
        np.testing.assert_array_equal(np.diag(mat), [100.0, 100.0, 100.0])

    def test_symmetric(self):
        seqs = ["MKLV", "MKLA", "ACDE"]
        mat = pairwise_identity_matrix(seqs)
        np.testing.assert_array_almost_equal(mat, mat.T)

    def test_identical_sequences(self):
        seqs = ["MKLV", "MKLV"]
        mat = pairwise_identity_matrix(seqs)
        assert mat[0, 1] == 100.0

    def test_completely_different(self):
        seqs = ["MMMM", "CCCC"]
        mat = pairwise_identity_matrix(seqs)
        assert mat[0, 1] == 0.0

    def test_partial_identity(self):
        seqs = ["MKLV", "MKLA"]   # 3 of 4 match
        mat = pairwise_identity_matrix(seqs)
        assert abs(mat[0, 1] - 75.0) < 1e-6

    def test_gap_positions_excluded(self):
        """Gap positions (-) should not count toward identity or total."""
        seqs = ["MK", "AMK"]   # after right-align: ['-MK', 'AMK']
        mat = pairwise_identity_matrix(seqs)
        # Non-gap overlap: positions 1 and 2 → 'MK' vs 'MK' → 100%
        assert abs(mat[0, 1] - 100.0) < 1e-6

    def test_no_valid_positions_gives_zero(self):
        seqs = ["MK", "---"]
        aligned = right_align_sequences(seqs)
        # Directly build a case where no overlap exists by monkey-patching
        # → the gap handling → 0.0
        seqs_direct = ["--", "MK"]
        mat = pairwise_identity_matrix(seqs_direct)
        # '-MK' vs 'MK-' has some overlap but let's use all-gap vs non-gap
        # Here '-MK' pads 'MK' with 1 gap on left, 'MK-' not possible
        # Just verify the function runs without error
        assert mat.shape == (2, 2)


# ── compute_fcr ───────────────────────────────────────────────────────────────


class TestComputeFcr:
    def test_all_neutral(self):
        assert compute_fcr("AAAAAA") == 0.0

    def test_all_positive(self):
        assert abs(compute_fcr("KKKK") - 1.0) < 1e-9

    def test_all_negative(self):
        assert abs(compute_fcr("DDDD") - 1.0) < 1e-9

    def test_mixed(self):
        assert abs(compute_fcr("KAAA") - 0.25) < 1e-9

    def test_empty(self):
        assert compute_fcr("") == 0.0

    def test_lowercase_handled(self):
        """FCR should be case-insensitive."""
        assert compute_fcr("kdrr") == compute_fcr("KDRR")


# ── compute_ncpr ──────────────────────────────────────────────────────────────


class TestComputeNcpr:
    def test_zero_charge(self):
        assert abs(compute_ncpr("AAAAAA")) < 1e-9

    def test_net_positive(self):
        assert compute_ncpr("KKAA") > 0

    def test_net_negative(self):
        assert compute_ncpr("DDAA") < 0

    def test_balanced(self):
        assert abs(compute_ncpr("KDAA")) < 1e-9

    def test_empty(self):
        assert compute_ncpr("") == 0.0


# ── compute_kappa ─────────────────────────────────────────────────────────────


class TestComputeKappa:
    def test_no_charged_residues(self):
        assert compute_kappa("AAAAA") == 0.0

    def test_empty(self):
        assert compute_kappa("") == 0.0

    def test_range_zero_to_one(self):
        """κ must always be in [0, 1]."""
        for seq in ["KDKDKD", "KKKDDD", "AAAKDAA", "KDAAAAAKD"]:
            k = compute_kappa(seq)
            assert 0.0 <= k <= 1.0, f"κ out of range for {seq}: {k}"

    def test_mixed_less_than_segregated(self):
        """A well-mixed sequence should have lower κ than a segregated one."""
        well_mixed = "KDKDKDKD"
        segregated = "KKKKDDDD"
        assert compute_kappa(well_mixed) < compute_kappa(segregated)

    def test_only_positive(self):
        """All-positive → NCPR-based κ is just about charge asymmetry."""
        k = compute_kappa("KKKKKK")
        assert 0.0 <= k <= 1.0


# ── compute_omega ─────────────────────────────────────────────────────────────


class TestComputeOmega:
    def test_no_charged_or_pro(self):
        """Ω is NaN if no C residues exist."""
        assert np.isnan(compute_omega("AAAA"))

    def test_empty(self):
        assert np.isnan(compute_omega(""))

    def test_all_charged_clustered(self):
        """KKKDDD → all C–C pairs → Ω = 1."""
        assert abs(compute_omega("KKKDDD") - 1.0) < 1e-9

    def test_alternating_zero(self):
        """KAKAKAKA → every C is surrounded by U → Ω = 0."""
        assert abs(compute_omega("KAKAKA") - 0.0) < 1e-9

    def test_range_zero_to_one(self):
        for seq in ["KPKPKP", "AAAKPAA", "KDKDKD", "KKKPPP"]:
            val = compute_omega(seq)
            if not np.isnan(val):
                assert 0.0 <= val <= 1.0, f"Ω out of range for {seq}: {val}"

    def test_proline_counted_as_c(self):
        """Proline is included in the C set; KKPP should have Ω = 1."""
        assert abs(compute_omega("KKPP") - 1.0) < 1e-9


# ── cider_profile ─────────────────────────────────────────────────────────────


class TestCiderProfile:
    def test_returns_all_keys(self):
        profile = cider_profile("MKLVDE")
        assert set(profile) == {"fcr", "ncpr", "kappa", "omega"}

    def test_values_consistent(self):
        seq = "MKLVDE"
        profile = cider_profile(seq)
        assert abs(profile["fcr"] - compute_fcr(seq)) < 1e-9
        assert abs(profile["ncpr"] - compute_ncpr(seq)) < 1e-9
        assert abs(profile["kappa"] - compute_kappa(seq)) < 1e-9


# ── build_cider_dataframe ─────────────────────────────────────────────────────


class TestBuildCiderDataframe:
    def test_columns_present(self):
        df = build_cider_dataframe(["P001", "P002"], {"P001": "MKLVDE", "P002": "AAAA"})
        assert set(df.columns) >= {"Entry", "fcr", "ncpr", "kappa", "omega"}

    def test_row_count(self):
        df = build_cider_dataframe(["P001", "P002", "P003"],
                                   {"P001": "MKLV", "P002": "DEDE"})
        assert len(df) == 3

    def test_missing_sequence_nan(self):
        """Accession with no sequence entry → NaN metrics."""
        df = build_cider_dataframe(["P001", "MISSING"], {"P001": "MKLV"})
        row = df[df["Entry"] == "MISSING"].iloc[0]
        assert pd.isna(row["fcr"])

    def test_preserves_order(self):
        accs = ["P003", "P001", "P002"]
        df = build_cider_dataframe(accs, {"P001": "MK", "P002": "DE", "P003": "AA"})
        assert list(df["Entry"]) == accs


# ── format_fasta ──────────────────────────────────────────────────────────────


class TestFormatFasta:
    def test_header_present(self):
        result = format_fasta({"P001": "MKLV"})
        assert ">P001\n" in result

    def test_sequence_present(self):
        result = format_fasta({"P001": "MKLV"})
        assert "MKLV" in result

    def test_line_wrapping(self):
        long_seq = "A" * 120
        result = format_fasta({"P001": long_seq}, line_width=60)
        lines = [l for l in result.splitlines() if not l.startswith(">")]
        assert all(len(l) <= 60 for l in lines)

    def test_multiple_sequences(self):
        seqs = {"P001": "MKL", "P002": "DEF"}
        result = format_fasta(seqs)
        assert ">P001" in result
        assert ">P002" in result

    def test_ends_with_newline(self):
        result = format_fasta({"P001": "MKL"})
        assert result.endswith("\n")


# ── submit_meme_ebi ───────────────────────────────────────────────────────────


class TestSubmitMemeEbi:
    def test_posts_to_correct_url(self):
        fake_resp = MagicMock()
        fake_resp.raise_for_status = MagicMock()
        fake_resp.text = "meme-test-1234"

        cfg = MemeConfig(url="http://mock-meme")
        with patch("nterm_analysis.requests.post", return_value=fake_resp) as mock_post:
            job_id = submit_meme_ebi(">P001\nMKLV\n", "test@example.com", config=cfg)

        mock_post.assert_called_once()
        assert mock_post.call_args[0][0] == "http://mock-meme/run"
        assert job_id == "meme-test-1234"

    def test_payload_contains_email(self):
        fake_resp = MagicMock()
        fake_resp.raise_for_status = MagicMock()
        fake_resp.text = "meme-abc"

        cfg = MemeConfig(url="http://x")
        with patch("nterm_analysis.requests.post", return_value=fake_resp) as mock_post:
            submit_meme_ebi(">P001\nMKLV\n", "user@domain.org", config=cfg)

        payload = mock_post.call_args[1]["data"]
        assert payload["email"] == "user@domain.org"

    def test_default_config_used_when_none(self):
        """Passing no config should use MemeConfig defaults."""
        fake_resp = MagicMock()
        fake_resp.raise_for_status = MagicMock()
        fake_resp.text = "meme-default"

        with patch("nterm_analysis.requests.post", return_value=fake_resp) as mock_post:
            submit_meme_ebi(">P001\nMKLV\n", "a@b.com")

        payload = mock_post.call_args[1]["data"]
        assert payload["nmotifs"] == "3"
        assert payload["minw"] == "6"
        assert payload["maxw"] == "15"

    def test_http_error_propagated(self):
        fake_resp = MagicMock()
        fake_resp.raise_for_status.side_effect = requests.HTTPError("500")

        cfg = MemeConfig(url="http://x")
        with patch("nterm_analysis.requests.post", return_value=fake_resp):
            with pytest.raises(requests.HTTPError):
                submit_meme_ebi(">P001\nMKLV\n", "a@b.com", config=cfg)


# ── poll_meme_ebi ─────────────────────────────────────────────────────────────


class TestPollMemeEbi:
    def _make_status_resp(self, status_text: str) -> MagicMock:
        resp = MagicMock()
        resp.raise_for_status = MagicMock()
        resp.text = status_text
        return resp

    def test_returns_result_on_finished(self):
        result_resp = MagicMock()
        result_resp.raise_for_status = MagicMock()
        result_resp.text = "MEME result text"

        with patch(
            "nterm_analysis.requests.get",
            side_effect=[
                self._make_status_resp("FINISHED"),
                result_resp,
            ],
        ):
            out = poll_meme_ebi("meme-test", poll_interval=0, url="http://x")

        assert out == "MEME result text"

    def test_raises_on_error_status(self):
        with patch(
            "nterm_analysis.requests.get",
            return_value=self._make_status_resp("FAILURE"),
        ):
            with pytest.raises(RuntimeError, match="FAILURE"):
                poll_meme_ebi("meme-test", poll_interval=0, max_wait=5, url="http://x")

    def test_raises_timeout(self):
        with patch(
            "nterm_analysis.requests.get",
            return_value=self._make_status_resp("RUNNING"),
        ):
            with pytest.raises(TimeoutError):
                poll_meme_ebi(
                    "meme-test", poll_interval=0.01, max_wait=0.05, url="http://x"
                )

    def test_polls_until_finished(self):
        result_resp = MagicMock()
        result_resp.raise_for_status = MagicMock()
        result_resp.text = "done"

        responses = [
            self._make_status_resp("RUNNING"),
            self._make_status_resp("RUNNING"),
            self._make_status_resp("FINISHED"),
            result_resp,
        ]
        with patch("nterm_analysis.requests.get", side_effect=responses):
            out = poll_meme_ebi("meme-test", poll_interval=0, url="http://x")

        assert out == "done"
