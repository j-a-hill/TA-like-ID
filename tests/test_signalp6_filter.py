"""Tests for analysis/signalp_6_filter.py."""

from __future__ import annotations

import io
import os
import sys
import tempfile

import pandas as pd
import pytest

# Make analysis module importable from the tests/ directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))

from signalp_6_filter import parse_signalp6


# ── Helpers ───────────────────────────────────────────────────────────────────


def _write_signalp_file(lines: list[str], tmp_path) -> str:
    """Write *lines* to a temp file and return its path."""
    path = str(tmp_path / "signalp_output.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


_VALID_LINE = (
    "sp|P12345|GENE1_HUMAN some protein GN=MYGENE\t"
    "OTHER\t0.95\t0.05\tCS pos: 25-26"
)
_SP_LINE = (
    "sp|P67890|GENE2_HUMAN another protein GN=ANOTHERGENE\t"
    "SP\t0.02\t0.98\tCS pos: 20-21"
)
_NO_GENE_LINE = (
    "SOMEID_no_uniprot_format\t"
    "OTHER\t0.80\t0.20\t"
)
_COMMENT_LINE = "# This is a comment"
_BLANK_LINE = ""
_SHORT_LINE = "too\tshort"


# ── parse_signalp6 ────────────────────────────────────────────────────────────


class TestParseSignalp6:
    def test_returns_dataframe(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert isinstance(result, pd.DataFrame)

    def test_required_columns_present(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        for col in ("UniProtID", "GeneName", "Prediction", "OTHER_Score", "SP_Score"):
            assert col in result.columns

    def test_parses_uniprot_id(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["UniProtID"] == "P12345"

    def test_parses_gene_name(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["GeneName"] == "MYGENE"

    def test_parses_prediction(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["Prediction"] == "OTHER"

    def test_parses_sp_prediction(self, tmp_path):
        path = _write_signalp_file([_SP_LINE], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["Prediction"] == "SP"

    def test_numeric_scores_parsed(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert abs(result.iloc[0]["OTHER_Score"] - 0.95) < 1e-6
        assert abs(result.iloc[0]["SP_Score"] - 0.05) < 1e-6

    def test_comment_lines_skipped(self, tmp_path):
        path = _write_signalp_file([_COMMENT_LINE, _VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 1

    def test_blank_lines_skipped(self, tmp_path):
        path = _write_signalp_file([_BLANK_LINE, _VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 1

    def test_short_lines_skipped(self, tmp_path):
        path = _write_signalp_file([_SHORT_LINE, _VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 1

    def test_no_uniprot_match_gives_none(self, tmp_path):
        path = _write_signalp_file([_NO_GENE_LINE], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 1
        assert result.iloc[0]["UniProtID"] is None
        assert result.iloc[0]["GeneName"] is None

    def test_multiple_records(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE, _SP_LINE], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 2

    def test_empty_file_returns_empty_df(self, tmp_path):
        path = _write_signalp_file([], tmp_path)
        result = parse_signalp6(path)
        assert len(result) == 0

    def test_cs_position_captured(self, tmp_path):
        path = _write_signalp_file([_VALID_LINE], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["CS_Position"] == "CS pos: 25-26"

    def test_missing_cs_position_is_none(self, tmp_path):
        """Lines with only 4 tab-columns should give CS_Position=None."""
        line = "sp|P11111|G_HUMAN foo GN=FOO\tOTHER\t0.9\t0.1"
        path = _write_signalp_file([line], tmp_path)
        result = parse_signalp6(path)
        assert result.iloc[0]["CS_Position"] is None
