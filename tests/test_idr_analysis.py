"""Tests for IDR/idr_analysis.py."""

import sys
import os

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

# Make IDR module importable from the tests/ directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "IDR"))

from idr_analysis import (
    filter_ta_proteins,
    filter_disprot_human,
    match_disprot,
    aiupred_summary,
    mobidb_lite_summary,
    run_disorder_predictions,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture()
def sample_protein_df() -> pd.DataFrame:
    """Small DataFrame that mirrors the main membrane-protein dataset schema."""
    return pd.DataFrame({
        "Entry":                ["P001", "P002", "P003", "P004", "P005"],
        "Entry.Name":           ["A_HUMAN", "B_HUMAN", "C_HUMAN", "D_HUMAN", "E_HUMAN"],
        "Length":               [200, 300, 150, 400, 250],
        "membrane_domain_count":[1,   1,   2,   1,   1],
        "N_term_md":            [0,   1,   0,   0,   0],
        "cterm_distance":       [10,  20,  15,  50,  25],
        "Prediction":           ["OTHER", "OTHER", "SP", "OTHER", "OTHER"],
        "Reduced.CC.Terms":     ["ER", "Golgi", "ER", "plasma membrane", "ER"],
    })


@pytest.fixture()
def sample_disprot_df() -> pd.DataFrame:
    """Small DisProt-formatted DataFrame with two entries: one human, one mouse."""
    return pd.DataFrame({
        "UniProt ACC":            ["P001", "P001", "P006"],
        "DisProt ID":             ["DP00001", "DP00001", "DP00002"],
        "Protein name":           ["Protein A", "Protein A", "Protein F"],
        "Gene name":              ["GENE_A", "GENE_A", "GENE_F"],
        "Sequence length":        ["200", "200", "350"],
        "Organism":               ["Homo sapiens", "Homo sapiens", "Mus musculus"],
        "NCBI Taxon ID":          ["9606", "9606", "10090"],
        "Protein Disorder Content": ["0.25", "0.25", "0.35"],
        "Region ID":              ["DP00001r001", "DP00001r002", "DP00002r001"],
        "Start":                  ["10", "80", "20"],
        "End":                    ["60", "120", "80"],
    })


# ── filter_ta_proteins ────────────────────────────────────────────────────────

class TestFilterTaProteins:
    def test_returns_single_tmd_c_terminal(self, sample_protein_df):
        """Proteins with 1 TMD, 0 N-term domains, and cterm_distance ≤ 30."""
        result = filter_ta_proteins(sample_protein_df)
        assert set(result["Entry"]) == {"P001", "P005"}

    def test_excludes_multi_tmd(self, sample_protein_df):
        """Proteins with more than 1 membrane domain are excluded."""
        result = filter_ta_proteins(sample_protein_df)
        assert "P003" not in result["Entry"].values

    def test_excludes_n_terminal_tmd(self, sample_protein_df):
        """Proteins with N-terminal TMDs are excluded (P002 has N_term_md=1)."""
        result = filter_ta_proteins(sample_protein_df)
        assert "P002" not in result["Entry"].values

    def test_excludes_long_cterm_extension(self, sample_protein_df):
        """Proteins with cterm_distance > max_cterm are excluded (P004 = 50)."""
        result = filter_ta_proteins(sample_protein_df)
        assert "P004" not in result["Entry"].values

    def test_custom_max_cterm(self, sample_protein_df):
        """Raising max_cterm to 50 should include P004."""
        result = filter_ta_proteins(sample_protein_df, max_cterm=50)
        assert "P004" in result["Entry"].values

    def test_custom_tmd_count(self, sample_protein_df):
        """Setting tmd_count=2 should return only P003."""
        result = filter_ta_proteins(sample_protein_df, tmd_count=2, max_n_term=0, max_cterm=30)
        assert set(result["Entry"]) == {"P003"}

    def test_returns_copy(self, sample_protein_df):
        """Modifying result must not affect the original DataFrame."""
        result = filter_ta_proteins(sample_protein_df)
        result["cterm_distance"] = -1
        assert sample_protein_df["cterm_distance"].min() > 0

    def test_numeric_coercion(self, sample_protein_df):
        """String values in numeric columns should be coerced without error."""
        df = sample_protein_df.copy()
        df["cterm_distance"] = df["cterm_distance"].astype(str)
        df["membrane_domain_count"] = df["membrane_domain_count"].astype(str)
        df["N_term_md"] = df["N_term_md"].astype(str)
        result = filter_ta_proteins(df)
        assert len(result) == 2


# ── filter_disprot_human ──────────────────────────────────────────────────────

class TestFilterDisprotHuman:
    def test_keeps_human_entries_only(self, sample_disprot_df):
        result = filter_disprot_human(sample_disprot_df)
        assert all(result["NCBI Taxon ID"].astype(str).str.strip() == "9606")

    def test_excludes_non_human(self, sample_disprot_df):
        result = filter_disprot_human(sample_disprot_df)
        assert "DP00002" not in result["DisProt ID"].values

    def test_empty_df_returns_empty(self):
        empty = pd.DataFrame(columns=["NCBI Taxon ID", "DisProt ID"])
        result = filter_disprot_human(empty)
        assert len(result) == 0


# ── match_disprot ─────────────────────────────────────────────────────────────

class TestMatchDisprot:
    def test_matching_entry_annotated(self, sample_protein_df, sample_disprot_df):
        human_dp = filter_disprot_human(sample_disprot_df)
        ta = filter_ta_proteins(sample_protein_df)
        result = match_disprot(ta, human_dp)
        p001 = result[result["Entry"] == "P001"].iloc[0]
        assert p001["disprot_id"] == "DP00001"

    def test_non_matching_entry_is_nan(self, sample_protein_df, sample_disprot_df):
        human_dp = filter_disprot_human(sample_disprot_df)
        ta = filter_ta_proteins(sample_protein_df)
        result = match_disprot(ta, human_dp)
        p005 = result[result["Entry"] == "P005"].iloc[0]
        assert pd.isna(p005["disprot_id"])

    def test_multiple_regions_aggregated(self, sample_protein_df, sample_disprot_df):
        """Two DisProt region rows for P001 should be joined with ';'."""
        human_dp = filter_disprot_human(sample_disprot_df)
        ta = filter_ta_proteins(sample_protein_df)
        result = match_disprot(ta, human_dp)
        p001 = result[result["Entry"] == "P001"].iloc[0]
        assert "DP00001r001" in p001["disprot_disorder_regions"]
        assert "DP00001r002" in p001["disprot_disorder_regions"]

    def test_added_columns_present(self, sample_protein_df, sample_disprot_df):
        human_dp = filter_disprot_human(sample_disprot_df)
        ta = filter_ta_proteins(sample_protein_df)
        result = match_disprot(ta, human_dp)
        for col in ("disprot_id", "disprot_disorder_content", "disprot_disorder_regions"):
            assert col in result.columns

    def test_row_count_unchanged(self, sample_protein_df, sample_disprot_df):
        """Left-join must not add extra rows to the TA DataFrame."""
        human_dp = filter_disprot_human(sample_disprot_df)
        ta = filter_ta_proteins(sample_protein_df)
        result = match_disprot(ta, human_dp)
        assert len(result) == len(ta)


# ── aiupred_summary ───────────────────────────────────────────────────────────

class TestAiupredSummary:
    def test_empty_scores_returns_none(self):
        result = aiupred_summary([])
        assert result["aiupred_mean_score"] is None
        assert result["aiupred_disordered_fraction"] is None

    def test_mean_score_correct(self):
        result = aiupred_summary([0.2, 0.4, 0.6, 0.8])
        assert abs(result["aiupred_mean_score"] - 0.5) < 1e-4

    def test_disordered_fraction_above_threshold(self):
        # 3 of 4 residues ≥ 0.5
        result = aiupred_summary([0.1, 0.5, 0.7, 0.9])
        assert abs(result["aiupred_disordered_fraction"] - 0.75) < 1e-4

    def test_custom_threshold(self):
        result = aiupred_summary([0.3, 0.6, 0.9], threshold=0.6)
        assert abs(result["aiupred_disordered_fraction"] - (2 / 3)) < 1e-4

    def test_all_ordered(self):
        result = aiupred_summary([0.1, 0.2, 0.3])
        assert result["aiupred_disordered_fraction"] == 0.0

    def test_all_disordered(self):
        result = aiupred_summary([0.8, 0.9, 1.0])
        assert result["aiupred_disordered_fraction"] == 1.0


# ── mobidb_lite_summary ───────────────────────────────────────────────────────

class TestMobidbLiteSummary:
    def test_empty_data_returns_none(self):
        result = mobidb_lite_summary({}, 200)
        assert result["mobidb_lite_disordered_fraction"] is None

    def test_zero_length_returns_none(self):
        data = {"prediction-disorder-mobidb_lite": {"regions": [[1, 50]]}}
        result = mobidb_lite_summary(data, 0)
        assert result["mobidb_lite_disordered_fraction"] is None

    def test_correct_fraction(self):
        # 50 disordered residues out of 200
        data = {"prediction-disorder-mobidb_lite": {"regions": [[1, 50]]}}
        result = mobidb_lite_summary(data, 200)
        assert abs(result["mobidb_lite_disordered_fraction"] - 0.25) < 1e-4

    def test_multiple_regions(self):
        # 20 + 30 = 50 residues disordered out of 100
        data = {"prediction-disorder-mobidb_lite": {"regions": [[1, 20], [50, 79]]}}
        result = mobidb_lite_summary(data, 100)
        assert abs(result["mobidb_lite_disordered_fraction"] - 0.5) < 1e-4

    def test_no_regions_key_returns_none(self):
        data = {"prediction-disorder-mobidb_lite": {}}
        result = mobidb_lite_summary(data, 200)
        assert result["mobidb_lite_disordered_fraction"] is None

    def test_missing_prediction_key_returns_none(self):
        data = {"some_other_key": {}}
        result = mobidb_lite_summary(data, 200)
        assert result["mobidb_lite_disordered_fraction"] is None

    def test_fraction_capped_at_one(self):
        # Regions exceed protein length → fraction should not exceed 1.0
        data = {"prediction-disorder-mobidb_lite": {"regions": [[1, 250]]}}
        result = mobidb_lite_summary(data, 200)
        assert result["mobidb_lite_disordered_fraction"] <= 1.0


# ── run_disorder_predictions (mocked API calls) ───────────────────────────────

class TestRunDisorderPredictions:
    def _make_ta_df(self):
        return pd.DataFrame({
            "Entry":  ["P001", "P002"],
            "Length": [100, 200],
        })

    def test_columns_added(self):
        ta_df = self._make_ta_df()
        with (
            patch("idr_analysis.fetch_sequence", return_value="MKLV"),
            patch("idr_analysis.predict_aiupred",
                  return_value={"scores": [0.8, 0.2, 0.7, 0.9]}),
            patch("idr_analysis.predict_mobidb_lite",
                  return_value={"prediction-disorder-mobidb_lite": {"regions": [[1, 25]]}}),
        ):
            result = run_disorder_predictions(ta_df)

        for col in ("aiupred_mean_score", "aiupred_disordered_fraction",
                    "mobidb_lite_disordered_fraction"):
            assert col in result.columns

    def test_handles_missing_sequence(self):
        """If fetch_sequence returns None, AIUPred columns should be None."""
        ta_df = self._make_ta_df()
        with (
            patch("idr_analysis.fetch_sequence", return_value=None),
            patch("idr_analysis.predict_mobidb_lite", return_value={}),
        ):
            result = run_disorder_predictions(ta_df)

        assert result["aiupred_mean_score"].isna().all()
        assert result["aiupred_disordered_fraction"].isna().all()

    def test_row_count_preserved(self):
        ta_df = self._make_ta_df()
        with (
            patch("idr_analysis.fetch_sequence", return_value="MKLV"),
            patch("idr_analysis.predict_aiupred",
                  return_value={"scores": [0.5, 0.5, 0.5, 0.5]}),
            patch("idr_analysis.predict_mobidb_lite", return_value={}),
        ):
            result = run_disorder_predictions(ta_df)

        assert len(result) == len(ta_df)
