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
    run_aiupred_local,
    run_mobidb_lite_local,
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


# ── run_aiupred_local ─────────────────────────────────────────────────────────

class TestRunAiupredLocal:
    """Tests for run_aiupred_local using a fake aiupred_lib module."""

    @staticmethod
    def _make_fake_aiupred_lib(tmp_path, scores: list):
        """Write a minimal fake aiupred_lib.py that returns *scores* from predict_disorder."""
        lib = tmp_path / "aiupred_lib.py"
        lib.write_text(
            "import os\n"
            "PATH = os.path.dirname(os.path.abspath(__file__))\n"
            f"_SCORES = {scores!r}\n"
            "def init_models(prediction_type, force_cpu=False, gpu_num=0):\n"
            "    return ('embed', 'reg', 'cpu')\n"
            "def predict_disorder(seq, emb, reg, dev, smoothing=True):\n"
            "    return _SCORES\n"
        )
        return lib

    def test_returns_scores_on_success(self, tmp_path):
        """Fake aiupred_lib with known scores → run_aiupred_local returns them."""
        # Clear cache so this test's tmp_path is freshly loaded
        import idr_analysis
        idr_analysis._aiupred_cache.clear()
        self._make_fake_aiupred_lib(tmp_path, [0.80, 0.23, 0.71, 0.55])
        result = run_aiupred_local("MKLV", aiupred_dir=str(tmp_path))
        assert result == {"scores": [0.80, 0.23, 0.71, 0.55]}

    def test_empty_dir_returns_empty(self):
        result = run_aiupred_local("MKLV", aiupred_dir="")
        assert result == {}

    def test_missing_dir_returns_empty(self):
        result = run_aiupred_local("MKLV", aiupred_dir="/nonexistent/path")
        assert result == {}

    def test_missing_lib_returns_empty(self, tmp_path):
        """No aiupred_lib.py in the dir → {}."""
        import idr_analysis
        idr_analysis._aiupred_cache.clear()
        result = run_aiupred_local("MKLV", aiupred_dir=str(tmp_path))
        assert result == {}

    def test_lib_import_error_returns_empty(self, tmp_path):
        """aiupred_lib.py raises ImportError (e.g. torch missing) → {}."""
        import idr_analysis
        idr_analysis._aiupred_cache.clear()
        (tmp_path / "aiupred_lib.py").write_text("raise ImportError('torch not found')\n")
        result = run_aiupred_local("MKLV", aiupred_dir=str(tmp_path))
        assert result == {}

    def test_models_cached_across_calls(self, tmp_path):
        """init_models should only be called once even for multiple sequences."""
        import idr_analysis
        idr_analysis._aiupred_cache.clear()

        lib = tmp_path / "aiupred_lib.py"
        lib.write_text(
            "import os\n"
            "PATH = os.path.dirname(os.path.abspath(__file__))\n"
            "_init_call_count = 0\n"
            "def init_models(t, force_cpu=False, gpu_num=0):\n"
            "    global _init_call_count\n"
            "    _init_call_count += 1\n"
            "    return ('e', 'r', 'cpu')\n"
            "def predict_disorder(s, e, r, d, smoothing=True):\n"
            "    return [0.5] * len(s)\n"
        )
        run_aiupred_local("AAA", aiupred_dir=str(tmp_path))
        run_aiupred_local("CCC", aiupred_dir=str(tmp_path))

        cached_mod = idr_analysis._aiupred_cache[str(tmp_path)]
        # init_models was called exactly once (during _load_aiupred_lib)
        assert cached_mod._init_call_count == 1


# ── run_mobidb_lite_local ─────────────────────────────────────────────────────

class TestRunMobidbLiteLocal:
    def test_empty_bindir_returns_empty(self):
        result = run_mobidb_lite_local("MKLV", bindir="")
        assert result == {}

    def test_nonexistent_bindir_returns_empty(self):
        result = run_mobidb_lite_local("MKLV", bindir="/nonexistent/bin")
        assert result == {}

    def test_returns_regions_on_success(self, tmp_path):
        """Mock mobidb_lite.consensus.run to yield a plausible result."""
        fake_regions = {"mobidblite": [(1, 40), (80, 120)]}
        fake_scores = {"mobidblite": [0.7] * 200}

        with patch("mobidb_lite.consensus.run",
                   return_value=iter([("query", fake_regions, fake_scores)])):
            result = run_mobidb_lite_local("M" * 200, bindir=str(tmp_path))

        assert "prediction-disorder-mobidb_lite" in result
        assert result["prediction-disorder-mobidb_lite"]["regions"] == [(1, 40), (80, 120)]

    def test_none_regions_returns_empty(self, tmp_path):
        """When the predictor returns None regions, return empty dict."""
        with patch("mobidb_lite.consensus.run",
                   return_value=iter([("query", None, {})])):
            result = run_mobidb_lite_local("MKLV", bindir=str(tmp_path))
        assert result == {}

    def test_import_error_returns_empty(self, tmp_path):
        """If mobidb_lite is not installed, return empty dict gracefully."""
        with patch.dict("sys.modules", {"mobidb_lite.consensus": None}):
            result = run_mobidb_lite_local("MKLV", bindir=str(tmp_path))
        assert result == {}


# ── run_disorder_predictions (mocked local tools) ─────────────────────────────

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
            patch("idr_analysis.run_aiupred_local",
                  return_value={"scores": [0.8, 0.2, 0.7, 0.9]}),
            patch("idr_analysis.run_mobidb_lite_local",
                  return_value={"prediction-disorder-mobidb_lite": {"regions": [[1, 25]]}}),
        ):
            result = run_disorder_predictions(ta_df)

        for col in ("aiupred_mean_score", "aiupred_disordered_fraction",
                    "mobidb_lite_disordered_fraction"):
            assert col in result.columns

    def test_handles_missing_sequence(self):
        """If fetch_sequence returns None, all prediction columns should be None."""
        ta_df = self._make_ta_df()
        with patch("idr_analysis.fetch_sequence", return_value=None):
            result = run_disorder_predictions(ta_df)

        assert result["aiupred_mean_score"].isna().all()
        assert result["aiupred_disordered_fraction"].isna().all()
        assert result["mobidb_lite_disordered_fraction"].isna().all()

    def test_row_count_preserved(self):
        ta_df = self._make_ta_df()
        with (
            patch("idr_analysis.fetch_sequence", return_value="MKLV"),
            patch("idr_analysis.run_aiupred_local",
                  return_value={"scores": [0.5, 0.5, 0.5, 0.5]}),
            patch("idr_analysis.run_mobidb_lite_local", return_value={}),
        ):
            result = run_disorder_predictions(ta_df)

        assert len(result) == len(ta_df)

    def test_sequence_fetched_once_per_protein(self):
        """fetch_sequence should be called exactly once per row (used for both tools)."""
        ta_df = self._make_ta_df()
        with (
            patch("idr_analysis.fetch_sequence", return_value="MKLV") as mock_fetch,
            patch("idr_analysis.run_aiupred_local", return_value={}),
            patch("idr_analysis.run_mobidb_lite_local", return_value={}),
        ):
            run_disorder_predictions(ta_df)

        assert mock_fetch.call_count == len(ta_df)
