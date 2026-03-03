"""Tests for analysis/non_srp_filter and analysis/srp_filter modules."""

import pytest
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from analysis.non_srp_filter import analyze_non_srp_by_cterm_distance, count_domains as non_srp_count_domains
from analysis.srp_filter import analyze_srp_by_cterm_distance, count_domains as srp_count_domains


@pytest.fixture()
def sample_protein_df() -> pd.DataFrame:
    """Minimal protein DataFrame with required columns."""
    return pd.DataFrame({
        'UniProtID': ['P001', 'P002', 'P003', 'P004', 'P005'],
        'GeneName': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
        'Length': [200, 300, 150, 400, 100],
        'Transmembrane': ['170..190', '50..70', float('nan'), '380..395', '80..95'],
        'Intramembrane': [float('nan'), '280..295', '130..145', float('nan'), float('nan')],
        'Membrane_Domain_Count': [1, 2, 1, 1, 1],
    })


class TestCountDomains:
    def test_non_srp_single_domain(self):
        assert non_srp_count_domains('50..100') == 1

    def test_non_srp_multiple_domains(self):
        assert non_srp_count_domains('50..100 150..200') == 2

    def test_non_srp_nan_returns_zero(self):
        assert non_srp_count_domains(float('nan')) == 0

    def test_srp_single_domain(self):
        assert srp_count_domains('50..100') == 1

    def test_srp_multiple_domains(self):
        assert srp_count_domains('50..100 150..200') == 2

    def test_srp_nan_returns_zero(self):
        assert srp_count_domains(float('nan')) == 0

    def test_srp_and_non_srp_consistent(self):
        test_str = '10..50 100..150 200..250'
        assert srp_count_domains(test_str) == non_srp_count_domains(test_str)


class TestAnalyzeNonSrpByCtermDistance:
    def test_returns_dataframe(self, sample_protein_df):
        result = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert isinstance(result, pd.DataFrame)

    def test_adds_cterm_distance_column(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert 'cterm_distance' in result.columns

    def test_filters_by_threshold(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert all(result['cterm_distance'] <= 30)

    def test_does_not_mutate_input(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        original_cols = list(sample_protein_df.columns)
        analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert list(sample_protein_df.columns) == original_cols

    def test_custom_threshold_returns_more_rows(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result_30 = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        result_100 = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=100)
        assert len(result_100) >= len(result_30)


class TestAnalyzeSrpByCtermDistance:
    def test_returns_dataframe(self, sample_protein_df):
        result = analyze_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert isinstance(result, pd.DataFrame)

    def test_adds_cterm_distance_column(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = analyze_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert 'cterm_distance' in result.columns

    def test_filters_by_threshold(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = analyze_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert all(result['cterm_distance'] <= 30)

    def test_does_not_mutate_input(self, sample_protein_df, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        original_cols = list(sample_protein_df.columns)
        analyze_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert list(sample_protein_df.columns) == original_cols

    def test_result_same_as_non_srp(self, sample_protein_df, tmp_path, monkeypatch):
        """Both functions compute the same C-terminal distance filter."""
        monkeypatch.chdir(tmp_path)
        result_srp = analyze_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        result_non_srp = analyze_non_srp_by_cterm_distance(sample_protein_df, cterm_threshold=30)
        assert len(result_srp) == len(result_non_srp)
