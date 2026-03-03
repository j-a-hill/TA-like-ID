"""Tests for analysis/Uniprot_TMD_search module."""

import pytest
import pandas as pd
import sys
from pathlib import Path

# Add repo root so analysis/ modules can import protein_analysis_utils
sys.path.insert(0, str(Path(__file__).parent.parent))

from analysis.Uniprot_TMD_search import (
    extract_near_c_terminus_domains,
    count_domains,
    build_gene_set,
    token_match,
    calculate_cterm_distance,
)
from analysis.tmd_count import count_domains as tmd_count_domains


@pytest.fixture()
def sample_uniprot_df() -> pd.DataFrame:
    """Minimal synthetic Uniprot DataFrame."""
    return pd.DataFrame({
        'Entry': ['P001', 'P002', 'P003', 'P004'],
        'Gene Names': ['GENE1 GN1', 'GENE2', 'GENE3', 'GENE4'],
        'Length': [200, 300, 150, 400],
        'Transmembrane': ['170..190', '50..70', float('nan'), '380..395'],
        'Intramembrane': [float('nan'), '280..295', '130..145', float('nan')],
    })


@pytest.fixture()
def sample_biogrid_df() -> pd.DataFrame:
    """Minimal synthetic BioGRID DataFrame."""
    return pd.DataFrame({
        'official_symbol_for_interactor_a': ['GENE1', 'GENE5'],
        'official_symbol_for_interactor_b': ['GENE2', 'GENE6'],
        'synonyms/aliases_for_interactor_a': ['GN1', 'GN5'],
        'synonyms/aliases_for_interactor_b': ['GN2', 'GN6'],
    })


class TestExtractNearCTerminusDomains:
    def test_filters_transmembrane_within_threshold(self, sample_uniprot_df):
        # P001: TM 200-190=10 ✓; P002: IM 300-295=5 ✓; P003: IM 150-145=5 ✓; P004: TM 400-395=5 ✓
        # P002's Transmembrane domain ends at 70: 300-70=230 > 30, but its Intramembrane at 295: 5 <= 30
        result = extract_near_c_terminus_domains(sample_uniprot_df, N=30)
        entries = set(result['Entry'])
        assert 'P001' in entries  # 200-190=10 <= 30
        assert 'P002' in entries  # 300-295=5 (intramembrane) <= 30
        assert 'P003' in entries  # 150-145=5 <= 30 (intramembrane)
        assert 'P004' in entries  # 400-395=5 <= 30

    def test_empty_df_returns_empty(self):
        empty_df = pd.DataFrame(columns=['Entry', 'Length', 'Transmembrane', 'Intramembrane'])
        result = extract_near_c_terminus_domains(empty_df, N=30)
        assert len(result) == 0

    def test_custom_threshold(self, sample_uniprot_df):
        result = extract_near_c_terminus_domains(sample_uniprot_df, N=5)
        # Only domains with distance <= 5
        entries = set(result['Entry'])
        assert 'P001' not in entries  # 200-190=10 > 5
        assert 'P003' in entries  # 150-145=5 <= 5
        assert 'P004' in entries  # 400-395=5 <= 5

    def test_na_domain_str_excluded(self, sample_uniprot_df):
        # P002 has no Transmembrane domain that's close (only P002's Intramembrane is 300-295=5)
        result = extract_near_c_terminus_domains(sample_uniprot_df, N=10)
        assert 'P002' in set(result['Entry'])  # via Intramembrane 300-295=5


class TestCountDomains:
    def test_single_domain(self):
        assert count_domains('50..100') == 1

    def test_multiple_domains(self):
        assert count_domains('50..100 150..200 250..300') == 3

    def test_nan_returns_zero(self):
        assert count_domains(float('nan')) == 0

    def test_empty_string(self):
        assert count_domains('') == 0

    def test_tmd_count_module_consistency(self):
        """tmd_count.count_domains should match Uniprot_TMD_search.count_domains."""
        test_str = '10..50 100..150'
        assert count_domains(test_str) == tmd_count_domains(test_str)


class TestBuildGeneSet:
    def test_returns_lowercase_set(self, sample_biogrid_df):
        gene_set = build_gene_set(sample_biogrid_df)
        assert 'gene1' in gene_set
        assert 'gene2' in gene_set
        assert 'gn1' in gene_set

    def test_missing_columns_ignored(self):
        df = pd.DataFrame({'official_symbol_for_interactor_a': ['BRCA1']})
        gene_set = build_gene_set(df)
        assert 'brca1' in gene_set

    def test_empty_df_returns_empty_set(self):
        df = pd.DataFrame({'official_symbol_for_interactor_a': pd.Series([], dtype=str)})
        gene_set = build_gene_set(df)
        assert isinstance(gene_set, set)


class TestTokenMatch:
    def test_match_found(self):
        assert token_match('GENE1 isoform 2', {'gene1', 'gene2'}) is True

    def test_no_match(self):
        assert token_match('GENE99 isoform 2', {'gene1', 'gene2'}) is False

    def test_nan_returns_false(self):
        assert token_match(float('nan'), {'gene1'}) is False

    def test_case_insensitive(self):
        assert token_match('Gene1', {'gene1'}) is True


class TestCalculateCtermDistance:
    def test_transmembrane_domain(self):
        row = pd.Series({'Length': 200, 'Transmembrane': '170..190', 'Intramembrane': float('nan')})
        assert calculate_cterm_distance(row) == 10  # 200 - 190

    def test_picks_minimum_across_domains(self):
        row = pd.Series({'Length': 200, 'Transmembrane': '50..80 170..190', 'Intramembrane': float('nan')})
        assert calculate_cterm_distance(row) == 10  # min(200-80=120, 200-190=10)

    def test_intramembrane_domain(self):
        row = pd.Series({'Length': 150, 'Transmembrane': float('nan'), 'Intramembrane': '130..145'})
        assert calculate_cterm_distance(row) == 5  # 150 - 145

    def test_both_domains_picks_minimum(self):
        row = pd.Series({'Length': 300, 'Transmembrane': '50..70', 'Intramembrane': '280..295'})
        # TM: 300-70=230, IM: 300-295=5 → min is 5
        assert calculate_cterm_distance(row) == 5

    def test_no_domain_returns_none(self):
        row = pd.Series({'Length': 200, 'Transmembrane': float('nan'), 'Intramembrane': float('nan')})
        assert calculate_cterm_distance(row) is None

    def test_missing_length_returns_none(self):
        row = pd.Series({'Length': float('nan'), 'Transmembrane': '50..70', 'Intramembrane': float('nan')})
        assert calculate_cterm_distance(row) is None
