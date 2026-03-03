"""Tests for protein_analysis_utils module."""

import pytest
import pandas as pd

from protein_analysis_utils import (
    FilterConfig,
    filter_and_compare,
    filter_by_cterm_distance,
    compare_categories,
    summary_report,
    analyze_cterm_distance_effects,
    cross_tabulate_categories,
    calc_min_cterm_distance,
)


@pytest.fixture()
def sample_df() -> pd.DataFrame:
    """Create a small synthetic DataFrame that mirrors the real data schema."""
    return pd.DataFrame({
        "cterm_distance": [5, 15, 25, 35, 50, 100],
        "membrane_domain_count": [1, 1, 2, 1, 3, 2],
        "Prediction": ["SP", "OTHER", "SP", "OTHER", "SP", "OTHER"],
        "in_biogrid": [True, False, True, False, True, False],
        "in_massspec": [False, True, False, True, False, True],
        "Reduced.CC.Terms": ["ER", "Golgi", "ER", "plasma membrane", "ER", "Golgi"],
    })


class TestFilterAndCompare:
    def test_lte_operator_returns_filtered_subset(self, sample_df):
        filtered, counts = filter_and_compare(sample_df, "cterm_distance", 30, "membrane_domain_count")
        assert len(filtered) == 3  # rows with cterm_distance <= 30: 5, 15, 25
        assert all(filtered["cterm_distance"] <= 30)

    def test_lt_operator(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "cterm_distance", 30, "membrane_domain_count", operator="<")
        assert len(filtered) == 3  # rows strictly < 30

    def test_gte_operator(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "cterm_distance", 50, "membrane_domain_count", operator=">=")
        assert len(filtered) == 2

    def test_gt_operator(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "cterm_distance", 50, "membrane_domain_count", operator=">")
        assert len(filtered) == 1

    def test_eq_operator(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "membrane_domain_count", 1, "Prediction", operator="==")
        assert len(filtered) == 3

    def test_ne_operator(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "membrane_domain_count", 1, "Prediction", operator="!=")
        assert len(filtered) == 3

    def test_unsupported_operator_raises(self, sample_df):
        with pytest.raises(ValueError, match="Unsupported operator"):
            filter_and_compare(sample_df, "cterm_distance", 30, "membrane_domain_count", operator="~")

    def test_counts_sorted_by_default(self, sample_df):
        _, counts = filter_and_compare(sample_df, "cterm_distance", 100, "membrane_domain_count")
        assert counts.iloc[0] >= counts.iloc[-1]

    def test_top_n_limits_results(self, sample_df):
        _, counts = filter_and_compare(sample_df, "cterm_distance", 100, "membrane_domain_count", top_n=1)
        assert len(counts) == 1

    def test_returns_copy_not_view(self, sample_df):
        filtered, _ = filter_and_compare(sample_df, "cterm_distance", 30, "membrane_domain_count")
        filtered["cterm_distance"] = -1
        assert sample_df["cterm_distance"].min() > 0


class TestFilterByCTermDistance:
    def test_default_threshold(self, sample_df):
        filtered, _ = filter_by_cterm_distance(sample_df)
        assert all(filtered["cterm_distance"] <= 30)

    def test_custom_threshold(self, sample_df):
        filtered, _ = filter_by_cterm_distance(sample_df, max_distance=20)
        assert all(filtered["cterm_distance"] <= 20)

    def test_compare_by_column(self, sample_df):
        _, counts = filter_by_cterm_distance(sample_df, max_distance=100, compare_by="Prediction")
        assert set(counts.index).issubset({"SP", "OTHER"})


class TestCompareCategories:
    def test_simple_value_counts(self, sample_df):
        result = compare_categories(sample_df, "Prediction")
        assert isinstance(result, pd.Series)
        assert result["SP"] == 3
        assert result["OTHER"] == 3

    def test_normalized(self, sample_df):
        result = compare_categories(sample_df, "Prediction", normalize=True)
        assert abs(result.sum() - 1.0) < 1e-9

    def test_group_by(self, sample_df):
        result = compare_categories(sample_df, "Prediction", group_by="membrane_domain_count")
        assert isinstance(result, pd.DataFrame)

    def test_group_by_normalized(self, sample_df):
        result = compare_categories(sample_df, "Prediction", group_by="membrane_domain_count", normalize=True)
        assert isinstance(result, pd.DataFrame)


class TestSummaryReport:
    def test_runs_without_error(self, sample_df, capsys):
        summary_report(sample_df, FilterConfig(filter_column="cterm_distance", filter_value=30, operator="<="))
        captured = capsys.readouterr()
        assert "Summary Report" in captured.out
        assert "Filter:" in captured.out

    def test_custom_compare_columns(self, sample_df, capsys):
        summary_report(
            sample_df,
            FilterConfig(
                filter_column="cterm_distance",
                filter_value=100,
                compare_columns=["Prediction"],
            ),
        )
        captured = capsys.readouterr()
        assert "Prediction" in captured.out


class TestAnalyzeCTermDistanceEffects:
    def test_returns_dataframe(self, sample_df):
        result = analyze_cterm_distance_effects(sample_df, distance_thresholds=[20, 50])
        assert isinstance(result, pd.DataFrame)

    def test_columns_match_thresholds(self, sample_df):
        result = analyze_cterm_distance_effects(sample_df, distance_thresholds=[10, 30])
        assert "<= 10" in result.columns
        assert "<= 30" in result.columns

    def test_default_thresholds(self, sample_df):
        result = analyze_cterm_distance_effects(sample_df)
        assert len(result.columns) == 5


class TestCrossTabulateCategories:
    def test_returns_crosstab(self, sample_df):
        result = cross_tabulate_categories(sample_df, "Prediction", "in_biogrid")
        assert isinstance(result, pd.DataFrame)
        assert "All" in result.index  # margins=True adds All row

    def test_with_filter(self, sample_df):
        result = cross_tabulate_categories(
            sample_df, "Prediction", "in_biogrid",
            filter_column="cterm_distance", filter_value=30,
        )
        assert isinstance(result, pd.DataFrame)

    def test_without_filter_uses_full_df(self, sample_df):
        result = cross_tabulate_categories(sample_df, "Prediction", "in_biogrid")
        # All rows total should equal len(sample_df)
        assert result.loc["All", "All"] == len(sample_df)


class TestFilterConfig:
    def test_defaults(self):
        config = FilterConfig()
        assert config.filter_column == 'cterm_distance'
        assert config.filter_value == 30
        assert config.operator == '<='
        assert 'membrane_domain_count' in config.compare_columns

    def test_custom_values(self):
        config = FilterConfig(filter_column="membrane_domain_count", filter_value=2, operator=">=")
        assert config.filter_column == "membrane_domain_count"
        assert config.filter_value == 2
        assert config.operator == ">="

    def test_compare_columns_are_independent_per_instance(self):
        c1 = FilterConfig()
        c2 = FilterConfig()
        c1.compare_columns.append("extra_col")
        assert "extra_col" not in c2.compare_columns


class TestSummaryReportWithConfig:
    def test_config_object_used(self, sample_df, capsys):
        config = FilterConfig(
            filter_column="cterm_distance",
            filter_value=100,
            operator="<=",
            compare_columns=["Prediction"],
        )
        summary_report(sample_df, config)
        captured = capsys.readouterr()
        assert "Prediction" in captured.out

    def test_default_config_used_when_none(self, sample_df, capsys):
        summary_report(sample_df)
        captured = capsys.readouterr()
        assert "Summary Report" in captured.out

    def test_config_filter_value_applied(self, sample_df, capsys):
        config = FilterConfig(filter_value=100)
        summary_report(sample_df, config)
        captured = capsys.readouterr()
        # filter_value=100 means all 6 rows are retained
        assert "6" in captured.out


class TestCalcMinCtermDistance:
    def test_single_domain(self):
        row = pd.Series({'Length': 200, 'Transmembrane': '170..190'})
        assert calc_min_cterm_distance(row) == 200 - 190

    def test_multiple_domains_picks_closest(self):
        row = pd.Series({'Length': 200, 'Transmembrane': '50..70 170..190'})
        assert calc_min_cterm_distance(row) == 200 - 190

    def test_no_domain_returns_none(self):
        row = pd.Series({'Length': 200, 'Transmembrane': float('nan')})
        assert calc_min_cterm_distance(row) is None

    def test_missing_length_returns_none(self):
        row = pd.Series({'Length': float('nan'), 'Transmembrane': '50..70'})
        assert calc_min_cterm_distance(row) is None
