#!/usr/bin/env python3
"""
Demonstration of Fixed Category Comparison Workflow

This script demonstrates how the improved analysis workflow fixes the 
category comparison issues mentioned in GitHub issue #3.

Before the fix: Users couldn't easily filter by cterm_distance and then 
compare within that set (e.g., MD count by localization).

After the fix: Users can easily perform multi-step filtering and comparisons.
"""

import pandas as pd
from protein_analysis_utils import (
    filter_and_compare,
    filter_by_cterm_distance,
    quick_cterm_analysis,
    summary_report
)


def _demo_cterm_filter_and_comparison(df: pd.DataFrame) -> None:
    """Demo 1: Filter by cterm_distance, compare MD count by localization."""
    print("DEMONSTRATION 1: Filter by cterm_distance, compare MD count by localization")
    print("=" * 80)

    filtered_df, md_counts = filter_by_cterm_distance(df, max_distance=30, compare_by='membrane_domain_count')

    print("Step 1: Filter by cterm_distance <= 30")
    print(f"        Result: {len(filtered_df):,} proteins ({len(filtered_df)/len(df)*100:.1f}% of total)")
    print()
    print("Step 2: Compare membrane domain counts within filtered set")
    print("        Top results:")
    for count, freq in md_counts.head(8).items():
        pct = freq / len(filtered_df) * 100
        print(f"          {count:2d} domains: {freq:3d} proteins ({pct:5.1f}%)")
    print()

    loc_counts = filtered_df['Reduced.CC.Terms'].value_counts()
    print("Step 3: Compare localization within the same filtered set")
    print("        Top localizations:")
    for loc, count in loc_counts.head(10).items():
        pct = count / len(filtered_df) * 100
        print(f"          {count:3d} ({pct:4.1f}%): {loc}")
    print()


def _demo_threshold_testing(df: pd.DataFrame) -> None:
    """Demo 2: Compare effects of different C-terminal distance thresholds."""
    print("DEMONSTRATION 2: Easy threshold testing")
    print("=" * 80)
    print("Compare effects of different C-terminal distance thresholds:")
    print()
    print("Threshold | Count | % of Total | Top MD Count")
    print("-" * 50)

    for threshold in [10, 20, 30, 50, 100]:
        subset_df, subset_md = filter_by_cterm_distance(df, threshold, 'membrane_domain_count')
        top_md = subset_md.index[0] if len(subset_md) > 0 else 'N/A'
        top_count = subset_md.iloc[0] if len(subset_md) > 0 else 0
        print(f"   <= {threshold:2d}  | {len(subset_df):5d} | {len(subset_df)/len(df)*100:6.1f}%   | {top_md} ({top_count} proteins)")
    print()


def _demo_advanced_comparisons(df: pd.DataFrame) -> None:
    """Demo 3: Advanced category comparisons for proteins near C-terminus."""
    print("DEMONSTRATION 3: Advanced category comparisons")
    print("=" * 80)

    close_proteins = df[df['cterm_distance'] <= 30].copy()
    print("For proteins with cterm_distance <= 30:")
    print()

    pred_dist = close_proteins['Prediction'].value_counts()
    print("Prediction type distribution:")
    for pred, count in pred_dist.items():
        pct = count / len(close_proteins) * 100
        print(f"  {pred:10s}: {count:4d} ({pct:5.1f}%)")
    print()

    biogrid_dist = close_proteins['in_biogrid'].value_counts()
    print("BioGRID interaction status:")
    for status, count in biogrid_dist.items():
        status_label = "In BioGRID" if status else "Not in BioGRID"
        pct = count / len(close_proteins) * 100
        print(f"  {status_label:15s}: {count:4d} ({pct:5.1f}%)")
    print()

    ms_dist = close_proteins['in_massspec'].value_counts()
    print("Mass spectrometry detection:")
    for status, count in ms_dist.items():
        status_label = "Detected in MS" if status else "Not in MS"
        pct = count / len(close_proteins) * 100
        print(f"  {status_label:15s}: {count:4d} ({pct:5.1f}%)")
    print()


def _demo_comprehensive_reporting(df: pd.DataFrame) -> None:
    """Demo 4: Comprehensive summary report."""
    print("DEMONSTRATION 4: Comprehensive reporting")
    print("=" * 80)

    summary_report(
        df,
        filter_column='cterm_distance',
        filter_value=30,
        operator='<=',
        compare_columns=['membrane_domain_count', 'Prediction', 'in_biogrid', 'in_massspec']
    )


def demonstrate_fixed_workflow() -> None:
    """Demonstrate the fixed category comparison workflow."""
    print("=" * 80)
    print("DEMONSTRATION: Fixed Category Comparison Workflow")
    print("=" * 80)
    print()

    try:
        df = pd.read_csv('raw_data/membrane_protein_analysis_with_reduced_cc.csv')
        print(f"✓ Loaded data: {len(df):,} membrane proteins")
    except FileNotFoundError:
        print("✗ Sample data not found. Please ensure membrane_protein_analysis_with_reduced_cc.csv is available.")
        return
    print()

    print("BEFORE THE FIX:")
    print("-" * 50)
    print("❌ Users had trouble filtering by cterm_distance and comparing within subset")
    print("❌ No easy way to compare MD count by localization for filtered proteins")
    print("❌ Difficult to test multiple distance thresholds")
    print()
    print("AFTER THE FIX:")
    print("-" * 50)
    print("✅ Easy filtering by cterm_distance with comparison within subset")
    print("✅ Simple MD count analysis for C-terminal proteins")
    print("✅ Streamlined category comparison workflow")
    print()

    _demo_cterm_filter_and_comparison(df)
    _demo_threshold_testing(df)
    _demo_advanced_comparisons(df)
    _demo_comprehensive_reporting(df)

    print("=" * 80)
    print("SUMMARY: Category Comparison Issues FIXED")
    print("=" * 80)
    print()
    print("✅ FIXED: Can now filter by cterm_distance and compare within subset")
    print("✅ FIXED: Easy MD count comparison by localization for filtered proteins")
    print("✅ FIXED: Simple threshold testing across multiple distance values")
    print("✅ FIXED: Streamlined workflow for complex category comparisons")
    print()
    print("The improved workflow provides:")
    print("  • filter_and_compare() function for general filtering and comparison")
    print("  • filter_by_cterm_distance() for the specific C-terminal use case")
    print("  • summary_report() for comprehensive analysis")
    print()
    print("This addresses all issues mentioned in GitHub issue #3.")


if __name__ == "__main__":
    demonstrate_fixed_workflow()

