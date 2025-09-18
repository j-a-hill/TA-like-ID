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
import numpy as np
from protein_analysis_utils import (
    filter_and_compare,
    filter_by_cterm_distance,
    quick_cterm_analysis,
    summary_report
)


def demonstrate_fixed_workflow():
    """Demonstrate the fixed category comparison workflow."""
    
    print("=" * 80)
    print("DEMONSTRATION: Fixed Category Comparison Workflow")
    print("=" * 80)
    print()
    
    # Load sample data
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
    print("❌ Complex workflow for category comparisons")
    print()
    
    print("AFTER THE FIX:")
    print("-" * 50)
    print("✅ Easy filtering by cterm_distance with comparison within subset")
    print("✅ Simple MD count analysis for C-terminal proteins")
    print("✅ Easy threshold testing and comparison")
    print("✅ Streamlined category comparison workflow")
    print()
    
    # ================================
    # DEMONSTRATION 1: Main issue fix
    # ================================
    
    print("DEMONSTRATION 1: Filter by cterm_distance, compare MD count by localization")
    print("=" * 80)
    
    # The main issue: Filter by cterm_distance <= 30, then compare MD counts
    filtered_df, md_counts = filter_by_cterm_distance(
        df, max_distance=30, compare_by='membrane_domain_count'
    )
    
    print(f"Step 1: Filter by cterm_distance <= 30")
    print(f"        Result: {len(filtered_df):,} proteins ({len(filtered_df)/len(df)*100:.1f}% of total)")
    print()
    
    print("Step 2: Compare membrane domain counts within filtered set")
    print("        Top results:")
    for count, freq in md_counts.head(8).items():
        pct = freq / len(filtered_df) * 100
        print(f"          {count:2d} domains: {freq:3d} proteins ({pct:5.1f}%)")
    print()
    
    # Now compare localization within the same filtered set
    loc_counts = filtered_df['Reduced.CC.Terms'].value_counts()
    print("Step 3: Compare localization within the same filtered set")
    print("        Top localizations:")
    for loc, count in loc_counts.head(10).items():
        pct = count / len(filtered_df) * 100
        print(f"          {count:3d} ({pct:4.1f}%): {loc}")
    print()
    
    # ================================
    # DEMONSTRATION 2: Easy threshold testing
    # ================================
    
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
    
    # ================================
    # DEMONSTRATION 3: Advanced comparisons  
    # ================================
    
    print("DEMONSTRATION 3: Advanced category comparisons")
    print("=" * 80)
    
    # Filter by one criterion, compare by multiple others
    close_proteins = df[df['cterm_distance'] <= 30].copy()
    
    print("For proteins with cterm_distance <= 30:")
    print()
    
    # Compare prediction types
    pred_dist = close_proteins['Prediction'].value_counts()
    print("Prediction type distribution:")
    for pred, count in pred_dist.items():
        pct = count / len(close_proteins) * 100
        print(f"  {pred:10s}: {count:4d} ({pct:5.1f}%)")
    print()
    
    # Compare BioGRID status
    biogrid_dist = close_proteins['in_biogrid'].value_counts()
    print("BioGRID interaction status:")
    for status, count in biogrid_dist.items():
        status_label = "In BioGRID" if status else "Not in BioGRID"
        pct = count / len(close_proteins) * 100
        print(f"  {status_label:15s}: {count:4d} ({pct:5.1f}%)")
    print()
    
    # Compare mass spec status
    ms_dist = close_proteins['in_massspec'].value_counts()
    print("Mass spectrometry detection:")
    for status, count in ms_dist.items():
        status_label = "Detected in MS" if status else "Not in MS"
        pct = count / len(close_proteins) * 100
        print(f"  {status_label:15s}: {count:4d} ({pct:5.1f}%)")
    print()
    
    # ================================
    # DEMONSTRATION 4: Comprehensive reporting
    # ================================
    
    print("DEMONSTRATION 4: Comprehensive reporting")
    print("=" * 80)
    
    # Generate summary report
    summary_report(
        df,
        filter_column='cterm_distance',
        filter_value=30,
        operator='<=',
        compare_columns=['membrane_domain_count', 'Prediction', 'in_biogrid', 'in_massspec']
    )
    
    # ================================
    # SUMMARY
    # ================================
    
    print("=" * 80)
    print("SUMMARY: Category Comparison Issues FIXED")
    print("=" * 80)
    print()
    print("✅ FIXED: Can now filter by cterm_distance and compare within subset")
    print("✅ FIXED: Easy MD count comparison by localization for filtered proteins")
    print("✅ FIXED: Simple threshold testing across multiple distance values")
    print("✅ FIXED: Streamlined workflow for complex category comparisons")
    print("✅ FIXED: Comprehensive reporting and data export functionality")
    print()
    print("The improved workflow provides:")
    print("  • filter_and_compare() function for general filtering and comparison")
    print("  • filter_by_cterm_distance() for the specific C-terminal use case")
    print("  • summary_report() for comprehensive analysis")
    print("  • Easy integration into existing analysis scripts")
    print()
    print("This addresses all issues mentioned in GitHub issue #3.")


if __name__ == "__main__":
    demonstrate_fixed_workflow()