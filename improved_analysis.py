#!/usr/bin/env python3
"""
Improved Membrane Protein Analysis

This script provides the fixed workflow for category comparisons that addresses
the issues mentioned in the GitHub issue. Users can now easily:
1. Filter by cterm_distance  
2. Compare MD count by localization within that filtered set
3. Perform other complex filtering and comparison operations

This replaces the problematic workflow from the Jupyter notebook.
"""

import pandas as pd
import sys
from pathlib import Path

# Import our custom utilities
from protein_analysis_utils import (
    filter_and_compare,
    filter_by_cterm_distance, 
    summary_report,
    analyze_cterm_distance_effects,
    cross_tabulate_categories,
    quick_cterm_analysis,
    quick_localization_analysis
)


def main():
    """Main analysis workflow with improved category comparisons."""
    
    # Load the data
    data_file = 'raw_data/membrane_protein_analysis_with_reduced_cc.csv'
    if not Path(data_file).exists():
        print(f"Error: Data file {data_file} not found.")
        print("Please ensure the membrane protein analysis data is available.")
        return
    
    df = pd.read_csv(data_file)
    print(f"Loaded membrane protein data: {len(df):,} proteins")
    print(f"Columns: {', '.join(df.columns)}")
    print()
    
    # ================================
    # FIXED WORKFLOW: Filter by cterm_distance, then compare within subset
    # ================================
    
    print("=== IMPROVED CATEGORY COMPARISON WORKFLOW ===")
    print()
    
    # Example 1: Filter by cterm_distance <= 30, then compare MD counts
    print("1. Filter by cterm_distance <= 30, compare membrane domain counts:")
    print("-" * 70)
    filtered_df, md_counts = filter_by_cterm_distance(df, max_distance=30, compare_by='membrane_domain_count')
    
    print(f"Found {len(filtered_df):,} proteins within 30 residues of C-terminus")
    print("Membrane domain count distribution:")
    for count, frequency in md_counts.head(10).items():
        percentage = frequency / len(filtered_df) * 100
        print(f"  {count:2d} domains: {frequency:4d} proteins ({percentage:5.1f}%)")
    print()
    
    # Example 2: Filter by cterm_distance <= 50, then compare by prediction type  
    print("2. Filter by cterm_distance <= 50, compare by prediction type:")
    print("-" * 70)
    filtered_df2, pred_counts = filter_and_compare(
        df, 'cterm_distance', 50, 'Prediction', operator='<='
    )
    
    print(f"Found {len(filtered_df2):,} proteins within 50 residues of C-terminus")
    print("Prediction type distribution:")
    for pred_type, count in pred_counts.items():
        percentage = count / len(filtered_df2) * 100
        print(f"  {pred_type:10s}: {count:4d} proteins ({percentage:5.1f}%)")
    print()
    
    # Example 3: Filter by membrane domain count, then compare localization
    print("3. Filter by membrane_domain_count == 1, compare localization:")
    print("-" * 70)
    filtered_df3, loc_counts = filter_and_compare(
        df, 'membrane_domain_count', 1, 'Reduced.CC.Terms', operator='=='
    )
    
    print(f"Found {len(filtered_df3):,} proteins with exactly 1 membrane domain")
    print("Top localizations:")
    for loc, count in loc_counts.head(15).items():
        percentage = count / len(filtered_df3) * 100
        print(f"  {count:3d} ({percentage:4.1f}%): {loc}")
    print()
    
    # ================================
    # ADVANCED COMPARISONS  
    # ================================
    
    print("=== ADVANCED ANALYSIS EXAMPLES ===")
    print()
    
    # Cross-tabulation: Compare prediction type vs in_biogrid for close C-term proteins
    print("4. Cross-tabulation: Prediction type vs BioGRID status (cterm_distance <= 30):")
    print("-" * 80)
    crosstab = cross_tabulate_categories(
        df, 'Prediction', 'in_biogrid', 
        filter_column='cterm_distance', filter_value=30
    )
    print(crosstab)
    print()
    
    # Effect of different distance thresholds on membrane domain distribution
    print("5. Effect of C-terminal distance thresholds on membrane domain counts:")
    print("-" * 80)
    threshold_analysis = analyze_cterm_distance_effects(
        df, distance_thresholds=[10, 20, 30, 50, 100], compare_column='membrane_domain_count'
    )
    print(threshold_analysis.head(15))
    print()
    
    # ================================
    # COMPREHENSIVE SUMMARY REPORTS
    # ================================
    
    print("=== COMPREHENSIVE SUMMARY REPORTS ===")
    print()
    
    # Summary for proteins close to C-terminus
    summary_report(
        df, 
        filter_column='cterm_distance', 
        filter_value=30, 
        operator='<=',
        compare_columns=['membrane_domain_count', 'Prediction', 'in_biogrid', 'in_massspec']
    )
    
    # Summary for high membrane domain count proteins
    print("Summary for proteins with many membrane domains:")
    print("-" * 60)
    summary_report(
        df,
        filter_column='membrane_domain_count',
        filter_value=5,
        operator='>=',
        compare_columns=['cterm_distance', 'Prediction', 'Reduced.CC.Terms']
    )
    
    # ================================
    # SAVE FILTERED DATASETS FOR FURTHER ANALYSIS
    # ================================
    
    print("=== SAVING FILTERED DATASETS ===")
    print()
    
    # Save commonly used filtered datasets
    output_dir = Path('filtered_datasets')
    output_dir.mkdir(exist_ok=True)
    
    # Dataset 1: Close to C-terminus
    close_cterm = df[df['cterm_distance'] <= 30].copy()
    output_file1 = output_dir / 'proteins_cterm_distance_30.csv'
    close_cterm.to_csv(output_file1, index=False)
    print(f"Saved {len(close_cterm):,} proteins with cterm_distance <= 30 to {output_file1}")
    
    # Dataset 2: High membrane domain count
    high_md = df[df['membrane_domain_count'] >= 5].copy()
    output_file2 = output_dir / 'proteins_high_membrane_domains.csv'
    high_md.to_csv(output_file2, index=False)
    print(f"Saved {len(high_md):,} proteins with membrane_domain_count >= 5 to {output_file2}")
    
    # Dataset 3: Close to C-terminus AND in mass spec data
    close_and_ms = df[(df['cterm_distance'] <= 30) & (df['in_massspec'] == True)].copy()
    output_file3 = output_dir / 'proteins_cterm_30_in_massspec.csv'
    close_and_ms.to_csv(output_file3, index=False)
    print(f"Saved {len(close_and_ms):,} proteins with cterm_distance <= 30 AND in mass spec to {output_file3}")
    
    print()
    print("=== ANALYSIS COMPLETE ===")
    print()
    print("The improved workflow now allows easy filtering and comparison:")
    print("1. Filter by any criterion (cterm_distance, membrane_domain_count, etc.)")
    print("2. Compare categories within the filtered subset")
    print("3. Generate comprehensive reports and cross-tabulations")
    print("4. Save filtered datasets for further analysis")
    print()
    print("This addresses the category comparison issues mentioned in the GitHub issue.")


if __name__ == "__main__":
    main()