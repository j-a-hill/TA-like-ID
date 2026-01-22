#!/usr/bin/env python3
"""
Protein Analysis Utilities

This module provides improved functionality for filtering and comparing
membrane protein data, specifically addressing the workflow issues where
users want to filter by cterm_distance and then compare categories within
the filtered subset.

Functions:
- filter_and_compare: Main function for filtering and comparing categories
- filter_by_cterm_distance: Specific function for C-terminal distance filtering  
- compare_categories: Compare categories within a dataset
- summary_report: Generate summary reports for filtered data
"""

import pandas as pd
import numpy as np
from typing import Union, Dict, List, Tuple, Optional


def filter_and_compare(
    df: pd.DataFrame,
    filter_column: str,
    filter_value: Union[int, float, str],
    compare_column: str,
    operator: str = '<=',
    sort_by_count: bool = True,
    top_n: Optional[int] = None
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Filter a DataFrame by one column and then compare counts within that filtered set.
    
    This addresses the main issue where users want to filter by cterm_distance 
    and then compare MD count by localization within that subset.
    
    Args:
        df: Input DataFrame
        filter_column: Column to filter on
        filter_value: Value to filter by
        compare_column: Column to compare categories within filtered data
        operator: Comparison operator ('<=', '>=', '==', '!=', '<', '>')
        sort_by_count: Whether to sort results by count (descending)
        top_n: Show only top N categories (None for all)
    
    Returns:
        Tuple of (filtered_dataframe, comparison_counts)
    """
    # Apply filter
    if operator == '<=':
        mask = df[filter_column] <= filter_value
    elif operator == '>=':
        mask = df[filter_column] >= filter_value
    elif operator == '==':
        mask = df[filter_column] == filter_value
    elif operator == '!=':
        mask = df[filter_column] != filter_value
    elif operator == '<':
        mask = df[filter_column] < filter_value
    elif operator == '>':
        mask = df[filter_column] > filter_value
    else:
        raise ValueError(f"Unsupported operator: {operator}")
    
    filtered_df = df[mask].copy()
    
    # Get comparison counts
    comparison_counts = filtered_df[compare_column].value_counts()
    
    if sort_by_count:
        comparison_counts = comparison_counts.sort_values(ascending=False)
    
    if top_n is not None:
        comparison_counts = comparison_counts.head(top_n)
    
    return filtered_df, comparison_counts


def filter_by_cterm_distance(
    df: pd.DataFrame,
    max_distance: int = 30,
    compare_by: str = 'membrane_domain_count'
) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Convenience function for the common use case of filtering by C-terminal distance
    and comparing within that subset.
    
    Args:
        df: Input DataFrame
        max_distance: Maximum C-terminal distance (default 30)
        compare_by: Column to compare categories by
    
    Returns:
        Tuple of (filtered_dataframe, comparison_counts)
    """
    return filter_and_compare(
        df=df,
        filter_column='cterm_distance',
        filter_value=max_distance,
        compare_column=compare_by,
        operator='<='
    )


def compare_categories(
    df: pd.DataFrame,
    category_column: str,
    group_by: Optional[str] = None,
    normalize: bool = False
) -> Union[pd.Series, pd.DataFrame]:
    """
    Compare categories within a dataset, optionally grouped by another column.
    
    Args:
        df: Input DataFrame
        category_column: Column with categories to compare
        group_by: Optional column to group by before comparing
        normalize: Whether to normalize counts to proportions
    
    Returns:
        Series or DataFrame with category counts/proportions
    """
    if group_by is None:
        result = df[category_column].value_counts(normalize=normalize)
    else:
        if normalize:
            result = df.groupby(group_by)[category_column].value_counts(normalize=True).unstack(fill_value=0)
        else:
            result = df.groupby(group_by)[category_column].value_counts().unstack(fill_value=0)
    
    return result


def summary_report(
    df: pd.DataFrame,
    filter_column: str = 'cterm_distance',
    filter_value: Union[int, float] = 30,
    operator: str = '<=',
    compare_columns: List[str] = None
) -> None:
    """
    Generate a comprehensive summary report for filtered data.
    
    Args:
        df: Input DataFrame
        filter_column: Column to filter on
        filter_value: Value to filter by
        operator: Comparison operator
        compare_columns: List of columns to compare within filtered data
    """
    if compare_columns is None:
        compare_columns = ['membrane_domain_count', 'Prediction', 'in_biogrid', 'in_massspec']
    
    # Filter the data
    filtered_df, _ = filter_and_compare(df, filter_column, filter_value, compare_columns[0], operator)
    
    print(f"=== Summary Report ===")
    print(f"Filter: {filter_column} {operator} {filter_value}")
    print(f"Original dataset size: {len(df):,}")
    print(f"Filtered dataset size: {len(filtered_df):,}")
    print(f"Percentage retained: {len(filtered_df)/len(df)*100:.1f}%")
    print()
    
    # Compare categories within filtered data
    for col in compare_columns:
        if col in filtered_df.columns:
            print(f"--- {col.replace('_', ' ').title()} Distribution ---")
            counts = filtered_df[col].value_counts()
            for category, count in counts.head(10).items():
                percentage = count / len(filtered_df) * 100
                print(f"  {category}: {count:,} ({percentage:.1f}%)")
            if len(counts) > 10:
                print(f"  ... and {len(counts) - 10} more categories")
            print()


def analyze_cterm_distance_effects(
    df: pd.DataFrame,
    distance_thresholds: List[int] = None,
    compare_column: str = 'membrane_domain_count'
) -> pd.DataFrame:
    """
    Analyze how different C-terminal distance thresholds affect category distributions.
    
    Args:
        df: Input DataFrame
        distance_thresholds: List of distance thresholds to test
        compare_column: Column to compare across thresholds
    
    Returns:
        DataFrame showing category distributions for each threshold
    """
    if distance_thresholds is None:
        distance_thresholds = [10, 20, 30, 50, 100]
    
    results = {}
    
    for threshold in distance_thresholds:
        filtered_df, counts = filter_by_cterm_distance(df, threshold, compare_column)
        results[f"<= {threshold}"] = counts
    
    # Combine into DataFrame
    result_df = pd.DataFrame(results).fillna(0).astype(int)
    
    return result_df


def cross_tabulate_categories(
    df: pd.DataFrame,
    column1: str,
    column2: str,
    filter_column: Optional[str] = None,
    filter_value: Optional[Union[int, float]] = None,
    operator: str = '<='
) -> pd.DataFrame:
    """
    Create a cross-tabulation of two categorical columns, optionally after filtering.
    
    Args:
        df: Input DataFrame
        column1: First categorical column
        column2: Second categorical column
        filter_column: Optional column to filter on first
        filter_value: Optional value to filter by
        operator: Comparison operator for filtering
    
    Returns:
        Cross-tabulation DataFrame
    """
    # Apply filter if specified
    if filter_column is not None and filter_value is not None:
        filtered_df, _ = filter_and_compare(df, filter_column, filter_value, column1, operator)
    else:
        filtered_df = df
    
    # Create cross-tabulation
    crosstab = pd.crosstab(filtered_df[column1], filtered_df[column2], margins=True)
    
    return crosstab


# Convenience functions for common analysis patterns
def quick_cterm_analysis(df: pd.DataFrame, distance: int = 30) -> None:
    """Quick analysis of proteins within C-terminal distance threshold."""
    print(f"Quick Analysis: Proteins with cterm_distance <= {distance}")
    print("=" * 60)
    
    filtered_df, md_counts = filter_by_cterm_distance(df, distance, 'membrane_domain_count')
    
    print(f"Found {len(filtered_df):,} proteins within {distance} residues of C-terminus")
    print(f"({len(filtered_df)/len(df)*100:.1f}% of total dataset)")
    print()
    
    print("Top membrane domain counts:")
    for count, frequency in md_counts.head(10).items():
        print(f"  {count} domains: {frequency:,} proteins")
    print()
    
    # Prediction distribution
    pred_counts = filtered_df['Prediction'].value_counts()
    print("Prediction distribution:")
    for pred, count in pred_counts.items():
        print(f"  {pred}: {count:,} proteins ({count/len(filtered_df)*100:.1f}%)")


def quick_localization_analysis(
    df: pd.DataFrame,
    distance: int = 30,
    localization_column: str = 'Reduced.CC.Terms'
) -> None:
    """Quick analysis of subcellular localization for filtered proteins."""
    print(f"Localization Analysis: Proteins with cterm_distance <= {distance}")
    print("=" * 60)
    
    filtered_df, loc_counts = filter_by_cterm_distance(df, distance, localization_column)
    
    print(f"Found {len(filtered_df):,} proteins within {distance} residues of C-terminus")
    print()
    
    print("Top localizations:")
    for loc, count in loc_counts.head(15).items():
        print(f"  {count:,}: {loc}")


if __name__ == "__main__":
    # Example usage
    print("Protein Analysis Utils - Example Usage")
    print("=" * 50)
    
    # Load sample data if available
    try:
        df = pd.read_csv('raw_data/membrane_protein_analysis_with_reduced_cc.csv')
        print(f"Loaded data with {len(df):,} proteins")
        
        # Quick analysis
        quick_cterm_analysis(df, 30)
        print()
        quick_localization_analysis(df, 30)
        
    except FileNotFoundError:
        print("Sample data not found. This module provides utilities for protein analysis.")
        print("Use the functions in your own analysis scripts.")