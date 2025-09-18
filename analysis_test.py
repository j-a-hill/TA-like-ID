#!/usr/bin/env python3
"""
Test script to understand the current category comparison issues
and implement the fix for filtering by cterm_distance and comparing within sets
"""

import pandas as pd
import numpy as np

def main():
    # Load the data
    df = pd.read_csv('raw_data/membrane_protein_analysis_with_reduced_cc.csv')
    print("Data loaded successfully")
    print(f"Shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    # Test current category comparisons that are supposedly broken
    print("\n=== Testing current category comparisons ===")
    
    # Test 1: cterm_distance filtering
    print("\n1. cterm_distance <= 30 filter:")
    cterm_filter = df['cterm_distance'] <= 30
    print("Filter result type:", type(cterm_filter))
    print("Value counts:")
    print(cterm_filter.value_counts())
    
    # Test 2: Show what users actually want - filter first, then compare within subset
    print("\n2. Current broken workflow:")
    print("   Problem: Users can't easily filter by cterm_distance AND THEN compare MD count by localization")
    
    # Show what they actually want to do
    print("\n3. What users want to do:")
    print("   Step 1: Filter by cterm_distance <= 30")
    filtered_df = df[df['cterm_distance'] <= 30]
    print(f"   Filtered dataset size: {len(filtered_df)}")
    
    print("   Step 2: Compare MD count by some other category within that filtered set")
    print("   MD count distribution in filtered set:")
    md_counts = filtered_df['membrane_domain_count'].value_counts().sort_index()
    print(md_counts)
    
    print("   Step 3: Compare by prediction type within cterm filtered set:")
    prediction_counts = filtered_df['Prediction'].value_counts()
    print(prediction_counts)
    
    # Test what might be the "newline" issue
    print("\n=== Testing for newline issues ===")
    print("Examining text fields for newlines...")
    
    # Check for newlines in key text fields
    text_columns = ['Protein.names', 'Gene.Names', 'Subcellular.location..CC.', 'Reduced.CC.Terms']
    for col in text_columns:
        if col in df.columns:
            # Check for newlines
            has_newlines = df[col].astype(str).str.contains('\n', na=False)
            newline_count = has_newlines.sum()
            if newline_count > 0:
                print(f"Column '{col}' has {newline_count} entries with newlines")
                # Show first few examples
                examples = df[has_newlines][col].head(3).tolist()
                for i, example in enumerate(examples):
                    print(f"  Example {i+1}: {repr(example)}")
    
    # Test the filtering and comparison function that users want
    print("\n=== Implementing the desired functionality ===")
    
    def filter_and_compare(df, filter_column, filter_value, compare_column, operator='<='):
        """
        Filter by one column and then compare counts within that filtered set
        """
        if operator == '<=':
            mask = df[filter_column] <= filter_value
        elif operator == '>=':
            mask = df[filter_column] >= filter_value
        elif operator == '==':
            mask = df[filter_column] == filter_value
        else:
            raise ValueError(f"Unsupported operator: {operator}")
        
        filtered_df = df[mask]
        comparison_counts = filtered_df[compare_column].value_counts()
        
        print(f"Filter: {filter_column} {operator} {filter_value}")
        print(f"Filtered dataset size: {len(filtered_df)} (from {len(df)})")
        print(f"Comparison by {compare_column}:")
        print(comparison_counts)
        print()
        
        return filtered_df, comparison_counts
    
    # Examples of what users want to do
    print("Example 1: Filter by cterm_distance <= 30, compare MD counts")
    filter_and_compare(df, 'cterm_distance', 30, 'membrane_domain_count')
    
    print("Example 2: Filter by cterm_distance <= 50, compare by prediction type")
    filter_and_compare(df, 'cterm_distance', 50, 'Prediction')
    
    print("Example 3: Filter by membrane_domain_count == 1, compare by localization")
    if 'Reduced.CC.Terms' in df.columns:
        filter_and_compare(df, 'membrane_domain_count', 1, 'Reduced.CC.Terms')

if __name__ == "__main__":
    main()