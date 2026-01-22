# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 12:35:07 2025

@author: Jake
"""

import pandas as pd
import re

# Load data
non_biogrid_df = pd.read_csv("C:/Users/Jake/Desktop/non_biogrid.tsv", sep="\t")
signalp6_sp_df = pd.read_csv("signalp6_SP.csv")

# Rename columns for consistency
non_biogrid_df = non_biogrid_df.rename(columns={
    'Entry': 'UniProtID',
    'Gene Names': 'GeneName'
})

# Merge on GeneName only
merged_df = pd.merge(
    non_biogrid_df,
    signalp6_sp_df[['GeneName', 'Prediction', 'OTHER_Score', 'SP_Score', 'CS_Position']],
    how='left',
    on='GeneName'
)

# Filter for OTHER predictions and make a copy to avoid warnings
srp_df = merged_df[merged_df['Prediction'] == 'SP'].copy()

def count_domains(domain_str):
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', domain_str))

# Count TMD and IMD domains
srp_df['TMD_count'] = srp_df['Transmembrane'].apply(count_domains)
srp_df['IMD_count'] = srp_df['Intramembrane'].apply(count_domains)

# Sum membrane domains
srp_df['Membrane_Domain_Count'] = srp_df['TMD_count'] + srp_df['IMD_count']

# Group by Membrane_Domain_Count and create dictionary of dfs
srp_membrane_dfs = {
    count: group_df.reset_index(drop=True)
    for count, group_df in srp_df.groupby('Membrane_Domain_Count')
}

# Example: Save one subset with a specific membrane domain count, e.g. 0
if 0 in srp_membrane_dfs:
    srp_membrane_dfs[0].to_csv("srp_cterm_0membrane.csv", index=False)

# Or save all groups individually
for count, df_group in srp_membrane_dfs.items():
    filename = f"srp_cterm_membrane_count_{count}.csv"
    df_group.to_csv(filename, index=False)

# Save full merged df as well
merged_df.to_csv("non_biogrid_SRP_df.csv", index=False)

# ================================
# IMPROVED CATEGORY COMPARISON FOR SRP PROTEINS
# ================================
# Fix for GitHub issue #3 - category comparisons

def analyze_srp_by_cterm_distance(df, cterm_threshold=30):
    """
    Analyze SRP proteins by C-terminal distance to address category comparison issues.
    Users want to filter by cterm_distance and then compare MD count by localization.
    """
    print(f"\n=== SRP Protein Analysis: C-terminal distance <= {cterm_threshold} ===")
    
    # Calculate C-terminal distance for SRP proteins
    def calc_cterm_dist(row):
        length = row.get('Length', 0)
        if pd.isna(length) or length == 0:
            return float('inf')
        
        min_dist = float('inf')
        # Check transmembrane domains
        if not pd.isna(row.get('Transmembrane', '')):
            ranges = re.findall(r'(\d+)\.\.(\d+)', str(row['Transmembrane']))
            for _, end in ranges:
                dist = length - int(end)
                min_dist = min(min_dist, dist)
        
        return min_dist if min_dist != float('inf') else None
    
    # Add C-terminal distance calculation
    df['cterm_distance'] = df.apply(calc_cterm_dist, axis=1)
    
    # Filter by C-terminal distance
    close_cterm = df[df['cterm_distance'] <= cterm_threshold].copy()
    
    print(f"Found {len(close_cterm):,} SRP proteins within {cterm_threshold} residues of C-terminus")
    print(f"({len(close_cterm)/len(df)*100:.1f}% of SRP proteins)")
    
    # Compare membrane domain counts within filtered SRP proteins
    print(f"\nMembrane domain distribution (SRP, cterm <= {cterm_threshold}):")
    md_dist = close_cterm['Membrane_Domain_Count'].value_counts().sort_index()
    for count, freq in md_dist.items():
        pct = freq / len(close_cterm) * 100
        print(f"  {count:2d} domains: {freq:3d} proteins ({pct:5.1f}%)")
    
    # Save filtered SRP proteins for further analysis
    close_cterm.to_csv(f"srp_cterm_{cterm_threshold}.csv", index=False)
    print(f"\nSaved filtered SRP proteins to srp_cterm_{cterm_threshold}.csv")
    
    return close_cterm

# Apply improved analysis to SRP proteins
if 'srp_df' in locals() and len(srp_df) > 0:
    srp_analyzed = analyze_srp_by_cterm_distance(srp_df, cterm_threshold=30)
    
    # Test multiple thresholds
    print(f"\n=== SRP Protein Count by C-terminal Distance Threshold ===")
    for threshold in [10, 20, 30, 50, 100]:
        if 'cterm_distance' in srp_df.columns:
            count = len(srp_df[srp_df['cterm_distance'] <= threshold])
            pct = count / len(srp_df) * 100
            print(f"Within {threshold:3d} residues: {count:3d} proteins ({pct:5.1f}%)")

print("\n=== SRP Filter Analysis Complete ===")
print("Improved workflow allows:")
print("- Filter SRP proteins by cterm_distance")
print("- Compare MD counts within filtered subset") 
print("- Easy threshold testing")
