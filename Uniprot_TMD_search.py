# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 11:16:21 2025

@author: Jake
"""

# Uniprot query = all transmembrane and intramembrane proteins in Humans. Returned 205205 total entries. 
# Downloaded as TSV with Header = Entry Entry Name Protein Name Length Transmembrane Topological domain Intramembrane

import pandas as pd
import re

# Load data
uniprot_df = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2025_06_25.tsv', sep='\t', dtype='string')
uniprot_df.columns = uniprot_df.columns.str.strip()  # Remove potential whitespace
uniprot_df['Length'] = pd.to_numeric(uniprot_df['Length'], errors='coerce')


def extract_near_c_terminus_domains(uniprot_df, N=30):
    def domain_near_c_terminus(domain_str, length, N):
        if pd.isna(domain_str) or pd.isna(length):
            return False
        # Find all ranges like 82..100
        ranges = re.findall(r'(\d+)\.\.(\d+)', domain_str)
        for _, end in ranges:
            if length - int(end) <= N:
                return True
        return False

    mask = (
        uniprot_df.apply(lambda row: domain_near_c_terminus(row['Transmembrane'], row['Length'], N), axis=1) |
        uniprot_df.apply(lambda row: domain_near_c_terminus(row['Intramembrane'], row['Length'], N), axis=1) 
        )
    
    return uniprot_df[mask]

# Run function
filtered_df = extract_near_c_terminus_domains(uniprot_df, N=30)

# Save result if desired
# filtered_df.to_csv('near_c_terminal_domains.tsv', sep='\t', index=False)

# Load the second gene name file
# Using relative path instead of hardcoded Windows path
try:
    biogrid_df = pd.read_csv("BIOGRID_SGTA.csv")
except FileNotFoundError:
    try:
        biogrid_df = pd.read_csv("raw_data/BIOGRID_SGTA.csv")
    except FileNotFoundError:
        print("Warning: BIOGRID_SGTA.csv not found. Using available SGTA data as fallback.")
        # Try using one of the available SGTA files
        try:
            biogrid_df = pd.read_csv("raw_data/SGTA and OP91 MS.csv")
        except FileNotFoundError:
            biogrid_df = pd.read_csv("raw_data/SGTA KD Prot Inhibition(ORIGINAL).csv")

# Columns that may contain gene names
columns_to_use = ['official_symbol_for_interactor_a', 'official_symbol_for_interactor_b', 'synonyms/aliases_for_interactor_a', 'synonyms/aliases_for_interactor_b' ]

# Combine all gene name entries into a single flat list
gene_list = pd.unique(
    pd.concat([biogrid_df[col].dropna().str.split('|', expand=True).stack()
               for col in columns_to_use], ignore_index=True)
    .str.strip()
)

# Lowercase and clean gene list
gene_set = set(g.lower() for g in gene_list if pd.notna(g))

# Define match function (whole-token match)
def token_match(name_field):
    if pd.isna(name_field):
        return False
    tokens = name_field.lower().split()
    return any(token in gene_set for token in tokens)

# Apply match
in_biogrid_df = filtered_df[filtered_df['Gene Names'].apply(token_match)].copy()
non_biogrid_df = filtered_df[~filtered_df['Gene Names'].apply(token_match)].copy()

# Optional: save to files
in_biogrid_df.to_csv('in_biogrid.tsv', sep='\t', index=False)
non_biogrid_df.to_csv('non_biogrid.tsv', sep='\t', index=False)

def count_domains(domain_str):
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', domain_str))

# Step 1: Count TMDs and intramembrane domains
non_biogrid_df['TMD_count'] = non_biogrid_df['Transmembrane'].apply(count_domains)
non_biogrid_df['IMD_count'] = non_biogrid_df['Intramembrane'].apply(count_domains)

# Step 2: Total membrane domains
non_biogrid_df['Membrane_Domain_Count'] = non_biogrid_df['TMD_count'] + non_biogrid_df['IMD_count']

# Step 3: Group by total count and split
membrane_dfs = {
    count: group_df.reset_index(drop=True)
    for count, group_df in non_biogrid_df.groupby('Membrane_Domain_Count')
}
# Optional: save each group to file
#for tmd_count, df in tmd_dfs.items():
    #df.to_csv(f'non_overlap_TMDcount_{tmd_count}.tsv', sep='\t', index=False)

# ================================
# IMPROVED CATEGORY COMPARISON WORKFLOW
# ================================
# Fix for the category comparison issues mentioned in GitHub issue #3

def analyze_filtered_proteins(df, cterm_threshold=30):
    """
    Improved function to filter by cterm_distance and compare categories within subset.
    This addresses the workflow issue where users want to filter by cterm_distance 
    and then compare MD count by localization.
    """
    print(f"\n=== Analysis for proteins with C-terminal domains within {cterm_threshold} residues ===")
    
    # Calculate C-terminal distance for each protein
    def calculate_cterm_distance(row):
        """Calculate distance from C-terminus for the closest membrane domain."""
        length = row['Length']
        if pd.isna(length):
            return float('inf')
        
        min_distance = float('inf')
        
        # Check transmembrane domains
        if not pd.isna(row['Transmembrane']):
            ranges = re.findall(r'(\d+)\.\.(\d+)', str(row['Transmembrane']))
            for _, end in ranges:
                distance = length - int(end)
                min_distance = min(min_distance, distance)
        
        # Check intramembrane domains  
        if not pd.isna(row['Intramembrane']):
            ranges = re.findall(r'(\d+)\.\.(\d+)', str(row['Intramembrane']))
            for _, end in ranges:
                distance = length - int(end)
                min_distance = min(min_distance, distance)
        
        return min_distance if min_distance != float('inf') else None
    
    # Add C-terminal distance column
    df['cterm_distance'] = df.apply(calculate_cterm_distance, axis=1)
    
    # Filter by C-terminal distance
    filtered_df = df[df['cterm_distance'] <= cterm_threshold].copy()
    
    print(f"Found {len(filtered_df):,} proteins with membrane domains within {cterm_threshold} residues of C-terminus")
    print(f"({len(filtered_df)/len(df)*100:.1f}% of total dataset)")
    
    # Compare membrane domain counts within filtered set
    print(f"\nMembrane domain count distribution (within {cterm_threshold} residues of C-term):")
    md_counts = filtered_df['Membrane_Domain_Count'].value_counts().sort_index()
    for count, frequency in md_counts.items():
        percentage = frequency / len(filtered_df) * 100
        print(f"  {count:2d} domains: {frequency:4d} proteins ({percentage:5.1f}%)")
    
    # Compare BioGRID status within filtered set
    print(f"\nBioGRID interaction status (within {cterm_threshold} residues of C-term):")
    biogrid_comparison = filtered_df.apply(token_match, axis=1).value_counts()
    for status, count in biogrid_comparison.items():
        status_label = "In BioGRID" if status else "Not in BioGRID"
        percentage = count / len(filtered_df) * 100
        print(f"  {status_label}: {count:4d} proteins ({percentage:5.1f}%)")
    
    return filtered_df

# Apply the improved analysis
if 'filtered_df' in locals():
    analyzed_df = analyze_filtered_proteins(filtered_df, cterm_threshold=30)
    
    # Save the analyzed results with C-terminal distance information
    analyzed_df.to_csv('analyzed_cterm_proteins.tsv', sep='\t', index=False)
    print(f"\nSaved analysis results to 'analyzed_cterm_proteins.tsv'")
    
    # Also test with different thresholds
    print(f"\n=== Comparison across different C-terminal distance thresholds ===")
    for threshold in [10, 20, 30, 50]:
        subset = analyzed_df[analyzed_df['cterm_distance'] <= threshold]
        print(f"Within {threshold:2d} residues: {len(subset):4d} proteins ({len(subset)/len(analyzed_df)*100:5.1f}%)")

print(f"\n=== Category comparison workflow improved ===")
print("Users can now:")
print("1. Filter by cterm_distance threshold") 
print("2. Compare MD count by localization within that filtered set")
print("3. Compare other categories (BioGRID status, etc.) within filtered set")
print("4. Test multiple thresholds easily")

