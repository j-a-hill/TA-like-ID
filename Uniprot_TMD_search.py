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
biogrid_df = pd.read_csv("C:/Users/Jake/Desktop/BIOGRID_SGTA.csv")

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

