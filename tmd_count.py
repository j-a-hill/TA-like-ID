# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 12:45:46 2025

@author: Jake"""


import pandas as pd
import re

uniprot_df = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2025_06_25.tsv', sep='\t', dtype='string')

def count_domains(domain_str):
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', domain_str))

# Step 1: Count TMDs and intramembrane domains
uniprot_df['TMD_count'] = uniprot_df['Transmembrane'].apply(count_domains)
uniprot_df['IMD_count'] = uniprot_df['Intramembrane'].apply(count_domains)

# Step 2: Total membrane domains
uniprot_df['Membrane_Domain_Count'] = uniprot_df['TMD_count'] + uniprot_df['IMD_count']