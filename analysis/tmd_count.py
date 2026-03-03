# -*- coding: utf-8 -*-
"""
TMD Count Analysis

Count transmembrane and intramembrane domains from Uniprot TSV data.
"""

import pandas as pd
import re


def count_domains(domain_str: str | float) -> int:
    """Count the number of domain ranges (e.g. ``82..100``) in a domain annotation string."""
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', str(domain_str)))


if __name__ == '__main__':
    uniprot_df = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2025_06_25.tsv', sep='\t', dtype='string')

    # Count TMDs and intramembrane domains
    uniprot_df['TMD_count'] = uniprot_df['Transmembrane'].apply(count_domains)
    uniprot_df['IMD_count'] = uniprot_df['Intramembrane'].apply(count_domains)

    # Total membrane domains
    uniprot_df['Membrane_Domain_Count'] = uniprot_df['TMD_count'] + uniprot_df['IMD_count']
