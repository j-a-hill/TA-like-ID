# -*- coding: utf-8 -*-
"""
Uniprot Transmembrane Domain Search

Filter Uniprot proteins by proximity of membrane domains to C-terminus,
cross-reference with BioGRID data, and analyze domain distributions.

Uniprot query: all transmembrane and intramembrane proteins in Humans.
Downloaded as TSV with Header: Entry Entry Name Protein Name Length Transmembrane Topological domain Intramembrane
"""

import pandas as pd
import re


def extract_near_c_terminus_domains(uniprot_df: pd.DataFrame, N: int = 30) -> pd.DataFrame:
    """Return rows where at least one membrane domain end is within N residues of the C-terminus."""
    def domain_near_c_terminus(domain_str: str | float, length: int | float, N: int) -> bool:
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


def count_domains(domain_str: str | float) -> int:
    """Count the number of domain ranges (e.g. ``82..100``) in a domain annotation string."""
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', str(domain_str)))


def build_gene_set(biogrid_df: pd.DataFrame) -> set[str]:
    """Build a lowercase set of gene names from a BioGRID DataFrame."""
    columns_to_use = [
        'official_symbol_for_interactor_a', 'official_symbol_for_interactor_b',
        'synonyms/aliases_for_interactor_a', 'synonyms/aliases_for_interactor_b'
    ]
    present_cols = [c for c in columns_to_use if c in biogrid_df.columns]
    series_list = []
    for col in present_cols:
        col_data = biogrid_df[col].dropna().astype(str)
        if len(col_data) > 0:
            series_list.append(col_data.str.split('|', expand=True).stack())
    if not series_list:
        return set()
    gene_list = pd.unique(pd.concat(series_list, ignore_index=True).str.strip())
    return {g.lower() for g in gene_list if pd.notna(g) and g}


def token_match(name_field: str | float, gene_set: set[str]) -> bool:
    """Return True if any whitespace-delimited token in ``name_field`` is in ``gene_set``."""
    if pd.isna(name_field):
        return False
    return any(token in gene_set for token in str(name_field).lower().split())


def calculate_cterm_distance(row: pd.Series) -> float | None:
    """Calculate minimum C-terminal distance across all membrane domains for a protein row."""
    length = row['Length']
    if pd.isna(length):
        return None

    min_distance: float = float('inf')
    for domain_col in ('Transmembrane', 'Intramembrane'):
        if not pd.isna(row.get(domain_col, float('nan'))):
            for _, end in re.findall(r'(\d+)\.\.(\d+)', str(row[domain_col])):
                min_distance = min(min_distance, int(length) - int(end))

    return min_distance if min_distance != float('inf') else None


def analyze_filtered_proteins(df: pd.DataFrame, gene_set: set[str], cterm_threshold: int = 30) -> pd.DataFrame:
    """
    Filter proteins by C-terminal distance and compare categories within the subset.

    Adds a ``cterm_distance`` column and returns the subset within ``cterm_threshold`` residues
    of the C-terminus, printing a summary of membrane domain counts and BioGRID status.
    """
    print(f"\n=== Analysis for proteins with C-terminal domains within {cterm_threshold} residues ===")

    df = df.copy()
    df['cterm_distance'] = df.apply(calculate_cterm_distance, axis=1)
    filtered_df = df[df['cterm_distance'] <= cterm_threshold].copy()

    print(f"Found {len(filtered_df):,} proteins with membrane domains within {cterm_threshold} residues of C-terminus")
    print(f"({len(filtered_df)/len(df)*100:.1f}% of total dataset)")

    print(f"\nMembrane domain count distribution (within {cterm_threshold} residues of C-term):")
    md_counts = filtered_df['Membrane_Domain_Count'].value_counts().sort_index()
    for count, frequency in md_counts.items():
        percentage = frequency / len(filtered_df) * 100
        print(f"  {count:2d} domains: {frequency:4d} proteins ({percentage:5.1f}%)")

    print(f"\nBioGRID interaction status (within {cterm_threshold} residues of C-term):")
    biogrid_comparison = filtered_df['Gene Names'].apply(lambda x: token_match(x, gene_set)).value_counts()
    for status, count in biogrid_comparison.items():
        status_label = "In BioGRID" if status else "Not in BioGRID"
        percentage = count / len(filtered_df) * 100
        print(f"  {status_label}: {count:4d} proteins ({percentage:5.1f}%)")

    return filtered_df


if __name__ == '__main__':
    # Load Uniprot data
    uniprot_df = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2025_06_25.tsv', sep='\t', dtype='string')
    uniprot_df.columns = uniprot_df.columns.str.strip()
    uniprot_df['Length'] = pd.to_numeric(uniprot_df['Length'], errors='coerce')

    filtered_df = extract_near_c_terminus_domains(uniprot_df, N=30)

    # Load BioGRID gene names
    biogrid_df = pd.read_csv("C:/Users/Jake/Desktop/BIOGRID_SGTA.csv")
    gene_set = build_gene_set(biogrid_df)

    in_biogrid_df = filtered_df[filtered_df['Gene Names'].apply(lambda x: token_match(x, gene_set))].copy()
    non_biogrid_df = filtered_df[~filtered_df['Gene Names'].apply(lambda x: token_match(x, gene_set))].copy()

    in_biogrid_df.to_csv('in_biogrid.tsv', sep='\t', index=False)
    non_biogrid_df.to_csv('non_biogrid.tsv', sep='\t', index=False)

    # Count membrane domains for non-BioGRID proteins
    non_biogrid_df['TMD_count'] = non_biogrid_df['Transmembrane'].apply(count_domains)
    non_biogrid_df['IMD_count'] = non_biogrid_df['Intramembrane'].apply(count_domains)
    non_biogrid_df['Membrane_Domain_Count'] = non_biogrid_df['TMD_count'] + non_biogrid_df['IMD_count']

    analyzed_df = analyze_filtered_proteins(non_biogrid_df, gene_set, cterm_threshold=30)
    analyzed_df.to_csv('analyzed_cterm_proteins.tsv', sep='\t', index=False)
    print(f"\nSaved analysis results to 'analyzed_cterm_proteins.tsv'")

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
