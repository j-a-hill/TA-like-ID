"""
SRP Protein Filter Analysis

Filter SRP (SP prediction) proteins from BioGRID-negative set,
classify by membrane domain count, and analyze C-terminal domain proximity.
"""

import pandas as pd
import re

from protein_analysis_utils import calc_min_cterm_distance


def count_domains(domain_str: str | float) -> int:
    """Count the number of domain ranges (e.g. ``82..100``) in a domain annotation string."""
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', str(domain_str)))


def analyze_srp_by_cterm_distance(df: pd.DataFrame, cterm_threshold: int = 30) -> pd.DataFrame:
    """
    Analyze SRP proteins by C-terminal distance to address category comparison issues.
    Users want to filter by cterm_distance and then compare MD count by localization.
    """
    print(f"\n=== SRP Protein Analysis: C-terminal distance <= {cterm_threshold} ===")

    df = df.copy()
    df['cterm_distance'] = df.apply(calc_min_cterm_distance, axis=1)

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


if __name__ == '__main__':
    # Load data
    non_biogrid_df = pd.read_csv("C:/Users/Jake/Desktop/non_biogrid.tsv", sep="\t")
    signalp6_sp_df = pd.read_csv("signalp6_SP.csv")

    # Rename columns for consistency
    non_biogrid_df = non_biogrid_df.rename(columns={'Entry': 'UniProtID', 'Gene Names': 'GeneName'})

    merged_df = pd.merge(
        non_biogrid_df,
        signalp6_sp_df[['GeneName', 'Prediction', 'OTHER_Score', 'SP_Score', 'CS_Position']],
        how='left', on='GeneName'
    )

    srp_df = merged_df[merged_df['Prediction'] == 'SP'].copy()

    srp_df['TMD_count'] = srp_df['Transmembrane'].apply(count_domains)
    srp_df['IMD_count'] = srp_df['Intramembrane'].apply(count_domains)
    srp_df['Membrane_Domain_Count'] = srp_df['TMD_count'] + srp_df['IMD_count']

    srp_membrane_dfs = {
        count: group_df.reset_index(drop=True)
        for count, group_df in srp_df.groupby('Membrane_Domain_Count')
    }

    for count, df_group in srp_membrane_dfs.items():
        df_group.to_csv(f"srp_cterm_membrane_count_{count}.csv", index=False)

    merged_df.to_csv("non_biogrid_SRP_df.csv", index=False)

    srp_analyzed = analyze_srp_by_cterm_distance(srp_df, cterm_threshold=30)

    print(f"\n=== SRP Protein Count by C-terminal Distance Threshold ===")
    for threshold in [10, 20, 30, 50, 100]:
        count = len(srp_analyzed[srp_analyzed['cterm_distance'] <= threshold])
        pct = count / len(srp_analyzed) * 100
        print(f"Within {threshold:3d} residues: {count:3d} proteins ({pct:5.1f}%)")

    print("\n=== SRP Filter Analysis Complete ===")

