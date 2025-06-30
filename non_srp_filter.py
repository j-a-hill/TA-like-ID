import pandas as pd
import re

# Load data
non_biogrid_df = pd.read_csv("C:/Users/Jake/Desktop/non_biogrid.tsv", sep="\t")
signalp6_other_df = pd.read_csv("signalp6_OTHER.csv")

# Rename columns for consistency
non_biogrid_df = non_biogrid_df.rename(columns={
    'Entry': 'UniProtID',
    'Gene Names': 'GeneName'
})

# Merge on GeneName only
merged_df = pd.merge(
    non_biogrid_df,
    signalp6_other_df[['GeneName', 'Prediction', 'OTHER_Score', 'SP_Score', 'CS_Position']],
    how='left',
    on='GeneName'
)

# Filter for OTHER predictions and make a copy to avoid warnings
non_srp_df = merged_df[merged_df['Prediction'] == 'OTHER'].copy()

def count_domains(domain_str):
    if pd.isna(domain_str):
        return 0
    return len(re.findall(r'\d+\.\.\d+', domain_str))

# Count TMD and IMD domains
non_srp_df['TMD_count'] = non_srp_df['Transmembrane'].apply(count_domains)
non_srp_df['IMD_count'] = non_srp_df['Intramembrane'].apply(count_domains)

# Sum membrane domains
non_srp_df['Membrane_Domain_Count'] = non_srp_df['TMD_count'] + non_srp_df['IMD_count']

# Group by Membrane_Domain_Count and create dictionary of dfs
non_srp_membrane_dfs = {
    count: group_df.reset_index(drop=True)
    for count, group_df in non_srp_df.groupby('Membrane_Domain_Count')
}

# Example: Save one subset with a specific membrane domain count, e.g. 0
if 0 in non_srp_membrane_dfs:
    non_srp_membrane_dfs[0].to_csv("non_srp_cterm_0membrane.csv", index=False)

# Or save all groups individually
for count, df_group in non_srp_membrane_dfs.items():
    filename = f"non_srp_cterm_membrane_count_{count}.csv"
    df_group.to_csv(filename, index=False)

# Save full merged df as well
merged_df.to_csv("non_biogrid_non_SRP_df.csv", index=False)
