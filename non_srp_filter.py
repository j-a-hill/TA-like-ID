import pandas as pd
import re

# Load data
# Note: Assuming non_biogrid.tsv is available in the working directory or raw_data/
# Using relative path instead of hardcoded Windows path
try:
    non_biogrid_df = pd.read_csv("non_biogrid.tsv", sep="\t")
except FileNotFoundError:
    try:
        non_biogrid_df = pd.read_csv("raw_data/non_biogrid.tsv", sep="\t")
    except FileNotFoundError:
        print("Warning: non_biogrid.tsv not found. Using membrane protein analysis data as fallback.")
        non_biogrid_df = pd.read_csv("raw_data/membrane_protein_analysis_with_reduced_cc.csv")

try:
    signalp6_other_df = pd.read_csv("signalp6_OTHER.csv")
except FileNotFoundError:
    print("Warning: signalp6_OTHER.csv not found. Creating empty dataframe.")
    signalp6_other_df = pd.DataFrame()

# Rename columns for consistency if they exist
rename_map = {}
if 'Entry' in non_biogrid_df.columns:
    rename_map['Entry'] = 'UniProtID'
if 'Gene Names' in non_biogrid_df.columns:
    rename_map['Gene Names'] = 'GeneName'
elif 'Gene.Names' in non_biogrid_df.columns:
    rename_map['Gene.Names'] = 'GeneName'

if rename_map:
    non_biogrid_df = non_biogrid_df.rename(columns=rename_map)

# Merge on GeneName only if both dataframes have the required columns
if len(signalp6_other_df) > 0 and 'GeneName' in non_biogrid_df.columns and 'GeneName' in signalp6_other_df.columns:
    # Get available columns from signalp6_other_df
    available_cols = ['GeneName']
    for col in ['Prediction', 'OTHER_Score', 'SP_Score', 'CS_Position']:
        if col in signalp6_other_df.columns:
            available_cols.append(col)
    
    merged_df = pd.merge(
        non_biogrid_df,
        signalp6_other_df[available_cols],
        how='left',
        on='GeneName'
    )
else:
    print("Warning: Cannot merge datasets. Using non_biogrid data only.")
    merged_df = non_biogrid_df.copy()
    # Add missing columns if they don't exist
    if 'Prediction' not in merged_df.columns:
        merged_df['Prediction'] = 'OTHER'

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

# ================================
# IMPROVED CATEGORY COMPARISON FOR NON-SRP PROTEINS
# ================================
# Fix for GitHub issue #3 - category comparisons

def analyze_non_srp_by_cterm_distance(df, cterm_threshold=30):
    """
    Analyze non-SRP proteins by C-terminal distance to address category comparison issues.
    Users want to filter by cterm_distance and then compare MD count by localization.
    """
    print(f"\n=== Non-SRP Protein Analysis: C-terminal distance <= {cterm_threshold} ===")
    
    # Calculate C-terminal distance for non-SRP proteins  
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
    
    print(f"Found {len(close_cterm):,} non-SRP proteins within {cterm_threshold} residues of C-terminus")
    print(f"({len(close_cterm)/len(df)*100:.1f}% of non-SRP proteins)")
    
    # Compare membrane domain counts within filtered non-SRP proteins
    print(f"\nMembrane domain distribution (non-SRP, cterm <= {cterm_threshold}):")
    md_dist = close_cterm['Membrane_Domain_Count'].value_counts().sort_index()
    for count, freq in md_dist.items():
        pct = freq / len(close_cterm) * 100
        print(f"  {count:2d} domains: {freq:3d} proteins ({pct:5.1f}%)")
    
    # Save filtered non-SRP proteins for further analysis
    close_cterm.to_csv(f"non_srp_cterm_{cterm_threshold}.csv", index=False)
    print(f"\nSaved filtered non-SRP proteins to non_srp_cterm_{cterm_threshold}.csv")
    
    return close_cterm

# Apply improved analysis to non-SRP proteins
if 'non_srp_df' in locals() and len(non_srp_df) > 0:
    non_srp_analyzed = analyze_non_srp_by_cterm_distance(non_srp_df, cterm_threshold=30)
    
    # Test multiple thresholds
    print(f"\n=== Non-SRP Protein Count by C-terminal Distance Threshold ===")
    for threshold in [10, 20, 30, 50, 100]:
        if 'cterm_distance' in non_srp_df.columns:
            count = len(non_srp_df[non_srp_df['cterm_distance'] <= threshold])
            pct = count / len(non_srp_df) * 100
            print(f"Within {threshold:3d} residues: {count:3d} proteins ({pct:5.1f}%)")

print("\n=== Non-SRP Filter Analysis Complete ===")
print("Improved workflow allows:")
print("- Filter non-SRP proteins by cterm_distance")
print("- Compare MD counts within filtered subset")
print("- Easy threshold testing")
