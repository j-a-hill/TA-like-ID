"""
SignalP-6.0 Output Parser

Parse SignalP-6.0 prediction results, extract UniProt IDs and gene names,
and split predictions into SP and OTHER groups.
"""

import pandas as pd
import re


def parse_signalp6(signalp_file: str) -> pd.DataFrame:
    """Parse a SignalP-6.0 prediction results file and return a DataFrame."""
    records = []
    with open(signalp_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            id_line = parts[0]
            prediction = parts[1]
            other_score = float(parts[2])
            sp_score = float(parts[3])
            cs_position = parts[4] if len(parts) > 4 else None

            match = re.search(r"sp\|(\w+)\|[^\s]+.*GN=([\w\-]+)", id_line)
            if match:
                uniprot_id = match.group(1)
                gene_name = match.group(2)
            else:
                uniprot_id = None
                gene_name = None

            records.append({
                "RawID": id_line,
                "UniProtID": uniprot_id,
                "GeneName": gene_name,
                "Prediction": prediction,
                "OTHER_Score": other_score,
                "SP_Score": sp_score,
                "CS_Position": cs_position
            })

    return pd.DataFrame(records)


if __name__ == '__main__':
    signalp_file = "sigp6_output/prediction_results.txt"
    signalp6_df = parse_signalp6(signalp_file)

    signalp6_other_df = signalp6_df[signalp6_df["Prediction"] == "OTHER"].copy()
    signalp6_sp_df = signalp6_df[signalp6_df["Prediction"] == "SP"].copy()

    signalp6_df.to_csv("signalp6_all.csv", index=False)
    signalp6_other_df.to_csv("signalp6_OTHER.csv", index=False)
    signalp6_sp_df.to_csv("signalp6_SP.csv", index=False)

    print("Total predictions:", len(signalp6_df))
    print("OTHER predictions:", len(signalp6_other_df))
    print("SP predictions:", len(signalp6_sp_df))

