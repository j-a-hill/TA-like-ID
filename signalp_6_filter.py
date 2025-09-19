import pandas as pd
import re
import os

# Step 1: Read and parse the SignalP-6.0 file
signalp_file = "sigp6_output/prediction_results.txt"

records = []

# Check if the input file exists
if os.path.exists(signalp_file):
    with open(signalp_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue  # Skip header/comments/blank lines

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue  # Skip malformed lines

            id_line = parts[0]
            prediction = parts[1]
            other_score = float(parts[2])
            sp_score = float(parts[3])
            cs_position = parts[4] if len(parts) > 4 else None

            # Extract UniProt ID and Gene Name
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
else:
    print(f"Warning: SignalP-6.0 input file '{signalp_file}' not found.")
    print("Creating empty dataframes as fallback.")

# Step 2: Create the main DataFrame
signalp6_df = pd.DataFrame(records)

# Step 3: Split based on Prediction
if len(signalp6_df) > 0:
    signalp6_other_df = signalp6_df[signalp6_df["Prediction"] == "OTHER"].copy()
    signalp6_sp_df = signalp6_df[signalp6_df["Prediction"] == "SP"].copy()
else:
    # Create empty dataframes with expected columns
    columns = ["RawID", "UniProtID", "GeneName", "Prediction", "OTHER_Score", "SP_Score", "CS_Position"]
    signalp6_df = pd.DataFrame(columns=columns)
    signalp6_other_df = pd.DataFrame(columns=columns)
    signalp6_sp_df = pd.DataFrame(columns=columns)

# (Optional) Save to files
signalp6_df.to_csv("signalp6_all.csv", index=False)
signalp6_other_df.to_csv("signalp6_OTHER.csv", index=False)
signalp6_sp_df.to_csv("signalp6_SP.csv", index=False)

# Print summary
print("Total predictions:", len(signalp6_df))
print("OTHER predictions:", len(signalp6_other_df))
print("SP predictions:", len(signalp6_sp_df))

if len(signalp6_df) == 0:
    print("Note: Empty CSV files created as placeholders since no input data was available.")
