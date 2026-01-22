# Web App User Guide

## How to Use the Membrane Protein Analysis Tool

### Getting Started

1. **Open the app** in your web browser (URL provided by your administrator)
2. You should see a clean interface with filters on the left and results on the right

### Filtering Your Data

#### Step 1: Set Distance Threshold
- **Filter by C-terminal Distance**: Select an operator and enter a value
  - `≤` (less than or equal): Proteins with C-terminal distance ≤ your value
  - `≥` (greater than or equal): Proteins with distance ≥ your value
  - `==` (exactly): Proteins with exact distance match
  - `<` (less than): Proteins with distance < your value
  - `>` (greater than): Proteins with distance > your value

**Example**: Setting `≤ 30` will show only proteins with C-terminal domain within 30 residues.

#### Step 2: Choose What to Compare
- **Compare by Category**: Select what you want to analyze in your filtered set:
  - **Membrane Domain Count**: Shows how many proteins have 1, 2, 3, etc. membrane domains
  - **Prediction**: Shows distribution of signal peptide predictions
  - **Localization**: Shows membrane protein localization patterns

#### Step 3: Apply Filters
- Click the blue **"Filter Data"** button to apply your criteria
- The app will show:
  - Summary statistics (total, filtered count, reduction %)
  - Bar chart of category distribution
  - Table of matching proteins

### Viewing Results

#### Statistics Box
Shows:
- **Total Proteins**: Your complete dataset
- **Filtered Results**: Proteins matching your criteria
- **Reduction %**: How much the data was filtered

#### Category Distribution Chart
A bar chart showing the breakdown of your selected comparison category across filtered proteins.

#### Data Table
Scroll through your filtered results. The table shows the first 100 rows (if more exist, you'll see a note).

### Downloading Results

1. Click the green **"📥 Download CSV"** button
2. Your filtered data will download as `filtered_proteins.csv`
3. Open in Excel, Google Sheets, or other spreadsheet software
4. Share with colleagues or perform additional analysis

### Common Workflows

#### Find proteins with C-terminal domains within 30 residues
1. Set distance to `≤ 30`
2. Compare by "Membrane Domain Count"
3. See how many have single vs multiple domains
4. Download results

#### Explore different distance thresholds
1. Try `≤ 20` → note the number of proteins
2. Reset filter
3. Try `≤ 50` → note the difference
4. Download the version you prefer

#### Identify specific localization patterns
1. Set your distance threshold
2. Compare by "Localization"
3. See where filtered proteins are typically located
4. Download for downstream analysis

### Tips & Tricks

✓ **Reset anytime**: Click "Reset" to clear filters and start over
✓ **Multiple criteria**: Each filter application builds on previous results
✓ **CSV reusable**: Downloaded CSVs can be used in other tools and analyses
✓ **Zoom charts**: Hover over the bar chart and use the Plotly toolbar to zoom/pan
✓ **Export charts**: Use the camera icon in the chart toolbar to save as PNG

### What Each Column Means

- **cterm_distance**: Distance (in residues) from protein C-terminus to nearest membrane domain
- **membrane_domain_count**: Number of transmembrane domains in the protein
- **Prediction**: Signal peptide prediction (SRP, non-SRP, etc.)
- **Localization**: Where the protein is localized in the cell

### Questions or Issues?

Contact your administrator with the URL and your filter criteria if you need help.
