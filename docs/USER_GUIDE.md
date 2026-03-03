# Web App User Guide

## How to Use the Membrane Protein Analysis Tool

### Getting Started

Open the app in your web browser (URL provided by your administrator). You will see a sidebar of filters on the left and an interactive data grid on the right.

---

### Column Visibility

At the top of the page is a **Column Visibility** panel. Tick or untick the checkboxes to show or hide individual columns in the data grid below.

---

### Filtering Your Data

The sidebar contains several filter groups. All active filters are applied simultaneously (AND logic).

#### Protein Length
Use the slider to set a minimum and maximum protein length.

#### C-terminal Distance
Select an operator (`<`, `>`, or `=`) and enter a threshold value.

**Example:** Operator `<`, value `30` → proteins with a C-terminal domain within 30 residues of the C-terminus.

#### Penultimate Distance
Same operator/value controls as C-terminal Distance, applied to the penultimate membrane domain distance column.

#### Membrane Domain Count
Operator/value filter on the number of transmembrane domains.

#### N-terminal Membrane Domains
Operator/value filter on N-terminal membrane domain count.

#### SRP Prediction
Select `ALL`, `SP` (signal peptide), or `OTHER` from the dropdown.

#### CC Terms
Type a search term to filter by subcellular localisation annotation (e.g. `ER`, `secretory`). The search is case-insensitive.

---

### Downloading Results

Click **📥 Download CSV** in the sidebar. The currently filtered data will download as `filtered_proteins.csv`, which opens directly in Excel or Google Sheets.

---

### Common Workflows

#### Find proteins with C-terminal domains within 30 residues
1. Set C-terminal Distance → Operator `<`, Value `30`
2. The grid updates immediately
3. Click **Download CSV** to save results

#### Explore different distance thresholds
1. Set C-terminal Distance to `< 20` — note the row count
2. Change to `< 50` — compare the difference
3. Download whichever threshold you prefer

#### Focus on SRP-independent proteins
1. Set SRP Prediction → `OTHER`
2. Combine with a C-terminal distance or domain count filter
3. Download filtered results

---

### What Each Column Means

| Column | Meaning |
|--------|---------|
| `cterm_distance` | Residues from the protein C-terminus to the nearest membrane domain end |
| `penultimate_distance` | Same metric for the second-closest domain |
| `membrane_domain_count` | Total number of transmembrane domains |
| `N_term_md` | Number of membrane domains in the N-terminal region |
| `Prediction` | SignalP-6 prediction: `SP` (signal peptide) or `OTHER` |
| `Reduced.CC.Terms` | Simplified subcellular localisation terms |

---

### Questions or Issues?

Contact your administrator with the URL and your filter criteria if you need help.

