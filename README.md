# TA-like-ID
Aims to identify proteins that have a membrane domain close to the C-terminal. 
Takes a tsv from uniprot search with location filters for trans- or intra-membrane domains

_Uniprot_TMD_search.py_ looks for membrane domains within N residues of the C-term. It then removes known interactors from a loaded csv from biogrid.

These are then filtered by _signalp_6_filter.py_ to identify predicted signal sequences.

Results are filtered by _srp_filter.py_ and _non-srp_filter.py_ to subset them into predicted SRP/non-SRP membrane proteins.

## Improved Category Comparison Workflow

**Fixed Issues (GitHub #3):**
- ✅ Users can now easily filter by cterm_distance and compare within that subset
- ✅ Simple MD count comparison by localization for filtered proteins  
- ✅ Easy threshold testing across multiple distance values
- ✅ Streamlined workflow for complex category comparisons

### New Analysis Utilities

**protein_analysis_utils.py** - Comprehensive analysis utilities:
- `filter_and_compare()` - Filter by one criterion, compare categories within subset
- `filter_by_cterm_distance()` - Specific function for C-terminal distance filtering
- `summary_report()` - Generate comprehensive analysis reports
- `analyze_cterm_distance_effects()` - Compare effects of different thresholds

**improved_analysis.py** - Complete analysis workflow demonstrating the fixed functionality

**demo_fixed_workflow.py** - Demonstration script showing before/after comparison

### Usage Examples

```python
from protein_analysis_utils import filter_and_compare, filter_by_cterm_distance

# Main issue fix: Filter by cterm_distance, compare MD counts
filtered_df, md_counts = filter_by_cterm_distance(df, max_distance=30, compare_by='membrane_domain_count')

# General filtering and comparison
filtered_df, counts = filter_and_compare(df, 'cterm_distance', 50, 'Prediction', operator='<=')

# Easy threshold testing
for threshold in [10, 20, 30, 50]:
    subset, md_dist = filter_by_cterm_distance(df, threshold, 'membrane_domain_count')
    print(f"<= {threshold}: {len(subset)} proteins")
```

### Files Generated

The improved workflow creates filtered datasets in `filtered_datasets/`:
- `proteins_cterm_distance_30.csv` - Proteins within 30 residues of C-terminus
- `proteins_high_membrane_domains.csv` - Proteins with >= 5 membrane domains
- `proteins_cterm_30_in_massspec.csv` - C-terminal proteins detected in mass spec

This addresses the category comparison issues and provides a streamlined workflow for membrane protein analysis.
