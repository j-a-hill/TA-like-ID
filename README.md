# TA-like-ID
Aims to identify proteins that have a membrane domain close to the C-termninal. 
Takes a tsv form uniprot search with location filters for trans- or intra-membrane domains

_Uniprot_TMD_search.py_ looks for membrane domains within N residues of the C-term. It then removes known interactors from a loaded csv from biogrid.

These are then filtered by _signalp_6_filter.py_ to identify predicted signal sequences.

Results are filtered by _srp_filter.py_ and _non-srp_filter.py_ to subset them into predicted SRP/non-SRP membrane proteins.
