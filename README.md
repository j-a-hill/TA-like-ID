# TA-like-ID: Membrane Protein Analysis

A browser-based interactive platform for membrane protein analysis built with **Shiny for Python**.

## Repository Layout

```
TA-like-ID/
├── app/                                    # Shiny web application
│   ├── app.py                             # Main app (column toggles, multi-filter)
│   ├── requirements.txt                   # Python dependencies
│   ├── start_local.sh                     # Local development script
│   └── wsgi_pythonanywhere.py             # PythonAnywhere deployment config
│
├── STANDALONE_APP.py                       # Single-file app (easy to share/gist)
│
├── protein_analysis_utils.py              # Core analysis utilities
│
├── analysis/                              # Data processing scripts
│   ├── Uniprot_TMD_search.py             # Uniprot transmembrane domain search
│   ├── tmd_count.py                       # TMD count helper
│   ├── signalp_6_filter.py               # SignalP-6 output parser
│   ├── srp_filter.py                      # SRP protein filter
│   ├── non_srp_filter.py                  # Non-SRP protein filter
│   ├── improved_analysis.py              # Category comparison workflow
│   ├── demo_fixed_workflow.py            # Demo of fixed workflow
│   └── analysis_test.py                  # Exploratory analysis script
│
├── tests/                                 # Unit tests (pytest)
│   ├── test_protein_analysis_utils.py
│   ├── test_uniprot_tmd_search.py
│   └── test_filter_analysis.py
│
├── notebooks/                             # Jupyter notebooks (reference)
│   ├── membrane_protein_analysis.ipynb
│   └── MS_QTA.ipynb
│
├── raw_data/                              # Input data
│   └── membrane_protein_analysis_with_reduced_cc.csv
│
├── filtered_datasets/                     # Pre-filtered output data
│
└── docs/                                  # Documentation
    ├── DEPLOYMENT.md                      # Deployment guide (Shinylive + PythonAnywhere)
    ├── SHINYLIVE_DEPLOYMENT.md            # Shinylive-specific guide
    └── USER_GUIDE.md                      # Guide for non-technical users
```

## Quick Start

```bash
cd app
pip install -r requirements.txt
bash start_local.sh
# Visit: http://localhost:8000
```

## Deployment

**Shinylive (browser-based, no server):**
```bash
pip install shinylive
shinylive export app site
# Deploy the site/ folder to GitHub Pages, Netlify, Vercel, etc.
```

**PythonAnywhere (server-based):**
See [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md) for step-by-step instructions.

## Documentation

- [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md) — Deployment guide
- [docs/SHINYLIVE_DEPLOYMENT.md](docs/SHINYLIVE_DEPLOYMENT.md) — Shinylive-specific details
- [docs/USER_GUIDE.md](docs/USER_GUIDE.md) — Guide for non-technical users (send to colleagues)

## Features

- Multi-column filtering: protein length, C-terminal distance, membrane domain count, N-terminal domains, SRP prediction, CC terms
- Column visibility toggles for the data grid
- One-click CSV export of filtered results
- Loads data directly from GitHub (no local files needed for deployment)

## Technology Stack

- **App framework:** Shiny for Python
- **Data processing:** pandas
- **Deployment:** Shinylive (GitHub Pages / Netlify) or PythonAnywhere

## Running Tests

```bash
python -m pytest tests/ -v
```

## Notes

- The main app (`app/app.py`) loads data from GitHub at startup.
- `STANDALONE_APP.py` is a single-file version with a simpler UI; useful for sharing via GitHub Gist or shinylive.io.
- Analysis scripts in `analysis/` must be run from the project root so they can import `protein_analysis_utils`.
