# TA-like-ID: Shiny Web Application

A browser-based interactive platform for membrane protein analysis using **Shiny for Python** with **Shinylive** deployment.

## 🌟 New: Shinylive Deployment

This app now uses **Shiny for Python** and can be deployed via **Shinylive**, which means:
- ✅ Runs entirely in the browser (no server needed!)
- ✅ Free hosting on GitHub Pages, Netlify, etc.
- ✅ Data stays private (never leaves the browser)
- ✅ Share via simple URL
- ✅ No backend maintenance required

```
TA-like-ID/
├── app/                                # 🚀 Shiny application
│   ├── app.py                         # Main Shiny app
│   ├── protein_analysis_utils.py      # Analysis utilities
│   ├── requirements.txt               # Python dependencies (Shiny)
│   └── start_local.sh                 # Local development script
│
├── docs/                              # 📚 Documentation
│   ├── SHINYLIVE_DEPLOYMENT.md        # Shinylive deployment guide
│   ├── USER_GUIDE.md                  # Guide for non-technical users
│   ├── FILE_INDEX.md                  # Guide to all files
│   └── (other docs...)
│
├── analysis/                          # 📊 Original analysis scripts
│   ├── protein_analysis_utils.py      # Shared analysis utilities
│   ├── improved_analysis.py           # Improved analysis workflow
│   └── (other analysis scripts...)
│
├── notebooks/                         # 📓 Jupyter notebooks
│   ├── membrane_protein_analysis.ipynb
│   └── MS_QTA.ipynb
│
├── raw_data/                          # 📦 Input data
│   ├── membrane_protein_analysis_with_reduced_cc.csv
│   ├── SGTA and OP91 MS.csv
│   └── ...
│
├── filtered_datasets/                 # 📁 Output data
│   ├── proteins_cterm_30_in_massspec.csv
│   ├── proteins_cterm_distance_30.csv
│   └── proteins_high_membrane_domains.csv
│
├── start_app.sh                       # 🚀 Run web app (convenience script)
├── README.md                          # This file
└── .gitignore                         # Git ignore rules (optional)
```

## 🚀 Quick Start

### Local Development

```bash
# From root directory:
bash start_app.sh

# Or from app directory:
cd app
bash start_local.sh

# Visit: http://localhost:8000
```

### Deploy to Shinylive

```bash
# Install shinylive
pip install shinylive

# Export app
shinylive export app site

# Deploy to GitHub Pages, Netlify, Vercel, etc.
# See docs/SHINYLIVE_DEPLOYMENT.md for details
```

## 📚 Documentation

**Start here:**
- [docs/SHINYLIVE_DEPLOYMENT.md](docs/SHINYLIVE_DEPLOYMENT.md) - Deployment guide

**For sharing with colleagues:**
- [docs/USER_GUIDE.md](docs/USER_GUIDE.md) - Non-technical user guide

## ✨ Features

✅ Interactive filtering by C-terminal distance (operators: ≤, ≥, ==, <, >)
✅ Category comparison (membrane domains, predictions)
✅ Real-time bar charts with matplotlib
✅ Interactive data tables
✅ One-click CSV export
✅ Beautiful responsive design
✅ **Runs in browser** - no server needed!
✅ **Free deployment** - GitHub Pages, Netlify, etc.

## 🎯 What to Send to Colleagues

When sharing with colleagues:

1. **The URL** (e.g., `https://yourusername.github.io/TA-like-ID/`)
2. **[docs/USER_GUIDE.md](docs/USER_GUIDE.md)** - How to use the app

That's all! No installation, no Python, no setup required on their end.

## 📊 Analysis Scripts

The original analysis scripts are in the `analysis/` directory:
- Use these for batch processing or command-line analysis
- The web app uses these utilities internally

## 📓 Notebooks

Jupyter notebooks are preserved in `notebooks/` for reference.

## 🔧 Technology Stack

- **Frontend & Backend:** Shiny for Python
- **Deployment:** Shinylive (browser-based, no server!)
- **Data:** pandas, numpy
- **Visualization:** matplotlib
- **Hosting:** Static site (GitHub Pages, Netlify, etc.)

## 📦 Installation

```bash
# Install dependencies
cd app
pip install -r requirements.txt
```

## ⚠️ Notes

- Raw data files are required in `raw_data/` directory
- The app caches data on startup for performance
- All filtering is done in-memory (fast, even with large datasets)

## 📞 Support

See the documentation in `docs/` for:
- Architecture details
- Deployment instructions
- Troubleshooting
- User guides

---

**Ready to share your analysis with colleagues!** 🚀
