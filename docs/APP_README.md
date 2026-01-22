# TA-like-ID Web Application

## Overview

Your membrane protein analysis project is now wrapped into a beautiful, shareable web application! 🚀

**What's included:**
- ✅ Interactive web interface (no coding required)
- ✅ Dynamic filtering by C-terminal distance and other criteria
- ✅ Real-time category distribution charts
- ✅ CSV export for downstream analysis
- ✅ Mobile-responsive design
- ✅ Ready for deployment on PythonAnywhere

## Files Created

### Core Application
- **`app.py`** - Main Bottle web application with all routes and logic
- **`requirements.txt`** - Python dependencies for deployment
- **`wsgi_pythonanywhere.py`** - WSGI configuration for PythonAnywhere

### Documentation
- **`DEPLOYMENT.md`** - Step-by-step deployment guide for PythonAnywhere
- **`USER_GUIDE.md`** - Guide for non-technical users
- **`start_local.sh`** - Bash script for local testing

## Features

### Data Filtering
Users can filter proteins using:
- **C-terminal distance operators**: ≤, ≥, ==, <, >
- **Category comparison**: Membrane domain count, predictions, localization
- **Real-time results**: Instant updates as filters are applied

### Data Visualization
- **Category Distribution Bar Charts**: See the breakdown of your filtered data
- **Interactive charts**: Zoom, pan, and export using Plotly
- **Summary statistics**: Total proteins, filtered count, reduction percentage

### Data Export
- **CSV Download**: Export filtered results in standard CSV format
- **Excel compatible**: Opens directly in Excel, Google Sheets, etc.
- **Full data fidelity**: All columns preserved during export

## Quick Start (Local Testing)

```bash
cd /workspaces/TA-like-ID

# Option 1: Use the start script
bash start_local.sh

# Option 2: Manual setup
pip install -r requirements.txt
python app.py
```

Visit `http://localhost:8080` in your browser.

## Deployment to PythonAnywhere

1. **Sign up** at https://www.pythonanywhere.com (free account available)
2. Follow the detailed steps in `DEPLOYMENT.md`
3. Share the URL with colleagues

The entire process takes ~15 minutes. See `DEPLOYMENT.md` for detailed instructions.

## How to Share with Colleagues

### Step 1: Deploy
- Deploy to PythonAnywhere following `DEPLOYMENT.md`
- Get your unique URL (e.g., `https://username.pythonanywhere.com`)

### Step 2: Share
- Send colleagues the URL
- Send them `USER_GUIDE.md` for instructions
- No Python installation needed on their computers!

### Step 3: Support
- Users can filter, visualize, and export data independently
- They can send you filtered CSV files for collaboration
- No need to re-run analysis scripts

## Architecture

```
┌─────────────────────────────────────┐
│   User's Browser                    │
│   (Beautiful web interface)         │
└─────────────────────────────────────┘
              ↕ HTTP
┌─────────────────────────────────────┐
│   app.py (Bottle Web Server)        │
│   - /          (serve HTML/CSS/JS)  │
│   - /api/data  (get dataset)        │
│   - /api/filter (apply filters)     │
│   - /api/chart (generate charts)    │
└─────────────────────────────────────┘
              ↕
┌─────────────────────────────────────┐
│   Data Processing                   │
│   - pandas (data manipulation)      │
│   - protein_analysis_utils.py       │
│   - matplotlib (chart generation)   │
└─────────────────────────────────────┘
              ↕
┌─────────────────────────────────────┐
│   CSV Data Files                    │
│   raw_data/                         │
│   filtered_datasets/                │
└─────────────────────────────────────┘
```

## Technology Stack

- **Backend**: Python 3.10+, Bottle framework (micro-web framework)
- **Frontend**: HTML5, CSS3, JavaScript (vanilla, no frameworks needed)
- **Data Processing**: pandas, numpy
- **Visualization**: Plotly.js (interactive charts)
- **Server**: WSGI-compatible (works with Gunicorn on PythonAnywhere)

## Customization Options

### Add New Filter Columns
Edit `app.py` in the `buildFilterControls()` JavaScript function:
```javascript
<option value="your_column">Your Column Display Name</option>
```

### Change Colors
Modify the CSS gradients (current: purple-blue):
```css
background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
/* Change #667eea and #764ba2 to your brand colors */
```

### Add Statistics
Extend the `updateStats()` function to show more metrics

### Different Chart Types
Plotly supports: bar, line, scatter, pie, histogram, and more. Edit `updateChart()` to change chart type.

## Performance Notes

- **Data loading**: CSV is loaded once on startup and cached
- **Filtering speed**: Real-time, even with thousands of rows
- **Chart generation**: Instant with Plotly
- **CSV export**: Handles large datasets efficiently

## Troubleshooting

### "Data file not found"
Ensure your raw CSV is at: `raw_data/membrane_protein_analysis_with_reduced_cc.csv`

### Charts not showing
Make sure Plotly CDN is accessible (requires internet connection)

### Filters not working
Check browser console (F12 → Console tab) for JavaScript errors

### App won't start locally
Run: `pip install --upgrade -r requirements.txt`

## Next Steps

1. **Test locally** with `bash start_local.sh`
2. **Share test URL** with one colleague for feedback
3. **Deploy to PythonAnywhere** following `DEPLOYMENT.md`
4. **Share with team** via the public URL

## Support Resources

- **Bottle documentation**: https://bottlepy.org
- **PythonAnywhere help**: https://help.pythonanywhere.com
- **User guide for colleagues**: See `USER_GUIDE.md`
- **Deployment help**: See `DEPLOYMENT.md`

---

**Ready to share your research tool with the world!** 🌍
