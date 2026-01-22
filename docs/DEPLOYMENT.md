# Deployment Guide for PythonAnywhere

## Quick Start

1. **Create a PythonAnywhere account** at https://www.pythonanywhere.com

2. **Upload files via Git**:
   ```bash
   # In PythonAnywhere console:
   git clone https://github.com/j-a-hill/TA-like-ID.git
   cd TA-like-ID/app
   ```

3. **Create a virtual environment**:
   ```bash
   mkvirtualenv --python=/usr/bin/python3.10 ta-like-id
   workon ta-like-id
   ```

4. **Install dependencies** (from the app directory):
   ```bash
   pip install -r requirements.txt
   ```

5. **Test locally** (in PythonAnywhere console):
   ```bash
   cd /home/username/TA-like-ID/app
   python app.py
   ```

## Configuring PythonAnywhere Web App

1. Go to **Web** tab → **Add a new web app**
2. Choose **Manual configuration** → **Python 3.10**
3. Set the **WSGI configuration file** using the one in the project:
   - Find your WSGI file path (usually `/var/www/username_pythonanywhere_com_wsgi.py`)
   - Copy contents from: `app/wsgi_pythonanywhere.py`
   - Just remember to **replace 'username' with your actual PythonAnywhere username**

4. **Set working directory** in Web tab:
   - Source code: `/home/username/TA-like-ID`
   - WSGI file: `/var/www/username_pythonanywhere_com_wsgi.py`

5. **Reload web app** and visit your URL (e.g., `https://username.pythonanywhere.com`)

## Environment Configuration

Make sure your data files are in the correct locations:
- Raw data: `raw_data/membrane_protein_analysis_with_reduced_cc.csv`
- Output: `filtered_datasets/` (will be created automatically)

## Features

✅ **Filter by C-terminal distance** with operators (≤, ≥, ==, <, >)
✅ **Compare by categories** (Membrane Domain Count, Prediction, Localization)
✅ **Interactive bar charts** showing category distributions
✅ **CSV download** of filtered results
✅ **Real-time statistics** showing filtering results

## Sharing with Colleagues

Simply send them the PythonAnywhere URL. No Python installation required!

They can:
1. View the data summary
2. Adjust filter criteria with simple controls
3. See category distributions in real-time
4. Download filtered results as CSV for use in Excel or other tools

## Customization

### Add More Filter Options

Edit the `buildFilterControls()` function in the JavaScript section to add more columns:

```javascript
<option value="your_column">Your Column Name</option>
```

### Change Chart Type

Modify the `updateChart()` function to use different Plotly chart types (line, pie, scatter, etc.)

### Add More Statistics

Extend the stats display by adding more stat items in `updateStats()`

## Troubleshooting

**"Data file not found"**: Ensure CSV files are uploaded to the correct paths
**"Module not found"**: Run `pip install -r requirements.txt` in the virtualenv
**"Port already in use"**: Change the port number in `app.run()` or restart web app

## Local Development

To test locally before deploying:

```bash
pip install -r requirements.txt
python app.py
# Visit http://localhost:8080
```
