#!/usr/bin/env python3
"""
QUICK REFERENCE - TA-like-ID Web App Deployment

Save this file for quick copy-paste commands
"""

# ============================================================================
# LOCAL TESTING (on your computer)
# ============================================================================

"""
STEP 1: Install dependencies
$ pip install -r requirements.txt

STEP 2: Run the app
$ python app.py

STEP 3: Open browser
Visit: http://localhost:8080

STEP 4: Test filtering
- Set distance to ≤ 30
- Select "Compare by: Membrane Domain Count"
- Click "Filter Data"
- See results appear
- Download CSV to test export
"""

# ============================================================================
# PYTHONANYWHERE DEPLOYMENT
# ============================================================================

"""
STEP 1: Create account
https://www.pythonanywhere.com (sign up)

STEP 2: Upload code via Git (in PythonAnywhere Bash console)
$ git clone https://github.com/j-a-hill/TA-like-ID.git
$ cd TA-like-ID

STEP 3: Create virtualenv (Bash console)
$ mkvirtualenv --python=/usr/bin/python3.10 ta-like-id
$ workon ta-like-id
$ pip install -r requirements.txt

STEP 4: Configure WSGI (in Web tab)
- Click "Add a new web app" → "Manual configuration" → Python 3.10
- Open the WSGI file shown
- Replace contents with code from wsgi_pythonanywhere.py (replace USERNAME)

STEP 5: Set paths in Web tab
- Source code: /home/USERNAME/TA-like-ID
- Working directory: /home/USERNAME/TA-like-ID  
- Virtualenv: /home/USERNAME/.virtualenvs/ta-like-id

STEP 6: Reload
- Click green "Reload" button

STEP 7: Share
- URL: https://USERNAME.pythonanywhere.com
- Send USER_GUIDE.md to colleagues
"""

# ============================================================================
# SHARING WITH NON-TECHNICAL COLLEAGUES
# ============================================================================

"""
1. Deploy to PythonAnywhere (steps above)
2. Get your URL: https://username.pythonanywhere.com
3. Send them:
   - The URL (paste in browser, no installation needed!)
   - USER_GUIDE.md for instructions
   
They can:
✓ Filter data with simple dropdowns
✓ See results in real-time charts
✓ Download filtered CSV for Excel
✓ NO PYTHON INSTALLATION REQUIRED
"""

# ============================================================================
# TROUBLESHOOTING QUICK FIXES
# ============================================================================

"""
App won't start:
$ pip install --upgrade -r requirements.txt
$ python app.py

CSV file not found:
$ ls raw_data/
(Check: membrane_protein_analysis_with_reduced_cc.csv exists)

Module errors on PythonAnywhere:
$ workon ta-like-id
$ pip install bottle pandas matplotlib

Charts not showing:
- Check browser console (F12)
- Try refreshing page
- Ensure internet connection (Plotly uses CDN)

Filter not working:
- Check browser console for errors
- Try resetting with Reset button
- Reload page
"""

# ============================================================================
# FILE LOCATIONS
# ============================================================================

"""
Local development:
/workspaces/TA-like-ID/app.py                    ← Main app
/workspaces/TA-like-ID/protein_analysis_utils.py ← Analysis logic
/workspaces/TA-like-ID/raw_data/*.csv            ← Data files

PythonAnywhere:
/home/USERNAME/TA-like-ID/app.py                 ← Main app
/home/USERNAME/TA-like-ID/raw_data/*.csv         ← Data files
/var/www/USERNAME_pythonanywhere_com_wsgi.py     ← WSGI config
"""

# ============================================================================
# WHAT USERS SEE
# ============================================================================

"""
Beautiful web interface with:

LEFT SIDE (Filters):
- Distance operator dropdown (≤, ≥, ==, <, >)
- Distance value input
- Category selector (what to compare)
- Filter button
- Reset button
- Download CSV button

RIGHT SIDE (Results):
- Summary statistics (total, filtered, reduction %)
- Bar chart (category distribution)
- Data table (first 100 rows)
"""

# ============================================================================
# FEATURES IMPLEMENTED
# ============================================================================

"""
✅ Web interface (HTML/CSS/JavaScript)
✅ Filter by C-terminal distance with operators
✅ Category comparison in filtered data
✅ Bar chart visualization (interactive)
✅ CSV export of filtered results
✅ Real-time statistics
✅ Mobile responsive design
✅ Ready for PythonAnywhere deployment
✅ User guide for non-technical users
"""

# ============================================================================
# TECH STACK
# ============================================================================

"""
Backend: Python 3.10+, Bottle framework
Frontend: HTML5, CSS3, Vanilla JavaScript
Data: pandas, numpy
Viz: Plotly.js, matplotlib
Server: WSGI-compatible (Gunicorn on PythonAnywhere)
"""

# ============================================================================
# DOCUMENTATION FILES
# ============================================================================

"""
APP_README.md              ← Overview and architecture
USER_GUIDE.md              ← For non-technical colleagues
DEPLOYMENT.md              ← Step-by-step deployment guide
PYTHONANYWHERE_CONFIG.md   ← Detailed config instructions
QUICK_START.py             ← This file
"""

# ============================================================================
# NEXT STEPS
# ============================================================================

"""
1. Test locally
   $ python app.py
   Visit http://localhost:8080

2. Deploy to PythonAnywhere
   Follow DEPLOYMENT.md or PYTHONANYWHERE_CONFIG.md

3. Share with colleagues
   - Send them your URL
   - Send them USER_GUIDE.md
   - They can start using it immediately!

4. Iterate
   - Get feedback
   - Add more columns/features if needed
   - Keep sharing link (data stays private on PythonAnywhere)
"""

print("✓ All files ready!")
print("✓ See DEPLOYMENT.md for PythonAnywhere setup")
print("✓ See USER_GUIDE.md for what to send colleagues")
