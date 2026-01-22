# PythonAnywhere Setup Configuration

## This file documents the exact configuration needed on PythonAnywhere

### Account Setup Checklist

- [ ] Create PythonAnywhere account at https://www.pythonanywhere.com
- [ ] Choose Python 3.10 or later
- [ ] Set username (you'll need this in steps below)
- [ ] Copy the repository or upload files

### Git Clone Method (Recommended)

In PythonAnywhere Bash console:
```bash
cd ~
git clone https://github.com/j-a-hill/TA-like-ID.git
cd TA-like-ID
```

### Virtual Environment Setup

In Bash console:
```bash
mkvirtualenv --python=/usr/bin/python3.10 ta-like-id
workon ta-like-id
cd ~/TA-like-ID/app
pip install -r requirements.txt
```

### Web App Configuration

Go to **Web** tab in PythonAnywhere dashboard:

#### Step 1: Add new web app
- Click "Add a new web app"
- Select "Manual configuration"
- Choose Python 3.10

#### Step 2: Configure WSGI

Find your WSGI configuration file path (shown in Web tab, typically):
```
/var/www/USERNAME_pythonanywhere_com_wsgi.py
```

Click on it and replace the entire content with:

```python
import sys
import os

# Add your project directory to the sys.path
path = '/home/USERNAME/TA-like-ID'
if path not in sys.path:
    sys.path.append(path)

# Add app directory to path
app_path = os.path.join(path, 'app')
if app_path not in sys.path:
    sys.path.append(app_path)

# Set working directory so relative paths work
os.chdir(path)

# Import and use the Bottle app from app/app.py
from app import app as application
```

**Replace USERNAME with your PythonAnywhere username**

#### Step 3: Configure source code

In Web tab, set:
- **Source code**: `/home/USERNAME/TA-like-ID`
- **Working directory**: `/home/USERNAME/TA-like-ID`
- **Virtualenv**: `/home/USERNAME/.virtualenvs/ta-like-id`

#### Step 4: Reload

Click the green "Reload" button at the top of the Web tab

#### Step 5: Visit your site

Your app will be at: `https://USERNAME.pythonanywhere.com`

### Data Files

Ensure these files exist in your project:
```
/home/USERNAME/TA-like-ID/
├── raw_data/
│   └── membrane_protein_analysis_with_reduced_cc.csv
├── filtered_datasets/  (created automatically on first use)
├── app.py
├── protein_analysis_utils.py
└── requirements.txt
```

### Testing the Setup

After reloading, visit your URL and test:
1. ✓ Page loads without errors
2. ✓ Data appears in the interface
3. ✓ Filter works (try ≤ 30)
4. ✓ Chart appears below
5. ✓ CSV download works

### Sharing with Colleagues

Once verified working, send them:
1. Your PythonAnywhere URL (e.g., https://myname.pythonanywhere.com)
2. Link to USER_GUIDE.md (in this repo)

### Common Issues and Fixes

#### "ModuleNotFoundError: No module named 'bottle'"
Solution:
```bash
workon ta-like-id
pip install bottle
```

#### "Data file not found"
Check that the CSV is in the right location:
```bash
ls -la raw_data/membrane_protein_analysis_with_reduced_cc.csv
```

#### Port/permission errors
Usually resolves after clicking "Reload" in Web tab

#### Static files not loading (styling broken)
Make sure static files are being served correctly. With Bottle, no special configuration needed.

### Upgrading Dependencies

If you need to update packages:
```bash
workon ta-like-id
pip install --upgrade -r requirements.txt
```

Then reload the web app.

### Logs

If something breaks, check logs in PythonAnywhere:
**Web** tab → scroll down to see error logs

Or use Bash console:
```bash
tail -50 /var/log/USERNAME.pythonanywhere.com.error.log
```

---

**Once deployed, your colleagues can use it immediately without any technical setup!**
