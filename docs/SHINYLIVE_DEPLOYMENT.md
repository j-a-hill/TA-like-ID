# Shinylive Deployment Guide

## Overview

This app is built with **Shiny for Python** and designed for deployment via **Shinylive**, which allows the app to run entirely in the browser using WebAssembly - no server needed!

## 🚀 Deployment Options

### Option 1: Shinylive.io (Quickest - Share via URL)

1. **Export your app to Shinylive format:**
   ```bash
   cd /path/to/TA-like-ID
   shinylive export app site
   ```

2. **Deploy to GitHub Pages:**
   ```bash
   # Create a gh-pages branch
   git checkout --orphan gh-pages
   git rm -rf .
   cp -r site/* .
   git add .
   git commit -m "Deploy Shinylive app"
   git push origin gh-pages
   ```

3. **Access your app:**
   `https://yourusername.github.io/TA-like-ID/`

### Option 2: Self-Hosted Static Site

1. **Export the app:**
   ```bash
   shinylive export app site
   ```

2. **Deploy to any static hosting:**
   - **Netlify**: Drag & drop the `site` folder
   - **Vercel**: Deploy the `site` folder
   - **GitHub Pages**: As shown above
   - **Any web server**: Upload `site` folder contents

### Option 3: Shinylive Editor (For Quick Testing/Sharing)

1. Visit https://shinylive.io/py/editor/

2. Copy your `app.py` content

3. Include the data:
   - Create a `data.csv` in the editor
   - Update the path in app.py to reference it

4. Share the URL that's generated

## 📦 Local Development

```bash
cd app
pip install -r requirements.txt
shiny run --reload --port 8000 app.py
```

Visit: http://localhost:8000

## 🌐 Features of Shinylive Deployment

✅ **No Server Required** - Runs entirely in the browser  
✅ **Free Hosting** - Use GitHub Pages, Netlify, etc.  
✅ **Fast** - No roundtrips to server  
✅ **Private** - Data never leaves user's browser  
✅ **Easy Sharing** - Just send a URL  

## 📝 Important Notes for Shinylive

### Data File Considerations

Shinylive apps run in the browser, so data files need to be handled specially:

**Option A: Embed data in the app**
```python
# Convert CSV to embedded Python dict or use pickle
import pickle
data = {...}  # Your data structure
```

**Option B: Load from URL**
```python
# Load data from a publicly accessible URL
df = pd.read_csv('https://raw.githubusercontent.com/user/repo/main/data.csv')
```

**Option C: User upload**
```python
# Add file upload widget
ui.input_file("upload", "Upload CSV", accept=".csv")
```

For this app, you'll want to either:
1. Host the CSV on GitHub and load it via URL
2. Embed a subset of the data in the app
3. Add a file upload feature

## 🔧 Modifying for Shinylive

Current app loads data from file system. For Shinylive deployment, update the data loading in `app.py`:

```python
# Instead of:
df = pd.read_csv(str(_data_file))

# Use one of:
# 1. Load from URL
df = pd.read_csv('https://raw.githubusercontent.com/j-a-hill/TA-like-ID/main/raw_data/membrane_protein_analysis_with_reduced_cc.csv')

# 2. Or add file upload
@reactive.Calc
def uploaded_data():
    file = input.upload()
    if file is None:
        return pd.DataFrame()  # Empty default
    return pd.read_csv(file[0]["datapath"])
```

## 🎯 Recommended Workflow

1. **Develop locally** with `shiny run`
2. **Test with data URL** or upload feature
3. **Export to Shinylive** with `shinylive export`
4. **Deploy to GitHub Pages** or other static host
5. **Share URL** with colleagues

No backend servers, no PythonAnywhere setup, no Docker - just pure browser-based Python!

## 📚 Resources

- Shinylive Documentation: https://shiny.posit.co/py/docs/shinylive.html
- Shiny for Python: https://shiny.posit.co/py/
- Example Apps: https://shinylive.io/py/examples/

## ✨ Benefits Over Server-Based Deployment

| Feature | Shinylive | Traditional Server |
|---------|-----------|-------------------|
| Cost | Free | Hosting fees |
| Setup | Minutes | Hours |
| Maintenance | None | Updates, security |
| Data Privacy | Browser only | Sent to server |
| Scalability | Infinite | Limited by server |
| Offline Use | Possible | No |

Perfect for sharing with colleagues who just need a simple interface!
