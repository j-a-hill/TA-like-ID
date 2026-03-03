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

## 📝 Data Loading

Shinylive apps run in the browser, so data must be loaded from a URL rather than the local filesystem.

The app already handles this — it loads data via the `_DATA_URL` constant in `app/app.py`:

```python
_DATA_URL = "https://raw.githubusercontent.com/j-a-hill/TA-like-ID/main/raw_data/membrane_protein_analysis_with_reduced_cc.csv"
df = pd.read_csv(_DATA_URL)
```

No changes to `app.py` are needed before exporting with `shinylive export`.

## 🎯 Recommended Workflow

1. **Develop locally** with `shiny run`
2. **Export to Shinylive** with `shinylive export app site`
3. **Deploy `site/`** to GitHub Pages or another static host
4. **Share the URL** with colleagues

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

Perfect for sharing with colleagues who just need a simple interface!
