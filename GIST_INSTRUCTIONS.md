# Share Your App via Gist 🚀

Your app is ready to share with colleagues! Here are two ways:

---

## **Option 1: Run Locally** (Easiest)

1. Copy the code from [app/app.py](app/app.py)
2. Create a new file `app.py` on their computer
3. Install dependencies:
   ```bash
   pip install shiny pandas matplotlib
   ```
4. Run locally:
   ```bash
   shiny run app.py
   ```
5. Visit `http://localhost:8000` in browser

---

## **Option 2: Create a GitHub Gist** (Best for Sharing)

1. Go to https://gist.github.com
2. Create a new gist with the code from [app/app.py](app/app.py)
3. Click "Create public gist"
4. Copy the gist URL: `https://gist.github.com/<username>/<gist-id>`

Share that link! Colleagues can either:
- **A) Run it locally** (copy paste the code, run `shiny run app.py`)
- **B) Run it in browser** (no installation needed):
  - Go to https://shinylive.io/py/editor/
  - Load gist using the gist URL
  - App runs entirely in their browser

---

## **Option 3: Shinylive + Static Hosting** (Zero-cost deployment)

Deploy to free hosting (GitHub Pages, Netlify, Vercel):

```bash
# Install shinylive
pip install shinylive

# Export the app
cd app
shinylive export . site

# Push 'site' folder to GitHub Pages or Netlify
```

Share the deployed URL with colleagues - app runs entirely in browser, no server needed!

---

## **Quick Sharing Checklist**

✅ Data loads from GitHub (no file paths needed)  
✅ All features work locally with `shiny run`  
✅ Can be deployed to Shinylive  
✅ CSV download works  
✅ Filtering & charting functional  

**Your app is ready to share!**
