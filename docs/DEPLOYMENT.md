# Deployment Guide

## Option 1: Shinylive (Recommended — no server needed)

The app can run entirely in the browser via [Shinylive](https://shiny.posit.co/py/docs/shinylive.html).

### Export and deploy

```bash
# From the project root
pip install shinylive
shinylive export app site
```

Then deploy the generated `site/` folder to any static host:

| Host | How |
|------|-----|
| **GitHub Pages** | Push `site/` contents to a `gh-pages` branch |
| **Netlify** | Drag-and-drop the `site/` folder |
| **Vercel** | Deploy the `site/` folder |
| **Any web server** | Upload `site/` folder contents |

After deployment your app is available at e.g. `https://yourusername.github.io/TA-like-ID/`.

### Data loading

The app loads its data from GitHub at startup using a public URL, so no additional file hosting is required.

---

## Option 2: PythonAnywhere (server-based)

### 1. Clone the repository

In the PythonAnywhere Bash console:

```bash
git clone https://github.com/j-a-hill/TA-like-ID.git
cd TA-like-ID
```

### 2. Create a virtual environment and install dependencies

```bash
mkvirtualenv --python=/usr/bin/python3.10 ta-like-id
workon ta-like-id
pip install -r app/requirements.txt
```

### 3. Configure the WSGI file

Go to **Web** tab → **Add a new web app** → **Manual configuration** → Python 3.10.

Find your WSGI file (usually `/var/www/username_pythonanywhere_com_wsgi.py`) and replace its contents with the contents of `app/wsgi_pythonanywhere.py`, then replace `username` with your PythonAnywhere username.

### 4. Set paths in the Web tab

| Setting | Value |
|---------|-------|
| Source code | `/home/username/TA-like-ID` |
| Working directory | `/home/username/TA-like-ID` |
| Virtualenv | `/home/username/.virtualenvs/ta-like-id` |

### 5. Reload and visit

Click the green **Reload** button. Your app will be at `https://username.pythonanywhere.com`.

---

## Local development

```bash
cd app
pip install -r requirements.txt
bash start_local.sh
# Visit http://localhost:8000
```

## Sharing with colleagues

Send them:
1. **The deployed URL**
2. **[docs/USER_GUIDE.md](USER_GUIDE.md)** — non-technical usage instructions

They need no Python installation — the app runs entirely in the browser (Shinylive) or on your server (PythonAnywhere).

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `ModuleNotFoundError` | Run `pip install -r app/requirements.txt` in the virtualenv |
| Data not loading | Check internet connection; the app loads from GitHub |
| App won't start locally | Ensure `shiny` is installed: `pip install shiny` |

