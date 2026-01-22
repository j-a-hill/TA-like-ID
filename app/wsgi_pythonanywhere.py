"""
WSGI wrapper for PythonAnywhere deployment.

On PythonAnywhere, point your WSGI config to this file.
Replace 'username' with your PythonAnywhere username.
"""

import sys
import os

# Add project directory to path
path = '/home/username/TA-like-ID'
if path not in sys.path:
    sys.path.append(path)

# Add app directory to path
app_path = os.path.join(path, 'app')
if app_path not in sys.path:
    sys.path.append(app_path)

# Change to app directory for relative file paths to work
os.chdir(path)

# Import the Bottle app from app/app.py
from app import app as application
