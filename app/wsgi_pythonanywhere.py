"""
WSGI wrapper for PythonAnywhere deployment.

On PythonAnywhere, point your WSGI config to this file.
Replace 'username' with your PythonAnywhere username.
"""

import sys
import os


def _configure_deployment_paths(project_root: str) -> None:
    """Add project directories to sys.path for PythonAnywhere deployment."""
    for path in (project_root, os.path.join(project_root, 'app')):
        if path not in sys.path:
            sys.path.insert(0, path)
    os.chdir(project_root)


_configure_deployment_paths('/home/username/TA-like-ID')

# Import the Shiny app from app/app.py
from app import app as application  # noqa: E402
