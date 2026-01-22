#!/bin/bash
# Quick start script for local Shiny development

echo "🧬 TA-like-ID Shiny App - Local Development"
echo "=========================================="
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed"
    exit 1
fi

echo "✓ Python 3 found"
python3 --version
echo ""

# Upgrade setuptools first
echo "📦 Preparing build tools..."
pip install -q --upgrade setuptools wheel 2>/dev/null

# Install dependencies
echo "📦 Installing dependencies..."
if pip install -q -r requirements.txt; then
    echo "✓ Dependencies installed"
else
    echo "❌ Failed to install dependencies"
    exit 1
fi

echo ""
echo "🚀 Starting Shiny app on http://localhost:8000"
echo "Press Ctrl+C to stop"
echo ""

shiny run --reload --port 8000 app.py
