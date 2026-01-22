#!/usr/bin/env python3
"""
DEPLOYMENT CHECKLIST
====================

Use this checklist to ensure everything is ready for deployment.
"""

CHECKLIST = """
╔═══════════════════════════════════════════════════════════════════════════╗
║                    DEPLOYMENT CHECKLIST                                  ║
║                                                                           ║
║  Follow this checklist to ensure successful deployment to PythonAnywhere ║
╚═══════════════════════════════════════════════════════════════════════════╝

PRE-DEPLOYMENT (Do this first)
─────────────────────────────────────────────────────────────────────────────

□ 1. Verify Data Files
   [ ] raw_data/membrane_protein_analysis_with_reduced_cc.csv exists
   [ ] File is readable (not corrupted)
   [ ] File size is reasonable (> 1 KB)
   Command: ls -lh raw_data/*.csv

□ 2. Test Local App
   [ ] Install dependencies: pip install -r requirements.txt
   [ ] Run app: python app.py
   [ ] Visit: http://localhost:8080 (opens in browser)
   [ ] Page loads without errors
   [ ] Try filtering with distance ≤ 30
   [ ] Chart appears below stats
   [ ] CSV download works
   [ ] Reset button clears filters
   [ ] Stop with Ctrl+C

□ 3. Verify Python Version
   [ ] Python 3.10 or higher installed locally
   [ ] Command: python --version

□ 4. Check All Files Created
   [ ] app.py (main application)
   [ ] requirements.txt (dependencies)
   [ ] wsgi_pythonanywhere.py (PythonAnywhere config)
   [ ] protein_analysis_utils.py (analysis functions)
   [ ] All documentation files (*.md)

PYTHONANYWHERE SETUP
─────────────────────────────────────────────────────────────────────────────

□ 5. Create PythonAnywhere Account
   [ ] Go to https://www.pythonanywhere.com
   [ ] Sign up for free account
   [ ] Verify email
   [ ] Choose username (remember this!)

□ 6. Upload Repository via Git
   [ ] In PythonAnywhere Bash console, run:
       git clone https://github.com/j-a-hill/TA-like-ID.git
   [ ] Verify cloned: ls TA-like-ID/app.py

□ 7. Create Virtual Environment
   [ ] In Bash console, run:
       mkvirtualenv --python=/usr/bin/python3.10 ta-like-id
   [ ] Activate: workon ta-like-id
   [ ] Install: pip install -r requirements.txt
   [ ] No errors during installation

□ 8. Edit WSGI Configuration File
   [ ] In PythonAnywhere Web tab, find WSGI file path
   [ ] Opens /var/www/USERNAME_pythonanywhere_com_wsgi.py
   [ ] Replace entire content with code from wsgi_pythonanywhere.py
   [ ] ⚠️  IMPORTANT: Replace "USERNAME" with your PythonAnywhere username
   [ ] Save file

□ 9. Configure Web App Settings
   [ ] In Web tab, set:
       Source code: /home/USERNAME/TA-like-ID
       Working directory: /home/USERNAME/TA-like-ID
       Virtualenv: /home/USERNAME/.virtualenvs/ta-like-id
   [ ] All paths filled in
   [ ] All paths use YOUR username

□ 10. Reload Web App
   [ ] Click green "Reload" button at top of Web tab
   [ ] Wait for status to show "running"

TESTING DEPLOYMENT
─────────────────────────────────────────────────────────────────────────────

□ 11. Visit Your Site
   [ ] Go to: https://USERNAME.pythonanywhere.com
   [ ] Page loads (no 500 error)
   [ ] Layout appears correct
   [ ] Purple/blue colors visible

□ 12. Test Core Features
   [ ] Filter by distance ≤ 30
   [ ] Chart appears showing categories
   [ ] Stats show filtered count
   [ ] Download CSV button works
   [ ] Downloaded file opens in Excel/spreadsheet
   [ ] Reset button clears filters
   [ ] Try different operators (≥, ==, <, >)

□ 13. Check Browser Console
   [ ] Press F12 to open developer tools
   [ ] Click Console tab
   [ ] No red errors shown
   [ ] Close developer tools

SHARING WITH COLLEAGUES
─────────────────────────────────────────────────────────────────────────────

□ 14. Prepare Sharing
   [ ] Your URL: https://USERNAME.pythonanywhere.com (keep this!)
   [ ] Download USER_GUIDE.md from repository
   [ ] Email to colleagues with:
       - The URL
       - USER_GUIDE.md attachment
       - Short explanation of what it does

□ 15. Send Test Email
   [ ] Send to one colleague first
   [ ] Have them test it
   [ ] Collect feedback
   [ ] Make any needed adjustments

□ 16. Full Release
   [ ] Satisfied with feedback?
   [ ] Send to all colleagues
   [ ] Include instructions from USER_GUIDE.md
   [ ] Provide your contact for questions

MAINTENANCE (After Deployment)
─────────────────────────────────────────────────────────────────────────────

□ 17. Monitor Performance
   [ ] Check error logs weekly
   [ ] In Web tab → scroll to logs section
   [ ] Look for any error messages
   [ ] Address issues promptly

□ 18. Keep Dependencies Updated (Optional)
   [ ] Periodically check for updates:
       pip list --outdated
   [ ] Update if needed:
       pip install --upgrade -r requirements.txt
   [ ] Test thoroughly before deploying
   [ ] Reload web app in PythonAnywhere

□ 19. Backup
   [ ] Your GitHub repo is your backup
   [ ] Keep .csv files backed up locally
   [ ] Consider keeping deployment notes

TROUBLESHOOTING
─────────────────────────────────────────────────────────────────────────────

If things go wrong:

□ Check error logs
  In PythonAnywhere, see tail of error log:
  /var/log/USERNAME.pythonanywhere.com.error.log

□ Verify virtualenv activated
  In Bash: workon ta-like-id
  Should show (ta-like-id) $ prompt

□ Reinstall dependencies
  pip install --force-reinstall -r requirements.txt

□ Check WSGI file
  Make sure USERNAME is replaced with actual username

□ Reload web app
  Click "Reload" button in Web tab (green button)

□ Still broken?
  Check APP_README.md troubleshooting section

FINAL CHECKLIST
─────────────────────────────────────────────────────────────────────────────

All done? Run through this final checklist:

□ Local testing works perfectly
□ All files uploaded to PythonAnywhere
□ WSGI configured correctly (USERNAME filled in)
□ Web app reloaded
□ Site accessible at https://USERNAME.pythonanywhere.com
□ Filters work (tried multiple operators)
□ Charts display
□ CSV download works
□ Console has no errors (F12)
□ Colleagues have URL and USER_GUIDE.md
□ One colleague tested successfully

═══════════════════════════════════════════════════════════════════════════════

✅ COMPLETE! Your app is deployed and ready to use!

Share the URL with your colleagues and enjoy!
Questions? See the documentation files (*.md) in the repo.

═══════════════════════════════════════════════════════════════════════════════
"""

if __name__ == '__main__':
    print(CHECKLIST)
