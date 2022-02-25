#!/usr/bin/env python

import sys,os,shutil

if sys.platform == 'win32':
  virtenv=os.path.join('buildwheel','Scripts','activate_this.py')
else:
  virtenv=os.path.join('buildwheel','bin','activate_this.py')

exec(open(virtenv).read(), {'__file__': virtenv})

os.system('python -m pip install --upgrade pip')
os.system('python -m pip install --upgrade -r wheel_requirements.txt')
print("Building new binary wheel")
os.system('python -m build -n --wheel -o .')
