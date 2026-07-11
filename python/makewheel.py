#!/usr/bin/env python

import sys,os,site

base = os.path.abspath('buildwheel')
if sys.platform == 'win32':
  bin_dir=os.path.join(base,'Scripts')
else:
  bin_dir=os.path.join(base,'bin')

# prepend bin to PATH, set venv path
os.environ["PATH"] = os.pathsep.join([bin_dir] + os.environ.get("PATH", "").split(os.pathsep))
os.environ["VIRTUAL_ENV"] = base

# add the virtual environments libraries to the host python import mechanism
prev_length = len(sys.path)
for lib in "__LIB_FOLDERS__".split(os.pathsep):
    path = os.path.realpath(os.path.join(bin_dir, lib))
    site.addsitedir(path)
sys.path[:] = sys.path[prev_length:] + sys.path[0:prev_length]

sys.real_prefix = sys.prefix
sys.prefix = base

# update pip and install all requirements to build the wheel
os.system('python -m pip install --upgrade pip')
os.system('python -m pip install --upgrade -r wheel_requirements.txt')

print("Building new binary wheel")
os.system('python -m build -n --wheel -o .')
