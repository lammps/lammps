
import sys

if len(sys.argv) < 2:
    sys.exit("usage: %s <directory/with/lammps.py>" % sys.argv[0])

sys.path.insert(0,sys.argv[1])

from lammps import *
from ctypes import *

lmp=lammps()

