"""
LAMMPS module global members:

.. data:: __version__

   Numerical representation of the LAMMPS version this
   module was taken from.  Has the same format as the
   result of :py:func:`lammps.version`.
"""

from .constants import *
from .core import *
from .data import *
from .pylammps import *

# convert installed module string version to numeric version
def get_version_number():
    import time
    from os.path import join
    from sys import version_info

    # must report 0 when inside LAMMPS source tree
    if __file__.find(join('python', 'lammps', '__init__.py')) > 0:
        return 0

    vstring = None
    if version_info.major == 3 and version_info.minor >= 8:
        from importlib.metadata import version
        try:
            vstring = version('lammps')
        except: pass
    else:
        from pkg_resources import get_distribution
        try:
            vstring = get_distribution('lammps').version
        except: pass

    if not vstring:
        return 0

    t = time.strptime(vstring, "%d%b%Y")
    return t.tm_year*10000 + t.tm_mon*100 + t.tm_mday

__version__ = get_version_number()
