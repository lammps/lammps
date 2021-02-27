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

# convert module string version to numeric version
def get_version_number():
    from datetime import datetime
    from sys import version_info
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

    d = datetime.strptime(vstring, "%d%b%Y")
    return d.year*10000 + d.month*100 + d.day

__version__ = get_version_number()
