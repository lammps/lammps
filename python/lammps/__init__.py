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
    import re
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

    vregex = re.compile(r"([0-9]+)([A-Za-z]+)(2[0-9]+)")
    m = vregex.match(vstring)

    if (m):
        month2num = { 'Jan' : 1, 'Feb' : 2, 'Mar' : 3, 'Apr' : 4, 'May' : 5, 'Jun' : 6,
                'Jul' : 7, 'Aug' : 8, 'Sep' : 9, 'Oct' : 10, 'Nov' : 11, 'Dec' : 12 }
        try:
            vernum = int(m.group(3))*10000
            vernum += month2num[m.group(2)]*100
            vernum += int(m.group(1))
        except:
            exit('Failure to parse version string: %s' % verstr)
    return vernum

__version__ = get_version_number()
