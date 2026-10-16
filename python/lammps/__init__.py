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

# convert installed module string version to numeric version
def get_version_number():
  """Extract LAMMPS version string and convert to number"""
  # pylint: disable=C0415
  from os import path
  from sys import version_info

  # must report 0 when inside LAMMPS source tree
  if __file__.find(path.join('python', 'lammps', '__init__.py')) > 0:
    return 0

  if version_info.major < 3 or (version_info.major == 3 and version_info.minor < 6):
    raise SystemError('LAMMPS only supports Python version 3.6 or later')

  vstring = None
  if version_info.major == 3 and version_info.minor >= 8:
    from importlib.metadata import version, PackageNotFoundError
    try:
      vstring = version('lammps')
    except PackageNotFoundError:
      # nothing to do, ignore
      pass

  else:
    from pkg_resources import get_distribution, DistributionNotFound
    try:
      vstring = get_distribution('lammps').version
    except DistributionNotFound:
      # nothing to do, ignore
      pass

  if not vstring:
    return 0

  # the day is incremented by one for development and maintenance versions (see setup.py),
  # so the version string is not always a valid date and must not be parsed as one
  year, month, day = (int(num) for num in vstring.split('.')[:3])
  return year*10000 + month*100 + day

__version__ = get_version_number()

# Local Variables:
# fill-column: 100
# python-indent-offset: 2
# End:
