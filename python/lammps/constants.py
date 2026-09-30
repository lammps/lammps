# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

"""
various symbolic constants to be used
in certain calls to select data formats
"""

# all exported symbols
__all__ = ['LAMMPS_AUTODETECT', 'LAMMPS_NONE', 'LAMMPS_INT', 'LAMMPS_INT_2D', 'LAMMPS_DOUBLE',
           'LAMMPS_DOUBLE_2D', 'LAMMPS_INT64', 'LAMMPS_INT64_2D', 'LAMMPS_STRING',
           'LMP_STYLE_GLOBAL', 'LMP_STYLE_ATOM', 'LMP_STYLE_LOCAL', 'LMP_TYPE_SCALAR',
           'LMP_TYPE_VECTOR', 'LMP_TYPE_ARRAY', 'LMP_SIZE_VECTOR', 'LMP_SIZE_ROWS',
           'LMP_SIZE_COLS', 'LMP_ERROR_WARNING', 'LMP_ERROR_ONE', 'LMP_ERROR_ALL',
           'LMP_ERROR_WORLD', 'LMP_ERROR_UNIVERSE', 'LMP_VAR_EQUAL', 'LMP_VAR_ATOM',
           'LMP_VAR_VECTOR', 'LMP_VAR_STRING', 'LMP_BUFSIZE', 'get_ctypes_int', 'LMP_MAX_GROUP']

# these must be kept in sync with the enums in src/library.h, src/lmptype.h,
# tools/swig/lammps.i, examples/COUPLE/plugin/liblammpsplugin.h,
# and the constants in fortran/lammps.f90

LAMMPS_AUTODETECT  = None
LAMMPS_NONE        = -1
LAMMPS_INT         = 0
LAMMPS_INT_2D      = 1
LAMMPS_DOUBLE      = 2
LAMMPS_DOUBLE_2D   = 3
LAMMPS_INT64       = 4
LAMMPS_INT64_2D    = 5
LAMMPS_STRING      = 6

LMP_STYLE_GLOBAL   = 0
LMP_STYLE_ATOM     = 1
LMP_STYLE_LOCAL    = 2

LMP_TYPE_SCALAR    = 0
LMP_TYPE_VECTOR    = 1
LMP_TYPE_ARRAY     = 2
LMP_SIZE_VECTOR    = 3
LMP_SIZE_ROWS      = 4
LMP_SIZE_COLS      = 5

LMP_ERROR_WARNING  = 0
LMP_ERROR_ONE      = 1
LMP_ERROR_ALL      = 2
LMP_ERROR_WORLD    = 4
LMP_ERROR_UNIVERSE = 8

LMP_VAR_EQUAL      = 0
LMP_VAR_ATOM       = 1
LMP_VAR_VECTOR     = 2
LMP_VAR_STRING     = 3

# default buffer size for string buffers
LMP_BUFSIZE        = 1024
LMP_MAX_GROUP      = 32
# -------------------------------------------------------------------------

def get_ctypes_int(size):
  """return ctypes type matching the configured C/C++ integer size in LAMMPS"""
  # pylint: disable=C0415
  from ctypes import c_int, c_int32, c_int64
  if size == 4:
    return c_int32
  if size == 8:
    return c_int64
  return c_int

# Local Variables:
# fill-column: 100
# python-indent-offset: 2
# End:
