# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------
# Python wrappers for the LAMMPS library via ctypes

# for python2/3 compatibility

from __future__ import print_function

# imports for simple LAMMPS python wrapper module "lammps"

import sys,traceback,types
import warnings
from ctypes import *
from os.path import dirname,abspath,join
from inspect import getsourcefile

# imports for advanced LAMMPS python wrapper modules "PyLammps" and "IPyLammps"

from collections import namedtuple
import os
import select
import re
import sys

# various symbolic constants to be used
# in certain calls to select data formats
LAMMPS_AUTODETECT = None
LAMMPS_INT    = 0
LAMMPS_INT_2D  = 1
LAMMPS_DOUBLE = 2
LAMMPS_DOUBLE_2D = 3
LAMMPS_INT64 = 4
LAMMPS_INT64_2D = 5
LAMMPS_STRING = 6

# these must be kept in sync with the enums in library.h
LMP_STYLE_GLOBAL = 0
LMP_STYLE_ATOM   = 1
LMP_STYLE_LOCAL  = 2

LMP_TYPE_SCALAR  = 0
LMP_TYPE_VECTOR  = 1
LMP_TYPE_ARRAY   = 2
LMP_SIZE_VECTOR  = 3
LMP_SIZE_ROWS    = 4
LMP_SIZE_COLS    = 5

LMP_VAR_EQUAL = 0
LMP_VAR_ATOM  = 1

# -------------------------------------------------------------------------

def get_ctypes_int(size):
  if size == 4:
    return c_int32
  elif size == 8:
    return c_int64
  return c_int

# -------------------------------------------------------------------------

class MPIAbortException(Exception):
  def __init__(self, message):
    self.message = message

  def __str__(self):
    return repr(self.message)

# -------------------------------------------------------------------------

class NeighList:
    """This is a wrapper class that exposes the contents of a neighbor list.

    It can be used like a regular Python list. Each element is a tuple of:

    * the atom local index
    * its number of neighbors
    * and a pointer to an c_int array containing local atom indices of its
      neighbors

    Internally it uses the lower-level LAMMPS C-library interface.

    :param lmp: reference to instance of :py:class:`lammps`
    :type  lmp: lammps
    :param idx: neighbor list index
    :type  idx: int
    """
    def __init__(self, lmp, idx):
        self.lmp = lmp
        self.idx = idx

    def __str__(self):
        return "Neighbor List ({} atoms)".format(self.size)

    def __repr__(self):
        return self.__str__()

    @property
    def size(self):
        """
        :return: number of elements in neighbor list
        """
        return self.lmp.get_neighlist_size(self.idx)

    def get(self, element):
        """
        :return: tuple with atom local index, numpy array of neighbor local atom indices
        :rtype:  (int, int, ctypes.POINTER(c_int))
        """
        iatom, numneigh, neighbors = self.lmp.get_neighlist_element_neighbors(self.idx, element)
        return iatom, numneigh, neighbors

    # the methods below implement the iterator interface, so NeighList can be used like a regular Python list

    def __getitem__(self, element):
        return self.get(element)

    def __len__(self):
        return self.size

    def __iter__(self):
        inum = self.size

        for ii in range(inum):
            yield self.get(ii)

# -------------------------------------------------------------------------

class NumPyNeighList(NeighList):
    """This is a wrapper class that exposes the contents of a neighbor list.

    It can be used like a regular Python list. Each element is a tuple of:

    * the atom local index
    * a NumPy array containing the local atom indices of its neighbors

    Internally it uses the lower-level LAMMPS C-library interface.

    :param lmp: reference to instance of :py:class:`lammps`
    :type  lmp: lammps
    :param idx: neighbor list index
    :type  idx: int
    """
    def __init__(self, lmp, idx):
      super(NumPyNeighList, self).__init__(lmp, idx)

    def get(self, element):
        """
        :return: tuple with atom local index, numpy array of neighbor local atom indices
        :rtype:  (int, numpy.array)
        """
        iatom, neighbors = self.lmp.numpy.get_neighlist_element_neighbors(self.idx, element)
        return iatom, neighbors


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

class lammps(object):
  """Create an instance of the LAMMPS Python class.

  .. _mpi4py_docs: https://mpi4py.readthedocs.io/

  This is a Python wrapper class that exposes the LAMMPS C-library
  interface to Python.  It either requires that LAMMPS has been compiled
  as shared library which is then dynamically loaded via the ctypes
  Python module or that this module called from a Python function that
  is called from a Python interpreter embedded into a LAMMPS executable,
  for example through the :doc:`python invoke <python>` command.
  When the class is instantiated it calls the :cpp:func:`lammps_open`
  function of the LAMMPS C-library interface, which in
  turn will create an instance of the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
  C++ class.  The handle to this C++ class is stored internally
  and automatically passed to the calls to the C library interface.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm
  """

  # -------------------------------------------------------------------------
  # create an instance of LAMMPS

  def __init__(self,name='',cmdargs=None,ptr=None,comm=None):
    self.comm = comm
    self.opened = 0

    # determine module file location

    modpath = dirname(abspath(getsourcefile(lambda:0)))
    # for windows installers the shared library is in a different folder
    winpath = abspath(os.path.join(modpath,'..','bin'))
    self.lib = None
    self.lmp = None

    # if a pointer to a LAMMPS object is handed in
    # when being called from a Python interpreter
    # embedded into a LAMMPS executable, all library
    # symbols should already be available so we do not
    # load a shared object.

    try:
      if ptr: self.lib = CDLL("",RTLD_GLOBAL)
    except:
      self.lib = None

    # load liblammps.so unless name is given
    #   if name = "g++", load liblammps_g++.so
    # try loading the LAMMPS shared object from the location
    #   of lammps.py with an absolute path,
    #   so that LD_LIBRARY_PATH does not need to be set for regular install
    # fall back to loading with a relative path,
    #   typically requires LD_LIBRARY_PATH to be set appropriately

    if any([f.startswith('liblammps') and f.endswith('.dylib')
            for f in os.listdir(modpath)]):
      lib_ext = ".dylib"
    elif any([f.startswith('liblammps') and f.endswith('.dll')
              for f in os.listdir(modpath)]):
      lib_ext = ".dll"
    elif os.path.exists(winpath) and any([f.startswith('liblammps') and f.endswith('.dll')
                  for f in os.listdir(winpath)]):
      lib_ext = ".dll"
      modpath = winpath
    else:
      lib_ext = ".so"

    if not self.lib:
      try:
        if not name:
          self.lib = CDLL(join(modpath,"liblammps" + lib_ext),RTLD_GLOBAL)
        else:
          self.lib = CDLL(join(modpath,"liblammps_%s" % name + lib_ext),
                          RTLD_GLOBAL)
      except:
        if not name:
          self.lib = CDLL("liblammps" + lib_ext,RTLD_GLOBAL)
        else:
          self.lib = CDLL("liblammps_%s" % name + lib_ext,RTLD_GLOBAL)


    # declare all argument and return types for all library methods here.
    # exceptions are where the arguments depend on certain conditions and
    # then are defined where the functions are used.
    self.lib.lammps_extract_setting.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_setting.restype = c_int

    # set default types
    # needed in later declarations
    self.c_bigint = get_ctypes_int(self.extract_setting("bigint"))
    self.c_tagint = get_ctypes_int(self.extract_setting("tagint"))
    self.c_imageint = get_ctypes_int(self.extract_setting("imageint"))

    self.lib.lammps_open.restype = c_void_p
    self.lib.lammps_open_no_mpi.restype = c_void_p
    self.lib.lammps_close.argtypes = [c_void_p]
    self.lib.lammps_free.argtypes = [c_void_p]

    self.lib.lammps_file.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_file.restype = None

    self.lib.lammps_command.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_command.restype = c_char_p
    self.lib.lammps_commands_list.restype = None
    self.lib.lammps_commands_string.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_commands_string.restype = None

    self.lib.lammps_get_natoms.argtypes = [c_void_p]
    self.lib.lammps_get_natoms.restype = c_double
    self.lib.lammps_extract_box.argtypes = \
      [c_void_p,POINTER(c_double),POINTER(c_double),
       POINTER(c_double),POINTER(c_double),POINTER(c_double),
       POINTER(c_int),POINTER(c_int)]
    self.lib.lammps_extract_box.restype = None

    self.lib.lammps_reset_box.argtypes = \
      [c_void_p,POINTER(c_double),POINTER(c_double),c_double,c_double,c_double]
    self.lib.lammps_reset_box.restype = None

    self.lib.lammps_gather_atoms.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms.restype = None

    self.lib.lammps_gather_atoms_concat.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms_concat.restype = None

    self.lib.lammps_gather_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_atoms_subset.restype = None

    self.lib.lammps_scatter_atoms.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_scatter_atoms.restype = None

    self.lib.lammps_scatter_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_scatter_atoms_subset.restype = None

    self.lib.lammps_gather.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather.restype = None

    self.lib.lammps_gather_concat.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_concat.restype = None

    self.lib.lammps_gather_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_subset.restype = None

    self.lib.lammps_scatter.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_scatter.restype = None

    self.lib.lammps_scatter_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_scatter_subset.restype = None


    self.lib.lammps_find_pair_neighlist.argtypes = [c_void_p, c_char_p, c_int, c_int, c_int]
    self.lib.lammps_find_pair_neighlist.restype  = c_int

    self.lib.lammps_find_fix_neighlist.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_find_fix_neighlist.restype  = c_int

    self.lib.lammps_find_compute_neighlist.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_find_compute_neighlist.restype  = c_int

    self.lib.lammps_neighlist_num_elements.argtypes = [c_void_p, c_int]
    self.lib.lammps_neighlist_num_elements.restype  = c_int

    self.lib.lammps_neighlist_element_neighbors.argtypes = [c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int))]
    self.lib.lammps_neighlist_element_neighbors.restype  = None

    self.lib.lammps_is_running.argtypes = [c_void_p]
    self.lib.lammps_is_running.restype = c_int

    self.lib.lammps_force_timeout.argtypes = [c_void_p]

    self.lib.lammps_has_error.argtypes = [c_void_p]
    self.lib.lammps_has_error.restype = c_int

    self.lib.lammps_get_last_error_message.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_get_last_error_message.restype = c_int

    self.lib.lammps_extract_global.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_global_datatype.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_global_datatype.restype = c_int
    self.lib.lammps_extract_compute.argtypes = [c_void_p, c_char_p, c_int, c_int]

    self.lib.lammps_get_thermo.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_get_thermo.restype = c_double

    self.lib.lammps_encode_image_flags.restype = self.c_imageint

    self.lib.lammps_config_package_name.argtypes = [c_int, c_char_p, c_int]

    self.lib.lammps_set_variable.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_has_style.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_style_count.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_style_name.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int]

    self.lib.lammps_has_id.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_id_count.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_id_name.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int]

    self.lib.lammps_version.argtypes = [c_void_p]

    self.lib.lammps_get_os_info.argtypes = [c_char_p, c_int]

    self.lib.lammps_get_mpi_comm.argtypes = [c_void_p]

    self.lib.lammps_decode_image_flags.argtypes = [self.c_imageint, POINTER(c_int*3)]

    self.lib.lammps_extract_atom.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_atom_datatype.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_atom_datatype.restype = c_int

    self.lib.lammps_extract_fix.argtypes = [c_void_p, c_char_p, c_int, c_int, c_int, c_int]

    self.lib.lammps_extract_variable.argtypes = [c_void_p, c_char_p, c_char_p]

    # TODO: NOT IMPLEMENTED IN PYTHON WRAPPER
    self.lib.lammps_fix_external_set_energy_global = [c_void_p, c_char_p, c_double]
    self.lib.lammps_fix_external_set_virial_global = [c_void_p, c_char_p, POINTER(c_double)]

    # detect if Python is using version of mpi4py that can pass a communicator

    self.has_mpi4py = False
    try:
      from mpi4py import __version__ as mpi4py_version
      # tested to work with mpi4py versions 2 and 3
      self.has_mpi4py = mpi4py_version.split('.')[0] in ['2','3']
    except:
      pass

    # if no ptr provided, create an instance of LAMMPS
    #   don't know how to pass an MPI communicator from PyPar
    #   but we can pass an MPI communicator from mpi4py v2.0.0 and later
    #   no_mpi call lets LAMMPS use MPI_COMM_WORLD
    #   cargs = array of C strings from args
    # if ptr, then are embedding Python in LAMMPS input script
    #   ptr is the desired instance of LAMMPS
    #   just convert it to ctypes ptr and store in self.lmp

    if not ptr:

      # with mpi4py v2, can pass MPI communicator to LAMMPS
      # need to adjust for type of MPI communicator object
      # allow for int (like MPICH) or void* (like OpenMPI)
      if self.has_mpi4py and self.has_mpi_support:
        from mpi4py import MPI
        self.MPI = MPI

      if comm:
        if not self.has_mpi4py:
          raise Exception('Python mpi4py version is not 2 or 3')
        if not self.has_mpi_support:
          raise Exception('LAMMPS not compiled with real MPI library')
        if self.MPI._sizeof(self.MPI.Comm) == sizeof(c_int):
          MPI_Comm = c_int
        else:
          MPI_Comm = c_void_p

        narg = 0
        cargs = None
        if cmdargs:
          cmdargs.insert(0,"lammps.py")
          narg = len(cmdargs)
          for i in range(narg):
            if type(cmdargs[i]) is str:
              cmdargs[i] = cmdargs[i].encode()
          cargs = (c_char_p*narg)(*cmdargs)
          self.lib.lammps_open.argtypes = [c_int, c_char_p*narg, \
                                           MPI_Comm, c_void_p]
        else:
          self.lib.lammps_open.argtypes = [c_int, c_char_p, \
                                           MPI_Comm, c_void_p]

        self.opened = 1
        comm_ptr = self.MPI._addressof(comm)
        comm_val = MPI_Comm.from_address(comm_ptr)
        self.lmp = c_void_p(self.lib.lammps_open(narg,cargs,comm_val,None))

      else:
        if self.has_mpi4py and self.has_mpi_support:
          self.comm = self.MPI.COMM_WORLD
        self.opened = 1
        if cmdargs:
          cmdargs.insert(0,"lammps.py")
          narg = len(cmdargs)
          for i in range(narg):
            if type(cmdargs[i]) is str:
              cmdargs[i] = cmdargs[i].encode()
          cargs = (c_char_p*narg)(*cmdargs)
          self.lib.lammps_open_no_mpi.argtypes = [c_int, c_char_p*narg, \
                                                  c_void_p]
          self.lmp = c_void_p(self.lib.lammps_open_no_mpi(narg,cargs,None))
        else:
          self.lib.lammps_open_no_mpi.argtypes = [c_int, c_char_p, c_void_p]
          self.lmp = c_void_p(self.lib.lammps_open_no_mpi(0,None,None))

    else:
      # magic to convert ptr to ctypes ptr
      if sys.version_info >= (3, 0):
        # Python 3 (uses PyCapsule API)
        pythonapi.PyCapsule_GetPointer.restype = c_void_p
        pythonapi.PyCapsule_GetPointer.argtypes = [py_object, c_char_p]
        self.lmp = c_void_p(pythonapi.PyCapsule_GetPointer(ptr, None))
      else:
        # Python 2 (uses PyCObject API)
        pythonapi.PyCObject_AsVoidPtr.restype = c_void_p
        pythonapi.PyCObject_AsVoidPtr.argtypes = [py_object]
        self.lmp = c_void_p(pythonapi.PyCObject_AsVoidPtr(ptr))

    # optional numpy support (lazy loading)
    self._numpy = None

    self._installed_packages = None
    self._available_styles = None

    # add way to insert Python callback for fix external
    self.callback = {}
    self.FIX_EXTERNAL_CALLBACK_FUNC = CFUNCTYPE(None, py_object, self.c_bigint, c_int, POINTER(self.c_tagint), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)))
    self.lib.lammps_set_fix_external_callback.argtypes = [c_void_p, c_char_p, self.FIX_EXTERNAL_CALLBACK_FUNC, py_object]
    self.lib.lammps_set_fix_external_callback.restype = None

  # -------------------------------------------------------------------------
  # shut-down LAMMPS instance

  def __del__(self):
    if self.lmp and self.opened:
      self.lib.lammps_close(self.lmp)
      self.opened = 0

  # -------------------------------------------------------------------------

  @property
  def numpy(self):
    """ Return object to access numpy versions of API

    It provides alternative implementations of API functions that
    return numpy arrays instead of ctypes pointers. If numpy is not installed,
    accessing this property will lead to an ImportError.

    :return: instance of numpy wrapper object
    :rtype: numpy_wrapper
    """
    if not self._numpy:
      self._numpy = numpy_wrapper(self)
    return self._numpy

  # -------------------------------------------------------------------------

  def close(self):
    """Explicitly delete a LAMMPS instance through the C-library interface.

    This is a wrapper around the :cpp:func:`lammps_close` function of the C-library interface.
    """
    if self.opened: self.lib.lammps_close(self.lmp)
    self.lmp = None
    self.opened = 0

  # -------------------------------------------------------------------------

  def finalize(self):
    """Shut down the MPI communication through the library interface by calling :cpp:func:`lammps_finalize`.
    """
    if self.opened: self.lib.lammps_close(self.lmp)
    self.lmp = None
    self.opened = 0
    self.lib.lammps_finalize()

  # -------------------------------------------------------------------------

  def version(self):
    """Return a numerical representation of the LAMMPS version in use.

    This is a wrapper around the :cpp:func:`lammps_version` function of the C-library interface.

    :return: version number
    :rtype:  int
    """
    return self.lib.lammps_version(self.lmp)

  # -------------------------------------------------------------------------

  def get_os_info(self):
    """Return a string with information about the OS and compiler runtime

    This is a wrapper around the :cpp:func:`lammps_get_os_info` function of the C-library interface.

    :return: OS info string
    :rtype:  string
    """

    sb = create_string_buffer(512)
    self.lib.lammps_get_os_info(sb,512)
    return sb

  # -------------------------------------------------------------------------

  def get_mpi_comm(self):
    """Get the MPI communicator in use by the current LAMMPS instance

    This is a wrapper around the :cpp:func:`lammps_get_mpi_comm` function
    of the C-library interface.  It will return ``None`` if either the
    LAMMPS library was compiled without MPI support or the mpi4py
    Python module is not available.

    :return: MPI communicator
    :rtype:  MPI_Comm
    """

    if self.has_mpi4py and self.has_mpi_support:
        from mpi4py import MPI
        f_comm = self.lib.lammps_get_mpi_comm(self.lmp)
        c_comm = MPI.Comm.f2py(f_comm)
        return c_comm
    else:
        return None

  # -------------------------------------------------------------------------

  @property
  def _lammps_exception(self):
    sb = create_string_buffer(100)
    error_type = self.lib.lammps_get_last_error_message(self.lmp, sb, 100)
    error_msg = sb.value.decode().strip()

    if error_type == 2:
      return MPIAbortException(error_msg)
    return Exception(error_msg)

  # -------------------------------------------------------------------------

  def file(self, path):
    """Read LAMMPS commands from a file.

    This is a wrapper around the :cpp:func:`lammps_file` function of the C-library interface.
    It will open the file with the name/path `file` and process the LAMMPS commands line by line until
    the end. The function will return when the end of the file is reached.

    :param path: Name of the file/path with LAMMPS commands
    :type path:  string
    """
    if path: path = path.encode()
    else: return
    self.lib.lammps_file(self.lmp, path)

    if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
      raise self._lammps_exception

  # -------------------------------------------------------------------------

  def command(self,cmd):
    """Process a single LAMMPS input command from a string.

    This is a wrapper around the :cpp:func:`lammps_command`
    function of the C-library interface.

    :param cmd: a single lammps command
    :type cmd:  string
    """
    if cmd: cmd = cmd.encode()
    else: return
    self.lib.lammps_command(self.lmp,cmd)

    if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
      raise self._lammps_exception

  # -------------------------------------------------------------------------

  def commands_list(self,cmdlist):
    """Process multiple LAMMPS input commands from a list of strings.

    This is a wrapper around the
    :cpp:func:`lammps_commands_list` function of
    the C-library interface.

    :param cmdlist: a single lammps command
    :type cmdlist:  list of strings
    """
    cmds = [x.encode() for x in cmdlist if type(x) is str]
    narg = len(cmdlist)
    args = (c_char_p * narg)(*cmds)
    self.lib.lammps_commands_list.argtypes = [c_void_p, c_int, c_char_p * narg]
    self.lib.lammps_commands_list(self.lmp,narg,args)

    if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
      raise self._lammps_exception

  # -------------------------------------------------------------------------

  def commands_string(self,multicmd):
    """Process a block of LAMMPS input commands from a string.

    This is a wrapper around the
    :cpp:func:`lammps_commands_string`
    function of the C-library interface.

    :param multicmd: text block of lammps commands
    :type multicmd:  string
    """
    if type(multicmd) is str: multicmd = multicmd.encode()
    self.lib.lammps_commands_string(self.lmp,c_char_p(multicmd))

    if self.has_exceptions and self.lib.lammps_has_error(self.lmp):
      raise self._lammps_exception

  # -------------------------------------------------------------------------

  def get_natoms(self):
    """Get the total number of atoms in the LAMMPS instance.

    Will be precise up to 53-bit signed integer due to the
    underlying :cpp:func:`lammps_get_natoms` function returning a double.

    :return: number of atoms
    :rtype: int
    """
    return int(self.lib.lammps_get_natoms(self.lmp))

  # -------------------------------------------------------------------------

  def extract_box(self):
    """Extract simulation box parameters

    This is a wrapper around the :cpp:func:`lammps_extract_box` function
    of the C-library interface.  Unlike in the C function, the result is
    returned as a list.

    :return: list of the extracted data: boxlo, boxhi, xy, yz, xz, periodicity, box_change
    :rtype: [ 3*double, 3*double, double, double, 3*int, int]
    """
    boxlo = (3*c_double)()
    boxhi = (3*c_double)()
    xy = c_double()
    yz = c_double()
    xz = c_double()
    periodicity = (3*c_int)()
    box_change = c_int()

    self.lib.lammps_extract_box(self.lmp,boxlo,boxhi,
                                byref(xy),byref(yz),byref(xz),
                                periodicity,byref(box_change))

    boxlo = boxlo[:3]
    boxhi = boxhi[:3]
    xy = xy.value
    yz = yz.value
    xz = xz.value
    periodicity = periodicity[:3]
    box_change = box_change.value

    return boxlo,boxhi,xy,yz,xz,periodicity,box_change

  # -------------------------------------------------------------------------

  def reset_box(self,boxlo,boxhi,xy,yz,xz):
    """Reset simulation box parameters

    This is a wrapper around the :cpp:func:`lammps_reset_box` function
    of the C-library interface.

    :param boxlo: new lower box boundaries
    :type boxlo: list of 3 floating point numbers
    :param boxhi: new upper box boundaries
    :type boxhi: list of 3 floating point numbers
    :param xy: xy tilt factor
    :type xy: float
    :param yz: yz tilt factor
    :type yz: float
    :param xz: xz tilt factor
    :type xz: float
    """
    cboxlo = (3*c_double)(*boxlo)
    cboxhi = (3*c_double)(*boxhi)
    self.lib.lammps_reset_box(self.lmp,cboxlo,cboxhi,xy,yz,xz)

  # -------------------------------------------------------------------------

  def get_thermo(self,name):
    """Get current value of a thermo keyword

    This is a wrapper around the :cpp:func:`lammps_get_thermo`
    function of the C-library interface.

    :param name: name of thermo keyword
    :type name: string
    :return: value of thermo keyword
    :rtype: double or None
    """
    if name: name = name.encode()
    else: return None
    return self.lib.lammps_get_thermo(self.lmp,name)

  # -------------------------------------------------------------------------

  def extract_setting(self, name):
    """Query LAMMPS about global settings that can be expressed as an integer.

    This is a wrapper around the :cpp:func:`lammps_extract_setting`
    function of the C-library interface.  Its documentation includes
    a list of the supported keywords.

    :param name: name of the setting
    :type name:  string
    :return: value of the setting
    :rtype: int
    """
    if name: name = name.encode()
    else: return None
    return int(self.lib.lammps_extract_setting(self.lmp,name))

  # -------------------------------------------------------------------------
  # extract global info datatype

  def extract_global_datatype(self, name):
    """Retrieve global property datatype from LAMMPS

    This is a wrapper around the :cpp:func:`lammps_extract_global_datatype`
    function of the C-library interface. Its documentation includes a
    list of the supported keywords.
    This function returns ``None`` if the keyword is not
    recognized. Otherwise it will return a positive integer value that
    corresponds to one of the :ref:`data type <py_datatype_constants>`
    constants define in the :py:mod:`lammps` module.

    :param name: name of the property
    :type name:  string
    :return: data type of global property, see :ref:`py_datatype_constants`
    :rtype: int
    """
    if name: name = name.encode()
    else: return None
    return self.lib.lammps_extract_global_datatype(self.lmp, name)

  # -------------------------------------------------------------------------
  # extract global info

  def extract_global(self, name, dtype=LAMMPS_AUTODETECT):
    """Query LAMMPS about global settings of different types.

    This is a wrapper around the :cpp:func:`lammps_extract_global`
    function of the C-library interface.  Unlike the C function
    this method returns the value and not a pointer and thus can
    only return the first value for keywords representing a list
    of values.  The :cpp:func:`lammps_extract_global` documentation
    includes a list of the supported keywords and their data types.
    Since Python needs to know the data type to be able to interpret
    the result, by default, this function will try to auto-detect the data type
    by asking the library. You can also force a specific data type.  For that
    purpose the :py:mod:`lammps` module contains :ref:`data type <py_datatype_constants>`
    constants. This function returns ``None`` if either the keyword is not recognized,
    or an invalid data type constant is used.

    :param name: name of the property
    :type name:  string
    :param dtype: data type of the returned data (see :ref:`py_datatype_constants`)
    :type dtype:  int, optional
    :return: value of the property or None
    :rtype: int, float, or NoneType
    """
    if dtype == LAMMPS_AUTODETECT:
      dtype = self.extract_global_datatype(name)

    if name: name = name.encode()
    else: return None

    if dtype == LAMMPS_INT:
      self.lib.lammps_extract_global.restype = POINTER(c_int32)
      target_type = int
    elif dtype == LAMMPS_INT64:
      self.lib.lammps_extract_global.restype = POINTER(c_int64)
      target_type = int
    elif dtype == LAMMPS_DOUBLE:
      self.lib.lammps_extract_global.restype = POINTER(c_double)
      target_type = float
    elif dtype == LAMMPS_STRING:
      self.lib.lammps_extract_global.restype = c_char_p
      target_type = lambda x: str(x, 'ascii')

    ptr = self.lib.lammps_extract_global(self.lmp, name)
    if ptr:
      return target_type(ptr[0])
    return None


  # -------------------------------------------------------------------------
  # extract per-atom info datatype

  def extract_atom_datatype(self, name):
    """Retrieve per-atom property datatype from LAMMPS

    This is a wrapper around the :cpp:func:`lammps_extract_atom_datatype`
    function of the C-library interface. Its documentation includes a
    list of the supported keywords.
    This function returns ``None`` if the keyword is not
    recognized. Otherwise it will return an integer value that
    corresponds to one of the :ref:`data type <py_datatype_constants>` constants
    defined in the :py:mod:`lammps` module.

    :param name: name of the property
    :type name:  string
    :return: data type of per-atom property (see :ref:`py_datatype_constants`)
    :rtype: int
    """
    if name: name = name.encode()
    else: return None
    return self.lib.lammps_extract_atom_datatype(self.lmp, name)

  # -------------------------------------------------------------------------
  # extract per-atom info

  def extract_atom(self, name, dtype=LAMMPS_AUTODETECT):
    """Retrieve per-atom properties from LAMMPS

    This is a wrapper around the :cpp:func:`lammps_extract_atom`
    function of the C-library interface. Its documentation includes a
    list of the supported keywords and their data types.
    Since Python needs to know the data type to be able to interpret
    the result, by default, this function will try to auto-detect the data type
    by asking the library. You can also force a specific data type by setting ``dtype``
    to one of the :ref:`data type <py_datatype_constants>` constants defined in the
    :py:mod:`lammps` module.
    This function returns ``None`` if either the keyword is not
    recognized, or an invalid data type constant is used.

    .. note::

       While the returned arrays of per-atom data are dimensioned
       for the range [0:nmax] - as is the underlying storage -
       the data is usually only valid for the range of [0:nlocal],
       unless the property of interest is also updated for ghost
       atoms.  In some cases, this depends on a LAMMPS setting, see
       for example :doc:`comm_modify vel yes <comm_modify>`.

    :param name: name of the property
    :type name:  string
    :param dtype: data type of the returned data (see :ref:`py_datatype_constants`)
    :type dtype:  int, optional
    :return: requested data or ``None``
    :rtype: ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.POINTER(ctypes.c_int64)),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
            or NoneType
    """
    if dtype == LAMMPS_AUTODETECT:
      dtype = self.extract_atom_datatype(name)

    if name: name = name.encode()
    else: return None

    if dtype == LAMMPS_INT:
      self.lib.lammps_extract_atom.restype = POINTER(c_int32)
    elif dtype == LAMMPS_INT_2D:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int32))
    elif dtype == LAMMPS_DOUBLE:
      self.lib.lammps_extract_atom.restype = POINTER(c_double)
    elif dtype == LAMMPS_DOUBLE_2D:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
    elif dtype == LAMMPS_INT64:
      self.lib.lammps_extract_atom.restype = POINTER(c_int64)
    elif dtype == LAMMPS_INT64_2D:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_int64))
    else: return None
    ptr = self.lib.lammps_extract_atom(self.lmp, name)
    if ptr: return ptr
    else:   return None


  # -------------------------------------------------------------------------

  def extract_compute(self,id,style,type):
    """Retrieve data from a LAMMPS compute

    This is a wrapper around the :cpp:func:`lammps_extract_compute`
    function of the C-library interface.
    This function returns ``None`` if either the compute id is not
    recognized, or an invalid combination of :ref:`style <py_style_constants>`
    and :ref:`type <py_type_constants>` constants is used. The
    names and functionality of the constants are the same as for
    the corresponding C-library function.  For requests to return
    a scalar or a size, the value is returned, otherwise a pointer.

    :param id: compute ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type style:  int
    :param type: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type type:  int
    :return: requested data as scalar, pointer to 1d or 2d double array, or None
    :rtype: c_double, ctypes.POINTER(c_double), ctypes.POINTER(ctypes.POINTER(c_double)), or NoneType
    """
    if id: id = id.encode()
    else: return None

    if type == LMP_TYPE_SCALAR:
      if style == LMP_STYLE_GLOBAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_double)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]
      elif style == LMP_STYLE_ATOM:
        return None
      elif style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    if type == LMP_TYPE_VECTOR:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr

    if type == LMP_TYPE_ARRAY:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr

    if type == LMP_SIZE_COLS:
      if style == LMP_STYLE_GLOBAL  \
         or style == LMP_STYLE_ATOM \
         or style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    if type == LMP_SIZE_VECTOR  \
       or type == LMP_SIZE_ROWS:
      if style == LMP_STYLE_GLOBAL  \
         or style == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
        return ptr[0]

    return None

  # -------------------------------------------------------------------------
  # extract fix info
  # in case of global data, free memory for 1 double via lammps_free()
  # double was allocated by library interface function

  def extract_fix(self,id,style,type,nrow=0,ncol=0):
    """Retrieve data from a LAMMPS fix

    This is a wrapper around the :cpp:func:`lammps_extract_fix`
    function of the C-library interface.
    This function returns ``None`` if either the fix id is not
    recognized, or an invalid combination of :ref:`style <py_style_constants>`
    and :ref:`type <py_type_constants>` constants is used. The
    names and functionality of the constants are the same as for
    the corresponding C-library function.  For requests to return
    a scalar or a size, the value is returned, also when accessing
    global vectors or arrays, otherwise a pointer.

    :param id: fix ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type style:  int
    :param type: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type type:  int
    :param nrow: index of global vector element or row index of global array element
    :type nrow:  int
    :param ncol: column index of global array element
    :type ncol:  int
    :return: requested data or None
    :rtype: c_double, ctypes.POINTER(c_double), ctypes.POINTER(ctypes.POINTER(c_double)), or NoneType

    """
    if id: id = id.encode()
    else: return None

    if style == LMP_STYLE_GLOBAL:
      if type in (LMP_TYPE_SCALAR, LMP_TYPE_VECTOR, LMP_TYPE_ARRAY):
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
        ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
        result = ptr[0]
        self.lib.lammps_free(ptr)
        return result
      elif type in (LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS):
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
        ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
        return ptr[0]
      else:
        return None

    elif style == LMP_STYLE_ATOM:
      if type == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif type == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif type == LMP_SIZE_COLS:
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
      if type == LMP_SIZE_COLS:
        return ptr[0]
      else:
        return ptr

    elif style == LMP_STYLE_LOCAL:
      if type == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif type == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif type in (LMP_TYPE_SCALAR, LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS):
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,nrow,ncol)
      if type in (LMP_TYPE_VECTOR, LMP_TYPE_ARRAY):
        return ptr
      else:
        return ptr[0]
    else:
      return None

  # -------------------------------------------------------------------------
  # extract variable info
  # free memory for 1 double or 1 vector of doubles via lammps_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function

  def extract_variable(self, name, group=None, vartype=LMP_VAR_EQUAL):
    """ Evaluate a LAMMPS variable and return its data

    This function is a wrapper around the function
    :cpp:func:`lammps_extract_variable` of the C-library interface,
    evaluates variable name and returns a copy of the computed data.
    The memory temporarily allocated by the C-interface is deleted
    after the data is copied to a Python variable or list.
    The variable must be either an equal-style (or equivalent)
    variable or an atom-style variable. The variable type has to
    provided as ``vartype`` parameter which may be one of two constants:
    ``LMP_VAR_EQUAL`` or ``LMP_VAR_STRING``; it defaults to
    equal-style variables.
    The group parameter is only used for atom-style variables and
    defaults to the group "all" if set to ``None``, which is the default.

    :param name: name of the variable to execute
    :type name: string
    :param group: name of group for atom-style variable
    :type group: string, only for atom-style variables
    :param vartype: type of variable, see :ref:`py_vartype_constants`
    :type vartype: int
    :return: the requested data
    :rtype: c_double, (c_double), or NoneType
    """
    if name: name = name.encode()
    else: return None
    if group: group = group.encode()
    if vartype == LMP_VAR_EQUAL:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      if ptr: result = ptr[0]
      else: return None
      self.lib.lammps_free(ptr)
      return result
    elif vartype == LMP_VAR_ATOM:
      nlocal = self.extract_global("nlocal")
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      if ptr:
        for i in range(nlocal): result[i] = ptr[i]
        self.lib.lammps_free(ptr)
      else: return None
      return result
    return None

  # -------------------------------------------------------------------------

  def set_variable(self,name,value):
    """Set a new value for a LAMMPS string style variable

    This is a wrapper around the :cpp:func:`lammps_set_variable`
    function of the C-library interface.

    :param name: name of the variable
    :type name: string
    :param value: new variable value
    :type value: any. will be converted to a string
    :return: either 0 on success or -1 on failure
    :rtype: int
    """
    if name: name = name.encode()
    else: return -1
    if value: value = str(value).encode()
    else: return -1
    return self.lib.lammps_set_variable(self.lmp,name,value)

  # -------------------------------------------------------------------------

  # return vector of atom properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def gather_atoms(self,name,type,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms(self.lmp,name,type,count,data)
    else: return None
    return data

  # -------------------------------------------------------------------------

  def gather_atoms_concat(self,name,type,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_atoms_concat(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_atoms_concat(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_atoms_subset(self,name,type,count,ndata,ids):
    if name: name = name.encode()
    if type == 0:
      data = ((count*ndata)*c_int)()
      self.lib.lammps_gather_atoms_subset(self.lmp,name,type,count,ndata,ids,data)
    elif type == 1:
      data = ((count*ndata)*c_double)()
      self.lib.lammps_gather_atoms_subset(self.lmp,name,type,count,ndata,ids,data)
    else: return None
    return data

  # -------------------------------------------------------------------------

  # scatter vector of atom properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter_atoms(self,name,type,count,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_atoms(self.lmp,name,type,count,data)

  # -------------------------------------------------------------------------

  def scatter_atoms_subset(self,name,type,count,ndata,ids,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_atoms_subset(self.lmp,name,type,count,ndata,ids,data)

  # return vector of atom/compute/fix properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes
  def gather(self,name,type,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_concat(self,name,type,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    if type == 0:
      data = ((count*natoms)*c_int)()
      self.lib.lammps_gather_concat(self.lmp,name,type,count,data)
    elif type == 1:
      data = ((count*natoms)*c_double)()
      self.lib.lammps_gather_concat(self.lmp,name,type,count,data)
    else: return None
    return data

  def gather_subset(self,name,type,count,ndata,ids):
    if name: name = name.encode()
    if type == 0:
      data = ((count*ndata)*c_int)()
      self.lib.lammps_gather_subset(self.lmp,name,type,count,ndata,ids,data)
    elif type == 1:
      data = ((count*ndata)*c_double)()
      self.lib.lammps_gather_subset(self.lmp,name,type,count,ndata,ids,data)
    else: return None
    return data

  # scatter vector of atom/compute/fix properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to insure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter(self,name,type,count,data):
    if name: name = name.encode()
    self.lib.lammps_scatter(self.lmp,name,type,count,data)

  def scatter_subset(self,name,type,count,ndata,ids,data):
    if name: name = name.encode()
    self.lib.lammps_scatter_subset(self.lmp,name,type,count,ndata,ids,data)

   # -------------------------------------------------------------------------

  def encode_image_flags(self,ix,iy,iz):
    """ convert 3 integers with image flags for x-, y-, and z-direction
    into a single integer like it is used internally in LAMMPS

    This method is a wrapper around the :cpp:func:`lammps_encode_image_flags`
    function of library interface.

    :param ix: x-direction image flag
    :type  ix: int
    :param iy: y-direction image flag
    :type  iy: int
    :param iz: z-direction image flag
    :type  iz: int
    :return: encoded image flags
    :rtype: lammps.c_imageint
    """
    return self.lib.lammps_encode_image_flags(ix,iy,iz)

  # -------------------------------------------------------------------------

  def decode_image_flags(self,image):
    """ Convert encoded image flag integer into list of three regular integers.

    This method is a wrapper around the :cpp:func:`lammps_decode_image_flags`
    function of library interface.

    :param image: encoded image flags
    :type image:  lammps.c_imageint
    :return: list of three image flags in x-, y-, and z- direction
    :rtype: list of 3 int
    """

    flags = (c_int*3)()
    self.lib.lammps_decode_image_flags(image,byref(flags))

    return [int(i) for i in flags]

  # -------------------------------------------------------------------------

  # create N atoms on all procs
  # N = global number of atoms
  # id = ID of each atom (optional, can be None)
  # type = type of each atom (1 to Ntypes) (required)
  # x = coords of each atom as (N,3) array (required)
  # v = velocity of each atom as (N,3) array (optional, can be None)
  # NOTE: how could we insure are passing correct type to LAMMPS
  #   e.g. for Python list or NumPy, etc
  #   ditto for gather_atoms() above

  def create_atoms(self,n,id,type,x,v=None,image=None,shrinkexceed=False):
    """
    Create N atoms from list of coordinates and properties

    This function is a wrapper around the :cpp:func:`lammps_create_atoms`
    function of the C-library interface, and the behavior is similar except
    that the *v*, *image*, and *shrinkexceed* arguments are optional and
    default to *None*, *None*, and *False*, respectively. With none being
    equivalent to a ``NULL`` pointer in C.

    The lists of coordinates, types, atom IDs, velocities, image flags can
    be provided in any format that may be converted into the required
    internal data types.  Also the list may contain more than *N* entries,
    but not fewer.  In the latter case, the function will return without
    attempting to create atoms.  You may use the :py:func:`encode_image_flags
    <lammps.encode_image_flags>` method to properly combine three integers
    with image flags into a single integer.

    :param n: number of atoms for which data is provided
    :type n: int
    :param id: list of atom IDs with at least n elements or None
    :type id: list of lammps.tagint
    :param type: list of atom types
    :type type: list of int
    :param x: list of coordinates for x-, y-, and z (flat list of 3n entries)
    :type x: list of float
    :param v: list of velocities for x-, y-, and z (flat list of 3n entries) or None (optional)
    :type v: list of float
    :param image: list of encoded image flags (optional)
    :type image: list of lammps.imageint
    :param shrinkexceed: whether to expand shrink-wrap boundaries if atoms are outside the box (optional)
    :type shrinkexceed: bool
    :return: number of atoms created. 0 if insufficient or invalid data
    :rtype: int
    """
    if id:
      id_lmp = (self.c_tagint*n)()
      try:
        id_lmp[:] = id[0:n]
      except:
        return 0
    else:
      id_lmp = None

    type_lmp = (c_int*n)()
    try:
      type_lmp[:] = type[0:n]
    except:
      return 0

    three_n = 3*n
    x_lmp = (c_double*three_n)()
    try:
      x_lmp[:] = x[0:three_n]
    except:
      return 0

    if v:
      v_lmp = (c_double*(three_n))()
      try:
        v_lmp[:] = v[0:three_n]
      except:
        return 0
    else:
      v_lmp = None

    if image:
      img_lmp = (self.c_imageint*n)()
      try:
        img_lmp[:] = image[0:n]
      except:
        return 0
    else:
      img_lmp = None

    if shrinkexceed:
      se_lmp = 1
    else:
      se_lmp = 0

    self.lib.lammps_create_atoms.argtypes = [c_void_p, c_int, POINTER(self.c_tagint*n),
                                     POINTER(c_int*n), POINTER(c_double*three_n),
                                     POINTER(c_double*three_n),
                                     POINTER(self.c_imageint*n), c_int]
    return self.lib.lammps_create_atoms(self.lmp, n, id_lmp, type_lmp, x_lmp, v_lmp, img_lmp, se_lmp)

  # -------------------------------------------------------------------------

  @property
  def has_mpi_support(self):
    """ Report whether the LAMMPS shared library was compiled with a
    real MPI library or in serial.

    This is a wrapper around the :cpp:func:`lammps_config_has_mpi_support`
    function of the library interface.

    :return: False when compiled with MPI STUBS, otherwise True
    :rtype: bool
    """
    return self.lib.lammps_config_has_mpi_support() != 0

  # -------------------------------------------------------------------------

  @property
  def is_running(self):
    """ Report whether being called from a function during a run or a minimization

    Various LAMMPS commands must not be called during an ongoing
    run or minimization.  This property allows to check for that.
    This is a wrapper around the :cpp:func:`lammps_is_running`
    function of the library interface.

    .. versionadded:: 9Oct2020

    :return: True when called during a run otherwise false
    :rtype: bool
    """
    return self.lib.lammps_is_running(self.lmp) == 1

  # -------------------------------------------------------------------------

  def force_timeout(self):
    """ Trigger an immediate timeout, i.e. a "soft stop" of a run.

    This function allows to cleanly stop an ongoing run or minimization
    at the next loop iteration.
    This is a wrapper around the :cpp:func:`lammps_force_timeout`
    function of the library interface.

    .. versionadded:: 9Oct2020
    """
    self.lib.lammps_force_timeout(self.lmp)

  # -------------------------------------------------------------------------

  @property
  def has_exceptions(self):
    """ Report whether the LAMMPS shared library was compiled with C++
    exceptions handling enabled

    This is a wrapper around the :cpp:func:`lammps_config_has_exceptions`
    function of the library interface.

    :return: state of C++ exception support
    :rtype: bool
    """
    return self.lib.lammps_config_has_exceptions() != 0

  # -------------------------------------------------------------------------

  @property
  def has_gzip_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for reading and writing compressed files through ``gzip``.

    This is a wrapper around the :cpp:func:`lammps_config_has_gzip_support`
    function of the library interface.

    :return: state of gzip support
    :rtype: bool
    """
    return self.lib.lammps_config_has_gzip_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_png_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for writing images in PNG format.

    This is a wrapper around the :cpp:func:`lammps_config_has_png_support`
    function of the library interface.

    :return: state of PNG support
    :rtype: bool
    """
    return self.lib.lammps_config_has_png_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_jpeg_support(self):
    """ Report whether the LAMMPS shared library was compiled with support
    for writing images in JPEG format.

    This is a wrapper around the :cpp:func:`lammps_config_has_jpeg_support`
    function of the library interface.

    :return: state of JPEG support
    :rtype: bool
    """
    return self.lib.lammps_config_has_jpeg_support() != 0

  # -------------------------------------------------------------------------

  @property
  def has_ffmpeg_support(self):
    """ State of support for writing movies with ``ffmpeg`` in the LAMMPS shared library

    This is a wrapper around the :cpp:func:`lammps_config_has_ffmpeg_support`
    function of the library interface.

    :return: state of ffmpeg support
    :rtype: bool
    """
    return self.lib.lammps_config_has_ffmpeg_support() != 0

  # -------------------------------------------------------------------------

  @property
  def installed_packages(self):
    """ List of the names of enabled packages in the LAMMPS shared library

    This is a wrapper around the functions :cpp:func:`lammps_config_package_count`
    and :cpp:func`lammps_config_package_name` of the library interface.

    :return
    """
    if self._installed_packages is None:
      self._installed_packages = []
      npackages = self.lib.lammps_config_package_count()
      sb = create_string_buffer(100)
      for idx in range(npackages):
        self.lib.lammps_config_package_name(idx, sb, 100)
        self._installed_packages.append(sb.value.decode())
    return self._installed_packages

  # -------------------------------------------------------------------------

  def has_style(self, category, name):
    """Returns whether a given style name is available in a given category

    This is a wrapper around the function :cpp:func:`lammps_has_style`
    of the library interface.

    :param category: name of category
    :type  category: string
    :param name: name of the style
    :type  name: string

    :return: true if style is available in given category
    :rtype:  bool
    """
    return self.lib.lammps_has_style(self.lmp, category.encode(), name.encode()) != 0

  # -------------------------------------------------------------------------

  def available_styles(self, category):
    """Returns a list of styles available for a given category

    This is a wrapper around the functions :cpp:func:`lammps_style_count()`
    and :cpp:func:`lammps_style_name()` of the library interface.

    :param category: name of category
    :type  category: string

    :return: list of style names in given category
    :rtype:  list
    """
    if self._available_styles is None:
      self._available_styles = {}

    if category not in self._available_styles:
      self._available_styles[category] = []
      nstyles = self.lib.lammps_style_count(self.lmp, category.encode())
      sb = create_string_buffer(100)
      for idx in range(nstyles):
        self.lib.lammps_style_name(self.lmp, category.encode(), idx, sb, 100)
        self._available_styles[category].append(sb.value.decode())
    return self._available_styles[category]

  # -------------------------------------------------------------------------

  def has_id(self, category, name):
    """Returns whether a given ID name is available in a given category

    This is a wrapper around the function :cpp:func:`lammps_has_id`
    of the library interface.

    .. versionadded:: 9Oct2020

    :param category: name of category
    :type  category: string
    :param name: name of the ID
    :type  name: string

    :return: true if ID is available in given category
    :rtype:  bool
    """
    return self.lib.lammps_has_id(self.lmp, category.encode(), name.encode()) != 0

  # -------------------------------------------------------------------------

  def available_ids(self, category):
    """Returns a list of IDs available for a given category

    This is a wrapper around the functions :cpp:func:`lammps_id_count()`
    and :cpp:func:`lammps_id_name()` of the library interface.

    .. versionadded:: 9Oct2020

    :param category: name of category
    :type  category: string

    :return: list of id names in given category
    :rtype:  list
    """

    categories = ['compute','dump','fix','group','molecule','region','variable']
    available_ids = []
    if category in categories:
      num = self.lib.lammps_id_count(self.lmp, category.encode())
      sb = create_string_buffer(100)
      for idx in range(num):
        self.lib.lammps_id_name(self.lmp, category.encode(), idx, sb, 100)
        available_ids.append(sb.value.decode())
    return available_ids

  # -------------------------------------------------------------------------

  def set_fix_external_callback(self, fix_name, callback, caller=None):
    import numpy as np

    def callback_wrapper(caller, ntimestep, nlocal, tag_ptr, x_ptr, fext_ptr):
      tag = self.numpy.iarray(self.c_tagint, tag_ptr, nlocal, 1)
      x   = self.numpy.darray(x_ptr, nlocal, 3)
      f   = self.numpy.darray(fext_ptr, nlocal, 3)
      callback(caller, ntimestep, nlocal, tag, x, f)

    cFunc   = self.FIX_EXTERNAL_CALLBACK_FUNC(callback_wrapper)
    cCaller = caller

    self.callback[fix_name] = { 'function': cFunc, 'caller': caller }
    self.lib.lammps_set_fix_external_callback(self.lmp, fix_name.encode(), cFunc, cCaller)


  # -------------------------------------------------------------------------

  def get_neighlist(self, idx):
    """Returns an instance of :class:`NeighList` which wraps access to the neighbor list with the given index

    See :py:meth:`lammps.numpy.get_neighlist() <lammps.numpy_wrapper.get_neighlist()>` if you want to use
    NumPy arrays instead of ``c_int`` pointers.

    :param idx: index of neighbor list
    :type  idx: int
    :return: an instance of :class:`NeighList` wrapping access to neighbor list data
    :rtype:  NeighList
    """
    if idx < 0:
        return None
    return NeighList(self, idx)

  # -------------------------------------------------------------------------

  def get_neighlist_size(self, idx):
    """Return the number of elements in neighbor list with the given index

    :param idx: neighbor list index
    :type  idx: int
    :return: number of elements in neighbor list with index idx
    :rtype:  int
     """
    return self.lib.lammps_neighlist_num_elements(self.lmp, idx)

  # -------------------------------------------------------------------------

  def get_neighlist_element_neighbors(self, idx, element):
    """Return data of neighbor list entry

    :param element: neighbor list index
    :type  element: int
    :param element: neighbor list element index
    :type  element: int
    :return: tuple with atom local index, number of neighbors and array of neighbor local atom indices
    :rtype:  (int, int, POINTER(c_int))
    """
    c_iatom = c_int()
    c_numneigh = c_int()
    c_neighbors = POINTER(c_int)()
    self.lib.lammps_neighlist_element_neighbors(self.lmp, idx, element, byref(c_iatom), byref(c_numneigh), byref(c_neighbors))
    return c_iatom.value, c_numneigh.value, c_neighbors

  # -------------------------------------------------------------------------

  def find_pair_neighlist(self, style, exact=True, nsub=0, request=0):
    """Find neighbor list index of pair style neighbor list

    Try finding pair instance that matches style. If exact is set, the pair must
    match style exactly. If exact is 0, style must only be contained. If pair is
    of style pair/hybrid, style is instead matched the nsub-th hybrid sub-style.

    Once the pair instance has been identified, multiple neighbor list requests
    may be found. Every neighbor list is uniquely identified by its request
    index. Thus, providing this request index ensures that the correct neighbor
    list index is returned.

    :param style: name of pair style that should be searched for
    :type  style: string
    :param exact: controls whether style should match exactly or only must be contained in pair style name, defaults to True
    :type  exact: bool, optional
    :param nsub:  match nsub-th hybrid sub-style, defaults to 0
    :type  nsub:  int, optional
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    style = style.encode()
    exact = int(exact)
    idx = self.lib.lammps_find_pair_neighlist(self.lmp, style, exact, nsub, request)
    return idx

  # -------------------------------------------------------------------------

  def find_fix_neighlist(self, fixid, request=0):
    """Find neighbor list index of fix neighbor list

    :param fixid: name of fix
    :type  fixid: string
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    fixid = fixid.encode()
    idx = self.lib.lammps_find_fix_neighlist(self.lmp, fixid, request)
    return idx

  # -------------------------------------------------------------------------

  def find_compute_neighlist(self, computeid, request=0):
    """Find neighbor list index of compute neighbor list

    :param computeid: name of compute
    :type  computeid: string
    :param request:   index of neighbor list request, in case there are more than one, defaults to 0
    :type  request:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int
     """
    computeid = computeid.encode()
    idx = self.lib.lammps_find_compute_neighlist(self.lmp, computeid, request)
    return idx

# -------------------------------------------------------------------------

class numpy_wrapper:
  """lammps API NumPy Wrapper

  This is a wrapper class that provides additional methods on top of an
  existing :py:class:`lammps` instance. The methods transform raw ctypes
  pointers into NumPy arrays, which give direct access to the
  original data while protecting against out-of-bounds accesses.

  There is no need to explicitly instantiate this class. Each instance
  of :py:class:`lammps` has a :py:attr:`numpy <lammps.numpy>` property
  that returns an instance.

  :param lmp: instance of the :py:class:`lammps` class
  :type  lmp: lammps
  """
  def __init__(self, lmp):
    self.lmp = lmp

  # -------------------------------------------------------------------------

  def _ctype_to_numpy_int(self, ctype_int):
    import numpy as np
    if ctype_int == c_int32:
      return np.int32
    elif ctype_int == c_int64:
      return np.int64
    return np.intc

  # -------------------------------------------------------------------------

  def extract_atom(self, name, dtype=LAMMPS_AUTODETECT, nelem=LAMMPS_AUTODETECT, dim=LAMMPS_AUTODETECT):
    """Retrieve per-atom properties from LAMMPS as NumPy arrays

    This is a wrapper around the :py:meth:`lammps.extract_atom()` method.
    It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

    .. note::

       While the returned arrays of per-atom data are dimensioned
       for the range [0:nmax] - as is the underlying storage -
       the data is usually only valid for the range of [0:nlocal],
       unless the property of interest is also updated for ghost
       atoms.  In some cases, this depends on a LAMMPS setting, see
       for example :doc:`comm_modify vel yes <comm_modify>`.

    :param name: name of the property
    :type name:  string
    :param dtype: type of the returned data (see :ref:`py_datatype_constants`)
    :type dtype:  int, optional
    :param nelem: number of elements in array
    :type nelem:  int, optional
    :param dim: dimension of each element
    :type dim:  int, optional
    :return: requested data as NumPy array with direct access to C data or None
    :rtype: numpy.array or NoneType
    """
    if dtype == LAMMPS_AUTODETECT:
      dtype = self.lmp.extract_atom_datatype(name)

    if nelem == LAMMPS_AUTODETECT:
      if name == "mass":
        nelem = self.lmp.extract_global("ntypes") + 1
      else:
        nelem = self.lmp.extract_global("nlocal")
    if dim == LAMMPS_AUTODETECT:
      if dtype in (LAMMPS_INT_2D, LAMMPS_DOUBLE_2D, LAMMPS_INT64_2D):
        # TODO add other fields
        if name in ("x", "v", "f", "angmom", "torque", "csforce", "vforce"):
          dim = 3
        else:
          dim = 2
      else:
        dim = 1

    raw_ptr = self.lmp.extract_atom(name, dtype)

    if dtype in (LAMMPS_DOUBLE, LAMMPS_DOUBLE_2D):
      return self.darray(raw_ptr, nelem, dim)
    elif dtype in (LAMMPS_INT, LAMMPS_INT_2D):
      return self.iarray(c_int32, raw_ptr, nelem, dim)
    elif dtype in (LAMMPS_INT64, LAMMPS_INT64_2D):
      return self.iarray(c_int64, raw_ptr, nelem, dim)
    return raw_ptr

  # -------------------------------------------------------------------------

  def extract_atom_iarray(self, name, nelem, dim=1):
    warnings.warn("deprecated, use extract_atom instead", DeprecationWarning)

    if name in ['id', 'molecule']:
      c_int_type = self.lmp.c_tagint
    elif name in ['image']:
      c_int_type = self.lmp.c_imageint
    else:
      c_int_type = c_int

    if dim == 1:
      raw_ptr = self.lmp.extract_atom(name, LAMMPS_INT)
    else:
      raw_ptr = self.lmp.extract_atom(name, LAMMPS_INT_2D)

    return self.iarray(c_int_type, raw_ptr, nelem, dim)

  # -------------------------------------------------------------------------

  def extract_atom_darray(self, name, nelem, dim=1):
    warnings.warn("deprecated, use extract_atom instead", DeprecationWarning)

    if dim == 1:
      raw_ptr = self.lmp.extract_atom(name, LAMMPS_DOUBLE)
    else:
      raw_ptr = self.lmp.extract_atom(name, LAMMPS_DOUBLE_2D)

    return self.darray(raw_ptr, nelem, dim)

  # -------------------------------------------------------------------------

  def extract_compute(self, cid, style, type):
    """Retrieve data from a LAMMPS compute

    This is a wrapper around the
    :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>` method.
    It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

    :param id: compute ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type style:  int
    :param type: type of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type type:  int
    :return: requested data either as float, as NumPy array with direct access to C data, or None
    :rtype: float, numpy.array, or NoneType
    """
    value = self.lmp.extract_compute(cid, style, type)

    if style in (LMP_STYLE_GLOBAL, LMP_STYLE_LOCAL):
      if type == LMP_TYPE_VECTOR:
        nrows = self.lmp.extract_compute(cid, style, LMP_SIZE_VECTOR)
        return self.darray(value, nrows)
      elif type == LMP_TYPE_ARRAY:
        nrows = self.lmp.extract_compute(cid, style, LMP_SIZE_ROWS)
        ncols = self.lmp.extract_compute(cid, style, LMP_SIZE_COLS)
        return self.darray(value, nrows, ncols)
    elif style == LMP_STYLE_ATOM:
      if type == LMP_TYPE_VECTOR:
        nlocal = self.lmp.extract_global("nlocal")
        return self.darray(value, nlocal)
      elif type == LMP_TYPE_ARRAY:
        nlocal = self.lmp.extract_global("nlocal")
        ncols = self.lmp.extract_compute(cid, style, LMP_SIZE_COLS)
        return self.darray(value, nlocal, ncols)
    return value

  # -------------------------------------------------------------------------

  def extract_fix(self, fid, style, type, nrow=0, ncol=0):
    """Retrieve data from a LAMMPS fix

    This is a wrapper around the :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>` method.
    It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

    :param id: fix ID
    :type id:  string
    :param style: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type style:  int
    :param type: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type type:  int
    :param nrow: index of global vector element or row index of global array element
    :type nrow:  int
    :param ncol: column index of global array element
    :type ncol:  int
    :return: requested data
    :rtype: integer or double value, pointer to 1d or 2d double array  or None

    """
    value = self.lmp.extract_fix(fid, style, type, nrow, ncol)
    if style == LMP_STYLE_ATOM:
      if type == LMP_TYPE_VECTOR:
        nlocal = self.lmp.extract_global("nlocal")
        return self.darray(value, nlocal)
      elif type == LMP_TYPE_ARRAY:
        nlocal = self.lmp.extract_global("nlocal")
        ncols = self.lmp.extract_fix(fid, style, LMP_SIZE_COLS, 0, 0)
        return self.darray(value, nlocal, ncols)
    elif style == LMP_STYLE_LOCAL:
      if type == LMP_TYPE_VECTOR:
        nrows = self.lmp.extract_fix(fid, style, LMP_SIZE_ROWS, 0, 0)
        return self.darray(value, nrows)
      elif type == LMP_TYPE_ARRAY:
        nrows = self.lmp.extract_fix(fid, style, LMP_SIZE_ROWS, 0, 0)
        ncols = self.lmp.extract_fix(fid, style, LMP_SIZE_COLS, 0, 0)
        return self.darray(value, nrows, ncols)
    return value

  # -------------------------------------------------------------------------

  def extract_variable(self, name, group=None, vartype=LMP_VAR_EQUAL):
    """ Evaluate a LAMMPS variable and return its data

    This function is a wrapper around the function
    :py:meth:`lammps.extract_variable() <lammps.lammps.extract_variable()>`
    method. It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

    :param name: name of the variable to execute
    :type name: string
    :param group: name of group for atom-style variable (ignored for equal-style variables)
    :type group: string
    :param vartype: type of variable, see :ref:`py_vartype_constants`
    :type vartype: int
    :return: the requested data or None
    :rtype: c_double, numpy.array, or NoneType
    """
    import numpy as np
    value = self.lmp.extract_variable(name, group, vartype)
    if vartype == LMP_VAR_ATOM:
      return np.ctypeslib.as_array(value)
    return value

  # -------------------------------------------------------------------------

  def get_neighlist(self, idx):
    """Returns an instance of :class:`NumPyNeighList` which wraps access to the neighbor list with the given index

    :param idx: index of neighbor list
    :type  idx: int
    :return: an instance of :class:`NumPyNeighList` wrapping access to neighbor list data
    :rtype:  NumPyNeighList
    """
    if idx < 0:
        return None
    return NumPyNeighList(self.lmp, idx)

  # -------------------------------------------------------------------------

  def get_neighlist_element_neighbors(self, idx, element):
    """Return data of neighbor list entry

    This function is a wrapper around the function
    :py:meth:`lammps.get_neighlist_element_neighbors() <lammps.lammps.get_neighlist_element_neighbors()>`
    method. It behaves the same as the original method, but returns a NumPy array containing the neighbors
    instead of a ``ctypes`` pointer.

    :param element: neighbor list index
    :type  element: int
    :param element: neighbor list element index
    :type  element: int
    :return: tuple with atom local index and numpy array of neighbor local atom indices
    :rtype:  (int, numpy.array)
    """
    iatom, numneigh, c_neighbors = self.lmp.get_neighlist_element_neighbors(idx, element)
    neighbors = self.iarray(c_int, c_neighbors, numneigh, 1)
    return iatom, neighbors

  # -------------------------------------------------------------------------

  def iarray(self, c_int_type, raw_ptr, nelem, dim=1):
    import numpy as np
    np_int_type = self._ctype_to_numpy_int(c_int_type)

    if dim == 1:
      ptr = cast(raw_ptr, POINTER(c_int_type * nelem))
    else:
      ptr = cast(raw_ptr[0], POINTER(c_int_type * nelem * dim))

    a = np.frombuffer(ptr.contents, dtype=np_int_type)
    a.shape = (nelem, dim)
    return a

  # -------------------------------------------------------------------------

  def darray(self, raw_ptr, nelem, dim=1):
    import numpy as np
    if dim == 1:
      ptr = cast(raw_ptr, POINTER(c_double * nelem))
    else:
      ptr = cast(raw_ptr[0], POINTER(c_double * nelem * dim))

    a = np.frombuffer(ptr.contents)
    a.shape = (nelem, dim)
    return a


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

################################################################################
# Alternative Python Wrapper
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

class OutputCapture(object):
  """ Utility class to capture LAMMPS library output """

  def __init__(self):
    self.stdout_pipe_read, self.stdout_pipe_write = os.pipe()
    self.stdout_fd = 1

  def __enter__(self):
    self.stdout = os.dup(self.stdout_fd)
    os.dup2(self.stdout_pipe_write, self.stdout_fd)
    return self

  def __exit__(self, type, value, tracebac):
    os.dup2(self.stdout, self.stdout_fd)
    os.close(self.stdout)
    os.close(self.stdout_pipe_read)
    os.close(self.stdout_pipe_write)

  # check if we have more to read from the pipe
  def more_data(self, pipe):
    r, _, _ = select.select([pipe], [], [], 0)
    return bool(r)

  # read the whole pipe
  def read_pipe(self, pipe):
    out = ""
    while self.more_data(pipe):
      out += os.read(pipe, 1024).decode()
    return out

  @property
  def output(self):
    return self.read_pipe(self.stdout_pipe_read)

# -------------------------------------------------------------------------

class Variable(object):
  def __init__(self, pylammps_instance, name, style, definition):
    self._pylmp = pylammps_instance
    self.name = name
    self.style = style
    self.definition = definition.split()

  @property
  def value(self):
    if self.style == 'atom':
      return list(self._pylmp.lmp.extract_variable(self.name, "all", 1))
    else:
      value = self._pylmp.lmp_print('"${%s}"' % self.name).strip()
      try:
        return float(value)
      except ValueError:
        return value

# -------------------------------------------------------------------------

class AtomList(object):
  """
  A dynamic list of atoms that returns either an :py:class:`Atom` or
  :py:class:`Atom2D` instance for each atom. Instances are only allocated
  when accessed.

  :ivar natoms: total number of atoms
  :ivar dimensions: number of dimensions in system
  """
  def __init__(self, pylammps_instance):
    self._pylmp = pylammps_instance
    self.natoms = self._pylmp.system.natoms
    self.dimensions = self._pylmp.system.dimensions
    self._loaded = {}

  def __getitem__(self, index):
    """
    Return Atom with given local index

    :param index: Local index of atom
    :type index: int
    :rtype: Atom or Atom2D
    """
    if index not in self._loaded:
        if self.dimensions == 2:
            atom = Atom2D(self._pylmp, index + 1)
        else:
            atom = Atom(self._pylmp, index + 1)
        self._loaded[index] = atom
    return self._loaded[index]

  def __len__(self):
    return self.natoms


# -------------------------------------------------------------------------

class Atom(object):
  """
  A wrapper class then represents a single atom inside of LAMMPS

  It provides access to properties of the atom and allows you to change some of them.
  """
  def __init__(self, pylammps_instance, index):
    self._pylmp = pylammps_instance
    self.index = index

  @property
  def id(self):
    """
    Return the atom ID

    :type: int
    """
    return int(self._pylmp.eval("id[%d]" % self.index))

  @property
  def type(self):
    """
    Return the atom type

    :type: int
    """
    return int(self._pylmp.eval("type[%d]" % self.index))

  @property
  def mol(self):
    """
    Return the atom molecule index

    :type: int
    """
    return self._pylmp.eval("mol[%d]" % self.index)

  @property
  def mass(self):
    """
    Return the atom mass

    :type: float
    """
    return self._pylmp.eval("mass[%d]" % self.index)

  @property
  def position(self):
    """
    :getter: Return position of atom
    :setter: Set position of atom
    :type: tuple (float, float, float)
    """
    return (self._pylmp.eval("x[%d]" % self.index),
            self._pylmp.eval("y[%d]" % self.index),
            self._pylmp.eval("z[%d]" % self.index))

  @position.setter
  def position(self, value):
    """
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: tuple (float, float, float)
    """
    self._pylmp.set("atom", self.index, "x", value[0])
    self._pylmp.set("atom", self.index, "y", value[1])
    self._pylmp.set("atom", self.index, "z", value[2])

  @property
  def velocity(self):
    return (self._pylmp.eval("vx[%d]" % self.index),
            self._pylmp.eval("vy[%d]" % self.index),
            self._pylmp.eval("vz[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self._pylmp.set("atom", self.index, "vx", value[0])
     self._pylmp.set("atom", self.index, "vy", value[1])
     self._pylmp.set("atom", self.index, "vz", value[2])

  @property
  def force(self):
    """
    Return the total force acting on the atom

    :type: tuple (float, float, float)
    """
    return (self._pylmp.eval("fx[%d]" % self.index),
            self._pylmp.eval("fy[%d]" % self.index),
            self._pylmp.eval("fz[%d]" % self.index))

  @property
  def charge(self):
    """
    Return the atom charge

    :type: float
    """
    return self._pylmp.eval("q[%d]" % self.index)

# -------------------------------------------------------------------------

class Atom2D(Atom):
  """
  A wrapper class then represents a single 2D atom inside of LAMMPS

  Inherits all properties from the :py:class:`Atom` class, but returns 2D versions
  of position, velocity, and force.

  It provides access to properties of the atom and allows you to change some of them.
  """
  def __init__(self, pylammps_instance, index):
    super(Atom2D, self).__init__(pylammps_instance, index)

  @property
  def position(self):
    """
    :getter: Return position of atom
    :setter: Set position of atom
    :type: tuple (float, float)
    """
    return (self._pylmp.eval("x[%d]" % self.index),
            self._pylmp.eval("y[%d]" % self.index))

  @position.setter
  def position(self, value):
     self._pylmp.set("atom", self.index, "x", value[0])
     self._pylmp.set("atom", self.index, "y", value[1])

  @property
  def velocity(self):
    """
    :getter: Return velocity of atom
    :setter: Set velocity of atom
    :type: tuple (float, float)
    """
    return (self._pylmp.eval("vx[%d]" % self.index),
            self._pylmp.eval("vy[%d]" % self.index))

  @velocity.setter
  def velocity(self, value):
     self._pylmp.set("atom", self.index, "vx", value[0])
     self._pylmp.set("atom", self.index, "vy", value[1])

  @property
  def force(self):
    """
    Return the total force acting on the atom

    :type: tuple (float, float)
    """
    return (self._pylmp.eval("fx[%d]" % self.index),
            self._pylmp.eval("fy[%d]" % self.index))

# -------------------------------------------------------------------------

class variable_set:
    def __init__(self, name, variable_dict):
        self._name = name
        array_pattern = re.compile(r"(?P<arr>.+)\[(?P<index>[0-9]+)\]")

        for key, value in variable_dict.items():
            m = array_pattern.match(key)
            if m:
                g = m.groupdict()
                varname = g['arr']
                idx = int(g['index'])
                if varname not in self.__dict__:
                    self.__dict__[varname] = {}
                self.__dict__[varname][idx] = value
            else:
                self.__dict__[key] = value

    def __str__(self):
        return "{}({})".format(self._name, ','.join(["{}={}".format(k, self.__dict__[k]) for k in self.__dict__.keys() if not k.startswith('_')]))

    def __repr__(self):
        return self.__str__()

# -------------------------------------------------------------------------

def get_thermo_data(output):
    """ traverse output of runs and extract thermo data columns """
    if isinstance(output, str):
        lines = output.splitlines()
    else:
        lines = output

    runs = []
    columns = []
    in_run = False
    current_run = {}

    for line in lines:
        if line.startswith("Per MPI rank memory allocation"):
            in_run = True
        elif in_run and len(columns) == 0:
            # first line after memory usage are column names
            columns = line.split()

            current_run = {}

            for col in columns:
                current_run[col] = []

        elif line.startswith("Loop time of "):
            in_run = False
            columns = None
            thermo_data = variable_set('ThermoData', current_run)
            r = {'thermo' : thermo_data }
            runs.append(namedtuple('Run', list(r.keys()))(*list(r.values())))
        elif in_run and len(columns) > 0:
            items = line.split()
            # Convert thermo output and store it.
            # It must have the same number of columns and
            # all of them must be convertible to floats.
            # Otherwise we ignore the line
            if len(items) == len(columns):
                try:
                    values = [float(x) for x in items]
                    for i, col in enumerate(columns):
                        current_run[col].append(values[i])
                except ValueError:
                  pass

    return runs

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

class PyLammps(object):
  """
  This is a Python wrapper class around the lower-level
  :py:class:`lammps` class, exposing a more Python-like,
  object-oriented interface for prototyping system inside of IPython and
  Jupyter notebooks.

  It either creates its own instance of :py:class:`lammps` or can be
  initialized with an existing instance. The arguments are the same of the
  lower-level interface. The original interface can still be accessed via
  :py:attr:`PyLammps.lmp`.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm

  :ivar lmp:  instance of original LAMMPS Python interface
  :vartype lmp: :py:class:`lammps`

  :ivar runs:  list of completed runs, each storing the thermo output
  :vartype run: list
  """

  def __init__(self, name="", cmdargs=None, ptr=None, comm=None):
    self.has_echo = False

    if cmdargs:
      if '-echo' in cmdargs:
        idx = cmdargs.index('-echo')
        # ensures that echo line is ignored during output capture
        self.has_echo = idx+1 < len(cmdargs) and cmdargs[idx+1] in ('screen', 'both')

    if ptr:
      if isinstance(ptr,PyLammps):
        self.lmp = ptr.lmp
      elif isinstance(ptr,lammps):
        self.lmp = ptr
      else:
        self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)
    else:
      self.lmp = lammps(name=name,cmdargs=cmdargs,ptr=None,comm=comm)
    print("LAMMPS output is captured by PyLammps wrapper")
    self._cmd_history = []
    self.runs = []

  def __del__(self):
    if self.lmp: self.lmp.close()
    self.lmp = None

  def close(self):
    """Explicitly delete a LAMMPS instance

    This is a wrapper around the :py:meth:`lammps.close` of the Python interface.
    """
    if self.lmp: self.lmp.close()
    self.lmp = None

  def version(self):
    """Return a numerical representation of the LAMMPS version in use.

    This is a wrapper around the :py:meth:`lammps.version` function of the Python interface.

    :return: version number
    :rtype:  int
    """
    return self.lmp.version()

  def file(self, file):
    """Read LAMMPS commands from a file.

    This is a wrapper around the :py:meth:`lammps.file` function of the Python interface.

    :param path: Name of the file/path with LAMMPS commands
    :type path:  string
    """
    self.lmp.file(file)

  def write_script(self, filepath):
    """
    Write LAMMPS script file containing all commands executed up until now

    :param filepath: path to script file that should be written
    :type filepath: string
    """
    with open(filepath, "w") as f:
      for cmd in self._cmd_history:
        print(cmd, file=f)

  def command(self, cmd):
    """
    Execute LAMMPS command

    All commands executed will be stored in a command history which can be
    written to a file using :py:meth:`PyLammps.write_script()`

    :param cmd: command string that should be executed
    :type: cmd: string
    """
    self.lmp.command(cmd)
    self._cmd_history.append(cmd)

  def run(self, *args, **kwargs):
    """
    Execute LAMMPS run command with given arguments

    All thermo output during the run is captured and saved as new entry in
    :py:attr:`PyLammps.runs`. The latest run can be retrieved by
    :py:attr:`PyLammps.last_run`.
    """
    output = self.__getattr__('run')(*args, **kwargs)

    comm = self.lmp.get_mpi_comm()
    if comm:
      output = self.lmp.comm.bcast(output, root=0)

    self.runs += get_thermo_data(output)
    return output

  @property
  def last_run(self):
    """
    Return data produced of last completed run command

    :getter: Returns an object containing information about the last run command
    :type: dict
    """
    if len(self.runs) > 0:
        return self.runs[-1]
    return None

  @property
  def atoms(self):
    """
    All atoms of this LAMMPS instance

    :getter: Returns a list of atoms currently in the system
    :type: AtomList
    """
    return AtomList(self)

  @property
  def system(self):
    """
    The system state of this LAMMPS instance

    :getter: Returns an object with properties storing the current system state
    :type: namedtuple
    """
    output = self.info("system")
    d = self._parse_info_system(output)
    return namedtuple('System', d.keys())(*d.values())

  @property
  def communication(self):
    """
    The communication state of this LAMMPS instance

    :getter: Returns an object with properties storing the current communication state
    :type: namedtuple
    """
    output = self.info("communication")
    d = self._parse_info_communication(output)
    return namedtuple('Communication', d.keys())(*d.values())

  @property
  def computes(self):
    """
    The list of active computes of this LAMMPS instance

    :getter: Returns a list of computes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("computes")
    return self._parse_element_list(output)

  @property
  def dumps(self):
    """
    The list of active dumps of this LAMMPS instance

    :getter: Returns a list of dumps that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("dumps")
    return self._parse_element_list(output)

  @property
  def fixes(self):
    """
    The list of active fixes of this LAMMPS instance

    :getter: Returns a list of fixes that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("fixes")
    return self._parse_element_list(output)

  @property
  def groups(self):
    """
    The list of active atom groups of this LAMMPS instance

    :getter: Returns a list of atom groups that are currently active in this LAMMPS instance
    :type: list
    """
    output = self.info("groups")
    return self._parse_groups(output)

  @property
  def variables(self):
    """
    Returns a dictionary of all variables defined in the current LAMMPS instance

    :getter: Returns a dictionary of all variables that are defined in this LAMMPS instance
    :type: dict
    """
    output = self.info("variables")
    vars = {}
    for v in self._parse_element_list(output):
      vars[v['name']] = Variable(self, v['name'], v['style'], v['def'])
    return vars

  def eval(self, expr):
    """
    Evaluate expression

    :param expr: the expression string that should be evaluated inside of LAMMPS
    :type expr: string

    :return: the value of the evaluated expression
    :rtype: float if numeric, string otherwise
    """
    value = self.lmp_print('"$(%s)"' % expr).strip()
    try:
      return float(value)
    except ValueError:
      return value

  def _split_values(self, line):
    return [x.strip() for x in line.split(',')]

  def _get_pair(self, value):
    return [x.strip() for x in value.split('=')]

  def _parse_info_system(self, output):
    lines = output[6:-2]
    system = {}

    for line in lines:
      if line.startswith("Units"):
        system['units'] = self._get_pair(line)[1]
      elif line.startswith("Atom style"):
        system['atom_style'] = self._get_pair(line)[1]
      elif line.startswith("Atom map"):
        system['atom_map'] = self._get_pair(line)[1]
      elif line.startswith("Atoms"):
        parts = self._split_values(line)
        system['natoms'] = int(self._get_pair(parts[0])[1])
        system['ntypes'] = int(self._get_pair(parts[1])[1])
        system['style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Kspace style"):
        system['kspace_style'] = self._get_pair(line)[1]
      elif line.startswith("Dimensions"):
        system['dimensions'] = int(self._get_pair(line)[1])
      elif line.startswith("Orthogonal box"):
        system['orthogonal_box'] = [float(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Boundaries"):
        system['boundaries'] = self._get_pair(line)[1]
      elif line.startswith("xlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("ylo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("zlo"):
        keys, values = [self._split_values(x) for x in self._get_pair(line)]
        for key, value in zip(keys, values):
          system[key] = float(value)
      elif line.startswith("Molecule type"):
        system['molecule_type'] = self._get_pair(line)[1]
      elif line.startswith("Bonds"):
        parts = self._split_values(line)
        system['nbonds'] = int(self._get_pair(parts[0])[1])
        system['nbondtypes'] = int(self._get_pair(parts[1])[1])
        system['bond_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Angles"):
        parts = self._split_values(line)
        system['nangles'] = int(self._get_pair(parts[0])[1])
        system['nangletypes'] = int(self._get_pair(parts[1])[1])
        system['angle_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Dihedrals"):
        parts = self._split_values(line)
        system['ndihedrals'] = int(self._get_pair(parts[0])[1])
        system['ndihedraltypes'] = int(self._get_pair(parts[1])[1])
        system['dihedral_style'] = self._get_pair(parts[2])[1]
      elif line.startswith("Impropers"):
        parts = self._split_values(line)
        system['nimpropers'] = int(self._get_pair(parts[0])[1])
        system['nimpropertypes'] = int(self._get_pair(parts[1])[1])
        system['improper_style'] = self._get_pair(parts[2])[1]

    return system

  def _parse_info_communication(self, output):
    lines = output[6:-3]
    comm = {}

    for line in lines:
      if line.startswith("MPI library"):
        comm['mpi_version'] = line.split(':')[1].strip()
      elif line.startswith("Comm style"):
        parts = self._split_values(line)
        comm['comm_style'] = self._get_pair(parts[0])[1]
        comm['comm_layout'] = self._get_pair(parts[1])[1]
      elif line.startswith("Processor grid"):
        comm['proc_grid'] = [int(x) for x in self._get_pair(line)[1].split('x')]
      elif line.startswith("Communicate velocities for ghost atoms"):
        comm['ghost_velocity'] = (self._get_pair(line)[1] == "yes")
      elif line.startswith("Nprocs"):
        parts = self._split_values(line)
        comm['nprocs'] = int(self._get_pair(parts[0])[1])
        comm['nthreads'] = int(self._get_pair(parts[1])[1])
    return comm

  def _parse_element_list(self, output):
    lines = output[6:-3]
    elements = []

    for line in lines:
      element_info = self._split_values(line.split(':')[1].strip())
      element = {'name': element_info[0]}
      for key, value in [self._get_pair(x) for x in element_info[1:]]:
        element[key] = value
      elements.append(element)
    return elements

  def _parse_groups(self, output):
    lines = output[6:-3]
    groups = []
    group_pattern = re.compile(r"(?P<name>.+) \((?P<type>.+)\)")

    for line in lines:
      m = group_pattern.match(line.split(':')[1].strip())
      group = {'name': m.group('name'), 'type': m.group('type')}
      groups.append(group)
    return groups

  def lmp_print(self, s):
    """ needed for Python2 compatibility, since print is a reserved keyword """
    return self.__getattr__("print")(s)

  def __dir__(self):
    return ['angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'atom_style',
    'bond_coeff', 'bond_style', 'boundary', 'change_box', 'communicate', 'compute',
    'create_atoms', 'create_box', 'delete_atoms', 'delete_bonds', 'dielectric',
    'dihedral_coeff', 'dihedral_style', 'dimension', 'dump', 'fix', 'fix_modify',
    'group', 'improper_coeff', 'improper_style', 'include', 'kspace_modify',
    'kspace_style', 'lattice', 'mass', 'minimize', 'min_style', 'neighbor',
    'neigh_modify', 'newton', 'nthreads', 'pair_coeff', 'pair_modify',
    'pair_style', 'processors', 'read', 'read_data', 'read_restart', 'region',
    'replicate', 'reset_timestep', 'restart', 'run', 'run_style', 'thermo',
    'thermo_modify', 'thermo_style', 'timestep', 'undump', 'unfix', 'units',
    'variable', 'velocity', 'write_restart']

  def __getattr__(self, name):
    """
    This method is where the Python 'magic' happens. If a method is not
    defined by the class PyLammps, it assumes it is a LAMMPS command. It takes
    all the arguments, concatinates them to a single string, and executes it using
    :py:meth:`lammps.PyLammps.command()`.

    :param verbose: Print output of command
    :type verbose:  bool
    :return: line or list of lines of output, None if no output
    :rtype: list or string
    """
    def handler(*args, **kwargs):
      cmd_args = [name] + [str(x) for x in args]

      with OutputCapture() as capture:
        cmd = ' '.join(cmd_args)
        self.command(cmd)
        output = capture.output

      if 'verbose' in kwargs and kwargs['verbose']:
        print(output)

      lines = output.splitlines()

      if self.has_echo:
        lines = lines[1:]

      if len(lines) > 1:
        return lines
      elif len(lines) == 1:
        return lines[0]
      return None

    return handler


class IPyLammps(PyLammps):
  """
  IPython wrapper for LAMMPS which adds embedded graphics capabilities to PyLammmps interface

  It either creates its own instance of :py:class:`lammps` or can be
  initialized with an existing instance. The arguments are the same of the
  lower-level interface. The original interface can still be accessed via
  :py:attr:`PyLammps.lmp`.

  :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
  :type  name: string
  :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
  :type  cmdargs: list
  :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
  :type  ptr: pointer
  :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
  :type  comm: MPI_Comm
  """

  def __init__(self,name="",cmdargs=None,ptr=None,comm=None):
    super(IPyLammps, self).__init__(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)

  def image(self, filename="snapshot.png", group="all", color="type", diameter="type",
            size=None, view=None, center=None, up=None, zoom=1.0, background_color="white"):
    """ Generate image using write_dump command and display it

    See :doc:`dump image <dump_image>` for more information.

    :param filename: Name of the image file that should be generated. The extension determines whether it is PNG or JPEG
    :type filename: string
    :param group: the group of atoms write_image should use
    :type group: string
    :param color: name of property used to determine color
    :type color: string
    :param diameter: name of property used to determine atom diameter
    :type diameter: string
    :param size: dimensions of image
    :type size: tuple (width, height)
    :param view: view parameters
    :type view: tuple (theta, phi)
    :param center: center parameters
    :type center: tuple (flag, center_x, center_y, center_z)
    :param up: vector pointing to up direction
    :type up: tuple (up_x, up_y, up_z)
    :param zoom: zoom factor
    :type zoom: float
    :param background_color: background color of scene
    :type background_color: string

    :return: Image instance used to display image in notebook
    :rtype: :py:class:`IPython.core.display.Image`
    """
    cmd_args = [group, "image", filename, color, diameter]

    if size:
      width = size[0]
      height = size[1]
      cmd_args += ["size", width, height]

    if view:
      theta = view[0]
      phi = view[1]
      cmd_args += ["view", theta, phi]

    if center:
      flag = center[0]
      Cx = center[1]
      Cy = center[2]
      Cz = center[3]
      cmd_args += ["center", flag, Cx, Cy, Cz]

    if up:
      Ux = up[0]
      Uy = up[1]
      Uz = up[2]
      cmd_args += ["up", Ux, Uy, Uz]

    if zoom:
      cmd_args += ["zoom", zoom]

    cmd_args.append("modify backcolor " + background_color)

    self.write_dump(*cmd_args)
    from IPython.core.display import Image
    return Image(filename)

  def video(self, filename):
    """
    Load video from file

    Can be used to visualize videos from :doc:`dump movie <dump_image>`.

    :param filename: Path to video file
    :type filename: string
    :return: HTML Video Tag used by notebook to embed a video
    :rtype: :py:class:`IPython.display.HTML`
    """
    from IPython.display import HTML
    return HTML("<video controls><source src=\"" + filename + "\"></video>")
