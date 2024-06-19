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
# Python wrapper for the LAMMPS library via ctypes

# for python2/3 compatibility

from __future__ import print_function

import os
import sys
from ctypes import *                    # lgtm [py/polluting-import]
from os.path import dirname,abspath,join
from inspect import getsourcefile

from .constants import *                # lgtm [py/polluting-import]
from .data import *                     # lgtm [py/polluting-import]

# -------------------------------------------------------------------------

class MPIAbortException(Exception):
  def __init__(self, message):
    self.message = message

  def __str__(self):
    return repr(self.message)

# -------------------------------------------------------------------------

class ExceptionCheck:
  """Utility class to rethrow LAMMPS C++ exceptions as Python exceptions"""
  def __init__(self, lmp):
    self.lmp = lmp

  def __enter__(self):
    pass

  def __exit__(self, exc_type, exc_value, traceback):
    if self.lmp.has_exceptions and self.lmp.lib.lammps_has_error(self.lmp.lmp):
      raise self.lmp._lammps_exception

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
    winpath = abspath(os.path.join(modpath,'..','..','bin'))
    # allow override for running tests on Windows
    if (os.environ.get("LAMMPSDLLPATH")):
      winpath = os.environ.get("LAMMPSDLLPATH")
    self.lib = None
    self.lmp = None

    # if a pointer to a LAMMPS object is handed in
    # when being called from a Python interpreter
    # embedded into a LAMMPS executable, all library
    # symbols should already be available so we do not
    # load a shared object.

    try:
      if ptr is not None: self.lib = CDLL("",RTLD_GLOBAL)
    except OSError:
      self.lib = None

    # load liblammps.so unless name is given
    #   if name = "g++", load liblammps_g++.so
    # try loading the LAMMPS shared object from the location
    #   of the lammps package with an absolute path,
    #   so that LD_LIBRARY_PATH does not need to be set for regular install
    # fall back to loading with a relative path,
    #   typically requires LD_LIBRARY_PATH to be set appropriately
    # guess shared library extension based on OS, if not inferred from actual file

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
    elif any([f.startswith('liblammps') and f.endswith('.so')
              for f in os.listdir(modpath)]):
      lib_ext = ".so"
    else:
      import platform
      if platform.system() == "Darwin":
        lib_ext = ".dylib"
      elif platform.system() == "Windows":
        lib_ext = ".dll"
      else:
        lib_ext = ".so"

    if not self.lib:
      if name:
        libpath = join(modpath,"liblammps_%s" % name + lib_ext)
      else:
        libpath = join(modpath,"liblammps" + lib_ext)
      if not os.path.isfile(libpath):
        if name:
          libpath = "liblammps_%s" % name + lib_ext
        else:
          libpath = "liblammps" + lib_ext
      self.lib = CDLL(libpath,RTLD_GLOBAL)

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
    self.lib.lammps_flush_buffers.argtypes = [c_void_p]
    self.lib.lammps_free.argtypes = [c_void_p]

    self.lib.lammps_error.argtypes = [c_void_p, c_int, c_char_p]

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

    self.lib.lammps_gather_atoms.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms.restype = None

    self.lib.lammps_gather_atoms_concat.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_atoms_concat.restype = None

    self.lib.lammps_gather_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_atoms_subset.restype = None

    self.lib.lammps_scatter_atoms.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_scatter_atoms.restype = None

    self.lib.lammps_scatter_atoms_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_scatter_atoms_subset.restype = None

    self.lib.lammps_gather_bonds.argtypes = [c_void_p,c_void_p]
    self.lib.lammps_gather_bonds.restype = None

    self.lib.lammps_gather_angles.argtypes = [c_void_p,c_void_p]
    self.lib.lammps_gather_angles.restype = None

    self.lib.lammps_gather_dihedrals.argtypes = [c_void_p,c_void_p]
    self.lib.lammps_gather_dihedrals.restype = None

    self.lib.lammps_gather_impropers.argtypes = [c_void_p,c_void_p]
    self.lib.lammps_gather_impropers.restype = None

    self.lib.lammps_gather.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather.restype = None

    self.lib.lammps_gather_concat.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
    self.lib.lammps_gather_concat.restype = None

    self.lib.lammps_gather_subset.argtypes = \
      [c_void_p,c_char_p,c_int,c_int,c_int,POINTER(c_int),c_void_p]
    self.lib.lammps_gather_subset.restype = None

    self.lib.lammps_scatter.argtypes = [c_void_p,c_char_p,c_int,c_int,c_void_p]
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

    self.lib.lammps_neighlist_element_neighbors.argtypes = \
      [c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int))]
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

    self.lib.lammps_map_atom.argtypes = [c_void_p, c_void_p]
    self.lib.lammps_map_atom.restype = c_int

    self.lib.lammps_get_thermo.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_get_thermo.restype = c_double

    self.lib.lammps_last_thermo.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_last_thermo.restype = c_void_p

    self.lib.lammps_encode_image_flags.restype = self.c_imageint

    self.lib.lammps_config_has_package.argtypes = [c_char_p]
    self.lib.lammps_config_package_name.argtypes = [c_int, c_char_p, c_int]
    self.lib.lammps_config_accelerator.argtypes = [c_char_p, c_char_p, c_char_p]

    self.lib.lammps_set_variable.argtypes = [c_void_p, c_char_p, c_char_p]
    self.lib.lammps_set_string_variable.argtypes = [c_void_p, c_char_p, c_char_p]
    self.lib.lammps_set_internal_variable.argtypes = [c_void_p, c_char_p, c_double]

    self.lib.lammps_has_style.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_style_count.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_style_name.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int]

    self.lib.lammps_has_id.argtypes = [c_void_p, c_char_p, c_char_p]

    self.lib.lammps_id_count.argtypes = [c_void_p, c_char_p]

    self.lib.lammps_id_name.argtypes = [c_void_p, c_char_p, c_int, c_char_p, c_int]

    self.lib.lammps_plugin_count.argtypes = [ ]
    self.lib.lammps_plugin_name.argtypes = [c_int, c_char_p, c_char_p, c_int]

    self.lib.lammps_version.argtypes = [c_void_p]

    self.lib.lammps_get_os_info.argtypes = [c_char_p, c_int]
    self.lib.lammps_get_gpu_device_info.argtypes = [c_char_p, c_int]

    self.lib.lammps_get_mpi_comm.argtypes = [c_void_p]

    self.lib.lammps_decode_image_flags.argtypes = [self.c_imageint, POINTER(c_int*3)]

    self.lib.lammps_extract_atom.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_atom_datatype.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_atom_datatype.restype = c_int

    self.lib.lammps_extract_fix.argtypes = [c_void_p, c_char_p, c_int, c_int, c_int, c_int]

    self.lib.lammps_extract_variable.argtypes = [c_void_p, c_char_p, c_char_p]
    self.lib.lammps_extract_variable_datatype.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_extract_variable_datatype.restype = c_int

    self.lib.lammps_fix_external_get_force.argtypes = [c_void_p, c_char_p]
    self.lib.lammps_fix_external_get_force.restype = POINTER(POINTER(c_double))

    self.lib.lammps_fix_external_set_energy_global.argtypes = [c_void_p, c_char_p, c_double]
    self.lib.lammps_fix_external_set_virial_global.argtypes = [c_void_p, c_char_p, POINTER(c_double)]
    self.lib.lammps_fix_external_set_energy_peratom.argtypes = [c_void_p, c_char_p, POINTER(c_double)]
    self.lib.lammps_fix_external_set_virial_peratom.argtypes = [c_void_p, c_char_p, POINTER(POINTER(c_double))]

    self.lib.lammps_fix_external_set_vector_length.argtypes = [c_void_p, c_char_p, c_int]
    self.lib.lammps_fix_external_set_vector.argtypes = [c_void_p, c_char_p, c_int, c_double]

    # detect if Python is using a version of mpi4py that can pass communicators
    # only needed if LAMMPS has been compiled with MPI support.
    self.has_mpi4py = False
    if self.has_mpi_support:
      try:
        from mpi4py import __version__ as mpi4py_version
        # tested to work with mpi4py versions 2 and 3
        self.has_mpi4py = mpi4py_version.split('.')[0] in ['2','3']
      except ImportError:
        # ignore failing import
        pass

    # if no ptr provided, create an instance of LAMMPS
    #   we can pass an MPI communicator from mpi4py v2.0.0 and later
    #   no_mpi call lets LAMMPS use MPI_COMM_WORLD
    #   cargs = array of C strings from args
    # if ptr, then are embedding Python in LAMMPS input script
    #   ptr is the desired instance of LAMMPS
    #   just convert it to ctypes ptr and store in self.lmp

    if ptr is None:

      # with mpi4py v2+, we can pass MPI communicators to LAMMPS
      # need to adjust for type of MPI communicator object
      # allow for int (like MPICH) or void* (like OpenMPI)
      if self.has_mpi_support and self.has_mpi4py:
        from mpi4py import MPI
        self.MPI = MPI

      if comm is not None:
        if not self.has_mpi_support:
          raise Exception('LAMMPS not compiled with real MPI library')
        if not self.has_mpi4py:
          raise Exception('Python mpi4py version is not 2 or 3')
        if self.MPI._sizeof(self.MPI.Comm) == sizeof(c_int):
          MPI_Comm = c_int
        else:
          MPI_Comm = c_void_p

        # Detect whether LAMMPS and mpi4py definitely use different MPI libs
        if sizeof(MPI_Comm) != self.lib.lammps_config_has_mpi_support():
          raise Exception('Inconsistent MPI library in LAMMPS and mpi4py')

        narg = 0
        cargs = None
        if cmdargs is not None:
          myargs = ["lammps".encode()]
          narg = len(cmdargs) + 1
          for arg in cmdargs:
            if type(arg) is str:
              myargs.append(arg.encode())
            elif type(arg) is bytes:
              myargs.append(arg)
            else:
              raise TypeError('Unsupported cmdargs type ', type(arg))
          cargs = (c_char_p*(narg+1))(*myargs)
          cargs[narg] = None
          self.lib.lammps_open.argtypes = [c_int, c_char_p*(narg+1), MPI_Comm, c_void_p]
        else:
          self.lib.lammps_open.argtypes = [c_int, c_char_p, MPI_Comm, c_void_p]

        self.opened = 1
        comm_ptr = self.MPI._addressof(comm)
        comm_val = MPI_Comm.from_address(comm_ptr)
        self.lmp = c_void_p(self.lib.lammps_open(narg,cargs,comm_val,None))

      else:
        if self.has_mpi4py and self.has_mpi_support:
          self.comm = self.MPI.COMM_WORLD
        self.opened = 1
        if cmdargs is not None:
          myargs = ["lammps".encode()]
          narg = len(cmdargs) + 1
          for arg in cmdargs:
            if type(arg) is str:
              myargs.append(arg.encode())
            elif type(arg) is bytes:
              myargs.append(arg)
            else:
              raise TypeError('Unsupported cmdargs type ', type(arg))
          cargs = (c_char_p*(narg+1))(*myargs)
          cargs[narg] = None
          self.lib.lammps_open_no_mpi.argtypes = [c_int, c_char_p*(narg+1), c_void_p]
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

    # check if library initilialization failed
    if not self.lmp:
      raise(RuntimeError("Failed to initialize LAMMPS object"))

    # optional numpy support (lazy loading)
    self._numpy = None

    self._installed_packages = None
    self._available_styles = None

    # check if liblammps version matches the installed python module version
    # but not for in-place usage, i.e. when the version is 0
    import lammps
    if lammps.__version__ > 0 and lammps.__version__ != self.lib.lammps_version(self.lmp):
        raise(AttributeError("LAMMPS Python module installed for LAMMPS version %d, but shared library is version %d" \
                % (lammps.__version__, self.lib.lammps_version(self.lmp))))

    # add way to insert Python callback for fix external
    self.callback = {}
    self.FIX_EXTERNAL_CALLBACK_FUNC = CFUNCTYPE(None, py_object, self.c_bigint, c_int, POINTER(self.c_tagint), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)))
    self.lib.lammps_set_fix_external_callback.argtypes = [c_void_p, c_char_p, self.FIX_EXTERNAL_CALLBACK_FUNC, py_object]
    self.lib.lammps_set_fix_external_callback.restype = None

  # -------------------------------------------------------------------------
  # shut-down LAMMPS instance

  def __del__(self):
    self.close()

  # -------------------------------------------------------------------------
  # context manager implementation

  def __enter__(self):
    return self

  def __exit__(self, ex_type, ex_value, ex_traceback):
    self.close()

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
      from .numpy_wrapper import numpy_wrapper
      self._numpy = numpy_wrapper(self)
    return self._numpy

  # -------------------------------------------------------------------------

  def close(self):
    """Explicitly delete a LAMMPS instance through the C-library interface.

    This is a wrapper around the :cpp:func:`lammps_close` function of the C-library interface.
    """
    if self.lmp and self.opened:
      self.lib.lammps_close(self.lmp)
    self.lmp = None
    self.opened = 0

  # -------------------------------------------------------------------------

  def finalize(self):
    """Shut down the MPI communication and Kokkos environment (if active) through the
       library interface by  calling :cpp:func:`lammps_mpi_finalize` and
       :cpp:func:`lammps_kokkos_finalize`.

       You cannot create or use any LAMMPS instances after this function is called
       unless LAMMPS was compiled without MPI and without Kokkos support.
    """
    self.close()
    self.lib.lammps_kokkos_finalize()
    self.lib.lammps_mpi_finalize()

  # -------------------------------------------------------------------------

  def error(self, error_type, error_text):
    """Forward error to the LAMMPS Error class.

    .. versionadded:: 3Nov2022

    This is a wrapper around the :cpp:func:`lammps_error` function of the C-library interface.

    :param error_type:
    :type error_type:  int
    :param error_text:
    :type error_text:  string
    """
    if error_text: error_text = error_text.encode()
    else: error_text = "(unknown error)".encode()

    with ExceptionCheck(self):
      self.lib.lammps_error(self.lmp, error_type, error_text)

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
    return sb.value.decode()

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

    with ExceptionCheck(self):
      self.lib.lammps_file(self.lmp, path)

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

    with ExceptionCheck(self):
      self.lib.lammps_command(self.lmp,cmd)

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

    with ExceptionCheck(self):
      self.lib.lammps_commands_list(self.lmp,narg,args)

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

    with ExceptionCheck(self):
      self.lib.lammps_commands_string(self.lmp,c_char_p(multicmd))

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

    with ExceptionCheck(self):
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
    with ExceptionCheck(self):
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

    with ExceptionCheck(self):
      return self.lib.lammps_get_thermo(self.lmp,name)

  # -------------------------------------------------------------------------
  @property
  def last_thermo_step(self):
    """ Get the last timestep where thermodynamic data was computed

    :return: the timestep or a negative number if there has not been any thermo output yet
    :rtype: int
    """
    with ExceptionCheck(self):
      ptr = self.lib.lammps_last_thermo(self.lmp, c_char_p("step".encode()), 0)
    return cast(ptr, POINTER(self.c_bigint)).contents.value

  def last_thermo(self):
    """Get a dictionary of the last thermodynamic output

    This is a wrapper around the :cpp:func:`lammps_last_thermo`
    function of the C-library interface.  It collects the cached thermo
    data from the last timestep into a dictionary.  The return value
    is None, if there has not been any thermo output yet.

    :return: a dictionary containing the last computed thermo output values
    :rtype: dict or None
    """

    rv = dict()
    mystep = self.last_thermo_step
    if mystep < 0:
      return None

    with ExceptionCheck(self):
      ptr = self.lib.lammps_last_thermo(self.lmp, c_char_p("num".encode()), 0)
    nfield = cast(ptr, POINTER(c_int)).contents.value

    for i in range(nfield):
      with ExceptionCheck(self):
        ptr = self.lib.lammps_last_thermo(self.lmp, c_char_p("keyword".encode()), i)
      kw = cast(ptr, c_char_p).value.decode()

      with ExceptionCheck(self):
        ptr = self.lib.lammps_last_thermo(self.lmp, c_char_p("type".encode()), i)
      typ = cast(ptr, POINTER(c_int)).contents.value

      with ExceptionCheck(self):
        ptr = self.lib.lammps_last_thermo(self.lmp, c_char_p("data".encode()), i)

      if typ == LAMMPS_DOUBLE:
        val = cast(ptr, POINTER(c_double)).contents.value
      elif typ == LAMMPS_INT:
        val = cast(ptr, POINTER(c_int)).contents.value
      elif typ == LAMMPS_INT64:
        val = cast(ptr, POINTER(c_int64)).contents.value
      else:
        # we should not get here
        raise TypeError("Unknown LAMMPS data type " + str(typ))
      rv[kw] = val

    return rv

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

    This is a wrapper around the :cpp:func:`lammps_extract_global` function
    of the C-library interface.  Since there are no pointers in Python, this
    method will - unlike the C function - return the value or a list of
    values.  The :cpp:func:`lammps_extract_global` documentation includes a
    list of the supported keywords and their data types.
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
    :return: value of the property or list of values or None
    :rtype: int, float, list, or NoneType
    """

    if dtype == LAMMPS_AUTODETECT:
      dtype = self.extract_global_datatype(name)

    # set length of vector for items that are not a scalar
    vec_dict = { 'boxlo':3, 'boxhi':3, 'sublo':3, 'subhi':3,
                 'sublo_lambda':3, 'subhi_lambda':3, 'periodicity':3,
                 'special_lj':4, 'special_coul':4, 'procgrid':3 }
    if name in vec_dict:
      veclen = vec_dict[name]
    elif name == 'respa_dt':
      veclen = self.extract_global('respa_levels',LAMMPS_INT)
    elif name == 'sametag':
      veclen = self.extract_setting('nall')
    else:
      veclen = 1

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
      target_type = str
    else:
      target_type = None

    ptr = self.lib.lammps_extract_global(self.lmp, name)
    if ptr:
      if dtype == LAMMPS_STRING:
        return ptr.decode('utf-8')
      if veclen > 1:
        result = []
        for i in range(0,veclen):
          result.append(target_type(ptr[i]))
        return result
      else: return target_type(ptr[0])
    return None

  # -------------------------------------------------------------------------
  # map global atom ID to local atom index

  def map_atom(self, id):
    """Map a global atom ID (aka tag) to the local atom index

    This is a wrapper around the :cpp:func:`lammps_map_atom`
    function of the C-library interface.

    :param id: atom ID
    :type id:  int
    :return: local index
    :rtype: int
    """

    tag = self.c_tagint(id)
    return self.lib.lammps_map_atom(self.lmp, pointer(tag))

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

  def extract_compute(self,cid,cstyle,ctype):
    """Retrieve data from a LAMMPS compute

    This is a wrapper around the :cpp:func:`lammps_extract_compute`
    function of the C-library interface.
    This function returns ``None`` if either the compute id is not
    recognized, or an invalid combination of :ref:`cstyle <py_style_constants>`
    and :ref:`ctype <py_type_constants>` constants is used. The
    names and functionality of the constants are the same as for
    the corresponding C-library function.  For requests to return
    a scalar or a size, the value is returned, otherwise a pointer.

    :param cid: compute ID
    :type cid:  string
    :param cstyle: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type cstyle:  int
    :param ctype: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type ctype:  int
    :return: requested data as scalar, pointer to 1d or 2d double array, or None
    :rtype: c_double, ctypes.POINTER(c_double), ctypes.POINTER(ctypes.POINTER(c_double)), or NoneType
    """
    if cid: cid = cid.encode()
    else: return None

    if ctype == LMP_TYPE_SCALAR:
      if cstyle == LMP_STYLE_GLOBAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_double)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
        return ptr[0]
      elif cstyle == LMP_STYLE_ATOM:
        return None
      elif cstyle == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
        return ptr[0]

    elif ctype == LMP_TYPE_VECTOR:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
      return ptr

    elif ctype == LMP_TYPE_ARRAY:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
      return ptr

    elif ctype == LMP_SIZE_COLS:
      if cstyle == LMP_STYLE_GLOBAL or cstyle == LMP_STYLE_ATOM or cstyle == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
        return ptr[0]

    elif ctype == LMP_SIZE_VECTOR or ctype == LMP_SIZE_ROWS:
      if cstyle == LMP_STYLE_GLOBAL or cstyle == LMP_STYLE_LOCAL:
        self.lib.lammps_extract_compute.restype = POINTER(c_int)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_compute(self.lmp,cid,cstyle,ctype)
        return ptr[0]

    return None

  # -------------------------------------------------------------------------
  # extract fix info
  # in case of global data, free memory for 1 double via lammps_free()
  # double was allocated by library interface function

  def extract_fix(self,fid,fstyle,ftype,nrow=0,ncol=0):
    """Retrieve data from a LAMMPS fix

    This is a wrapper around the :cpp:func:`lammps_extract_fix` function
    of the C-library interface.  This function returns ``None`` if
    either the fix id is not recognized, or an invalid combination of
    :ref:`fstyle <py_style_constants>` and :ref:`ftype
    <py_type_constants>` constants is used. The names and functionality
    of the constants are the same as for the corresponding C-library
    function.  For requests to return a scalar or a size, the value is
    returned, also when accessing global vectors or arrays, otherwise a
    pointer.

    .. note::

       When requesting global data, the fix data can only be accessed
       one item at a time without access to the whole vector or array.
       Thus this function will always return a scalar.  To access vector
       or array elements the "nrow" and "ncol" arguments need to be set
       accordingly (they default to 0).

    :param fid: fix ID
    :type fid:  string
    :param fstyle: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type fstyle:  int
    :param ftype: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type ftype:  int
    :param nrow: index of global vector element or row index of global array element
    :type nrow:  int
    :param ncol: column index of global array element
    :type ncol:  int
    :return: requested data or None
    :rtype: c_double, ctypes.POINTER(c_double), ctypes.POINTER(ctypes.POINTER(c_double)), or NoneType

    """
    if fid: fid = fid.encode()
    else: return None

    if fstyle == LMP_STYLE_GLOBAL:
      if ftype in (LMP_TYPE_SCALAR, LMP_TYPE_VECTOR, LMP_TYPE_ARRAY):
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_fix(self.lmp,fid,fstyle,ftype,nrow,ncol)
        result = ptr[0]
        self.lib.lammps_free(ptr)
        return result
      elif ftype in (LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS):
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
        with ExceptionCheck(self):
          ptr = self.lib.lammps_extract_fix(self.lmp,fid,fstyle,ftype,nrow,ncol)
        return ptr[0]
      else:
        return None

    elif fstyle == LMP_STYLE_ATOM:
      if ftype == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif ftype == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif ftype == LMP_SIZE_COLS:
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_fix(self.lmp,fid,fstyle,ftype,nrow,ncol)
      if ftype == LMP_SIZE_COLS:
        return ptr[0]
      else:
        return ptr

    elif fstyle == LMP_STYLE_LOCAL:
      if ftype == LMP_TYPE_VECTOR:
        self.lib.lammps_extract_fix.restype = POINTER(c_double)
      elif ftype == LMP_TYPE_ARRAY:
        self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      elif ftype in (LMP_TYPE_SCALAR, LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS):
        self.lib.lammps_extract_fix.restype = POINTER(c_int)
      else:
        return None
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_fix(self.lmp,fid,fstyle,ftype,nrow,ncol)
      if ftype in (LMP_TYPE_VECTOR, LMP_TYPE_ARRAY):
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

  def extract_variable(self, name, group=None, vartype=None):
    """ Evaluate a LAMMPS variable and return its data

    This function is a wrapper around the function
    :cpp:func:`lammps_extract_variable` of the C library interface,
    evaluates variable name and returns a copy of the computed data.
    The memory temporarily allocated by the C-interface is deleted
    after the data is copied to a Python variable or list.
    The variable must be either an equal-style (or equivalent)
    variable or an atom-style variable. The variable type can be
    provided as the ``vartype`` parameter, which may be one of several
    constants: ``LMP_VAR_EQUAL``, ``LMP_VAR_ATOM``, ``LMP_VAR_VECTOR``,
    or ``LMP_VAR_STRING``. If omitted or ``None``, LAMMPS will determine its
    value for you based on a call to
    :cpp:func:`lammps_extract_variable_datatype` from the C library interface.
    The group parameter is only used for atom-style variables and defaults to
    the group "all".

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
    if vartype is None :
      vartype = self.lib.lammps_extract_variable_datatype(self.lmp, name)
    if vartype == LMP_VAR_EQUAL:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      if ptr: result = ptr[0]
      else: return None
      self.lib.lammps_free(ptr)
      return result
    elif vartype == LMP_VAR_ATOM:
      nlocal = self.extract_global("nlocal")
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      with ExceptionCheck(self):
        ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      if ptr:
        for i in range(nlocal): result[i] = ptr[i]
        self.lib.lammps_free(ptr)
      else: return None
      return result
    elif vartype == LMP_VAR_VECTOR :
      nvector = 0
      self.lib.lammps_extract_variable.restype = POINTER(c_int)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,
              'LMP_SIZE_VECTOR'.encode())
      if ptr :
        nvector = ptr[0]
        self.lib.lammps_free(ptr)
      else :
        return None
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      result = (c_double*nvector)()
      values = self.lib.lammps_extract_variable(self.lmp,name,group)
      if values :
        for i in range(nvector) :
          result[i] = values[i]
        # do NOT free the values pointer (points to internal vector data)
        return result
      else :
        return None
    elif vartype == LMP_VAR_STRING :
      self.lib.lammps_extract_variable.restype = c_char_p
      with ExceptionCheck(self) :
        ptr = self.lib.lammps_extract_variable(self.lmp, name, group)
        return ptr.decode('utf-8')
    return None

  # -------------------------------------------------------------------------

  def flush_buffers(self):
    """Flush output buffers

    This is a wrapper around the :cpp:func:`lammps_flush_buffers`
    function of the C-library interface.
    """
    self.lib.lammps_flush_buffers(self.lmp)

  # -------------------------------------------------------------------------

  def set_variable(self,name,value):
    """Set a new value for a LAMMPS string style variable

    .. deprecated:: 7Feb2024

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
    with ExceptionCheck(self):
      return self.lib.lammps_set_variable(self.lmp,name,value)

  # -------------------------------------------------------------------------

  def set_string_variable(self,name,value):
    """Set a new value for a LAMMPS string style variable

    .. versionadded:: 7Feb2024

    This is a wrapper around the :cpp:func:`lammps_set_string_variable`
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
    with ExceptionCheck(self):
      return self.lib.lammps_set_string_variable(self.lmp,name,value)

  # -------------------------------------------------------------------------

  def set_internal_variable(self,name,value):
    """Set a new value for a LAMMPS internal style variable

    .. versionadded:: 7Feb2024

    This is a wrapper around the :cpp:func:`lammps_set_internal_variable`
    function of the C-library interface.

    :param name: name of the variable
    :type name: string
    :param value: new variable value
    :type value: float or compatible. will be converted to float
    :return: either 0 on success or -1 on failure
    :rtype: int
    """
    if name: name = name.encode()
    else: return -1
    with ExceptionCheck(self):
      return self.lib.lammps_set_internal_variable(self.lmp,name,value)

  # -------------------------------------------------------------------------

  # return vector of atom properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # dtype = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to ensure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def gather_atoms(self,name,dtype,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*natoms)*c_int)()
        self.lib.lammps_gather_atoms(self.lmp,name,dtype,count,data)
      elif dtype == 1:
        data = ((count*natoms)*c_double)()
        self.lib.lammps_gather_atoms(self.lmp,name,dtype,count,data)
      else:
        return None
    return data

  # -------------------------------------------------------------------------

  def gather_atoms_concat(self,name,dtype,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*natoms)*c_int)()
        self.lib.lammps_gather_atoms_concat(self.lmp,name,dtype,count,data)
      elif dtype == 1:
        data = ((count*natoms)*c_double)()
        self.lib.lammps_gather_atoms_concat(self.lmp,name,dtype,count,data)
      else:
          return None
    return data

  def gather_atoms_subset(self,name,dtype,count,ndata,ids):
    if name: name = name.encode()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*ndata)*c_int)()
        self.lib.lammps_gather_atoms_subset(self.lmp,name,dtype,count,ndata,ids,data)
      elif dtype == 1:
        data = ((count*ndata)*c_double)()
        self.lib.lammps_gather_atoms_subset(self.lmp,name,dtype,count,ndata,ids,data)
      else:
        return None
    return data

  # -------------------------------------------------------------------------

  # scatter vector of atom properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to ensure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter_atoms(self,name,dtype,count,data):
    if name: name = name.encode()
    with ExceptionCheck(self):
      self.lib.lammps_scatter_atoms(self.lmp,name,dtype,count,data)

  # -------------------------------------------------------------------------

  def scatter_atoms_subset(self,name,dtype,count,ndata,ids,data):
    if name: name = name.encode()
    with ExceptionCheck(self):
      self.lib.lammps_scatter_atoms_subset(self.lmp,name,dtype,count,ndata,ids,data)


  # -------------------------------------------------------------------------

  def gather_bonds(self):
    """Retrieve global list of bonds

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_gather_bonds`
    function of the C-library interface.

    This function returns a tuple with the number of bonds and a
    flat list of ctypes integer values with the bond type, bond atom1,
    bond atom2 for each bond.

    :return: a tuple with the number of bonds and a list of c_int or c_long
    :rtype: (int, 3*nbonds*c_tagint)
    """
    nbonds = self.extract_global("nbonds")
    with ExceptionCheck(self):
        data = ((3*nbonds)*self.c_tagint)()
        self.lib.lammps_gather_bonds(self.lmp,data)
        return nbonds,data

  # -------------------------------------------------------------------------

  def gather_angles(self):
    """Retrieve global list of angles

    .. versionadded:: 8Feb2023

    This is a wrapper around the :cpp:func:`lammps_gather_angles`
    function of the C-library interface.

    This function returns a tuple with the number of angles and a
    flat list of ctypes integer values with the angle type, angle atom1,
    angle atom2, angle atom3 for each angle.

    :return: a tuple with the number of angles and a list of c_int or c_long
    :rtype: (int, 4*nangles*c_tagint)
    """
    nangles = self.extract_global("nangles")
    with ExceptionCheck(self):
        data = ((4*nangles)*self.c_tagint)()
        self.lib.lammps_gather_angles(self.lmp,data)
        return nangles,data

  # -------------------------------------------------------------------------

  def gather_dihedrals(self):
    """Retrieve global list of dihedrals

    .. versionadded:: 8Feb2023

    This is a wrapper around the :cpp:func:`lammps_gather_dihedrals`
    function of the C-library interface.

    This function returns a tuple with the number of dihedrals and a
    flat list of ctypes integer values with the dihedral type, dihedral atom1,
    dihedral atom2, dihedral atom3, dihedral atom4 for each dihedral.

    :return: a tuple with the number of dihedrals and a list of c_int or c_long
    :rtype: (int, 5*ndihedrals*c_tagint)
    """
    ndihedrals = self.extract_global("ndihedrals")
    with ExceptionCheck(self):
        data = ((5*ndihedrals)*self.c_tagint)()
        self.lib.lammps_gather_dihedrals(self.lmp,data)
        return ndihedrals,data

  # -------------------------------------------------------------------------

  def gather_impropers(self):
    """Retrieve global list of impropers

    .. versionadded:: 8Feb2023

    This is a wrapper around the :cpp:func:`lammps_gather_impropers`
    function of the C-library interface.

    This function returns a tuple with the number of impropers and a
    flat list of ctypes integer values with the improper type, improper atom1,
    improper atom2, improper atom3, improper atom4 for each improper.

    :return: a tuple with the number of impropers and a list of c_int or c_long
    :rtype: (int, 5*nimpropers*c_tagint)
    """
    nimpropers = self.extract_global("nimpropers")
    with ExceptionCheck(self):
        data = ((5*nimpropers)*self.c_tagint)()
        self.lib.lammps_gather_impropers(self.lmp,data)
        return nimpropers,data

  # -------------------------------------------------------------------------

  # return vector of atom/compute/fix properties gathered across procs
  # 3 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # returned data is a 1d vector - doc how it is ordered?
  # NOTE: need to ensure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes
  def gather(self,name,dtype,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*natoms)*c_int)()
        self.lib.lammps_gather(self.lmp,name,dtype,count,data)
      elif dtype == 1:
        data = ((count*natoms)*c_double)()
        self.lib.lammps_gather(self.lmp,name,dtype,count,data)
      else:
        return None
    return data

  def gather_concat(self,name,dtype,count):
    if name: name = name.encode()
    natoms = self.get_natoms()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*natoms)*c_int)()
        self.lib.lammps_gather_concat(self.lmp,name,dtype,count,data)
      elif dtype == 1:
        data = ((count*natoms)*c_double)()
        self.lib.lammps_gather_concat(self.lmp,name,dtype,count,data)
      else:
        return None
    return data

  def gather_subset(self,name,dtype,count,ndata,ids):
    if name: name = name.encode()
    with ExceptionCheck(self):
      if dtype == 0:
        data = ((count*ndata)*c_int)()
        self.lib.lammps_gather_subset(self.lmp,name,dtype,count,ndata,ids,data)
      elif dtype == 1:
        data = ((count*ndata)*c_double)()
        self.lib.lammps_gather_subset(self.lmp,name,dtype,count,ndata,ids,data)
      else:
        return None
    return data

  # scatter vector of atom/compute/fix properties across procs
  # 2 variants to match src/library.cpp
  # name = atom property recognized by LAMMPS in atom->extract()
  # type = 0 for integer values, 1 for double values
  # count = number of per-atom valus, 1 for type or charge, 3 for x or f
  # assume data is of correct type and length, as created by gather_atoms()
  # NOTE: need to ensure are converting to/from correct Python type
  #   e.g. for Python list or NumPy or ctypes

  def scatter(self,name,dtype,count,data):
    if name: name = name.encode()
    with ExceptionCheck(self):
      self.lib.lammps_scatter(self.lmp,name,dtype,count,data)

  def scatter_subset(self,name,dtype,count,ndata,ids,data):
    if name: name = name.encode()
    with ExceptionCheck(self):
      self.lib.lammps_scatter_subset(self.lmp,name,dtype,count,ndata,ids,data)

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
  # NOTE: how could we ensure are passing correct type to LAMMPS
  #   e.g. for Python list or NumPy, etc
  #   ditto for gather_atoms() above

  def create_atoms(self,n,id,type,x,v=None,image=None,shrinkexceed=False):
    """
    Create N atoms from list of coordinates and properties

    This function is a wrapper around the :cpp:func:`lammps_create_atoms`
    function of the C-library interface, and the behavior is similar except
    that the *v*, *image*, and *shrinkexceed* arguments are optional and
    default to *None*, *None*, and *False*, respectively. With *None* being
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
    if id is not None:
      id_lmp = (self.c_tagint*n)()
      try:
        id_lmp[:] = id[0:n]
      except ValueError:
        return 0
    else:
      id_lmp = None

    type_lmp = (c_int*n)()
    try:
      type_lmp[:] = type[0:n]
    except ValueError:
      return 0

    three_n = 3*n
    x_lmp = (c_double*three_n)()
    try:
      x_lmp[:] = x[0:three_n]
    except ValueError:
      return 0

    if v is not None:
      v_lmp = (c_double*(three_n))()
      try:
        v_lmp[:] = v[0:three_n]
      except ValueError:
        return 0
    else:
      v_lmp = None

    if image is not None:
      img_lmp = (self.c_imageint*n)()
      try:
        img_lmp[:] = image[0:n]
      except ValueError:
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
    with ExceptionCheck(self):
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

    .. versionadded:: 9Oct2020

    Various LAMMPS commands must not be called during an ongoing
    run or minimization.  This property allows to check for that.
    This is a wrapper around the :cpp:func:`lammps_is_running`
    function of the library interface.

    :return: True when called during a run otherwise false
    :rtype: bool
    """
    return self.lib.lammps_is_running(self.lmp) == 1

  # -------------------------------------------------------------------------

  def force_timeout(self):
    """ Trigger an immediate timeout, i.e. a "soft stop" of a run.

    .. versionadded:: 9Oct2020

    This function allows to cleanly stop an ongoing run or minimization
    at the next loop iteration.
    This is a wrapper around the :cpp:func:`lammps_force_timeout`
    function of the library interface.

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

  def has_package(self, name):
    """ Report if the named package has been enabled in the LAMMPS shared library.

    .. versionadded:: 3Nov2022

    This is a wrapper around the :cpp:func:`lammps_config_has_package`
    function of the library interface.

    :param name: name of the package
    :type  name: string

    :return: state of package availability
    :rtype: bool
    """
    return self.lib.lammps_config_has_package(name.encode()) != 0
  # -------------------------------------------------------------------------

  @property
  def accelerator_config(self):
    """ Return table with available accelerator configuration settings.

    This is a wrapper around the :cpp:func:`lammps_config_accelerator`
    function of the library interface which loops over all known packages
    and categories and returns enabled features as a nested dictionary
    with all enabled settings as list of strings.

    :return: nested dictionary with all known enabled settings as list of strings
    :rtype: dictionary
    """

    result = {}
    for p in ['GPU', 'KOKKOS', 'INTEL', 'OPENMP']:
      result[p] = {}
      c = 'api'
      result[p][c] = []
      for s in ['cuda', 'hip', 'phi', 'pthreads', 'opencl', 'openmp', 'serial']:
        if self.lib.lammps_config_accelerator(p.encode(),c.encode(),s.encode()):
          result[p][c].append(s)
      c = 'precision'
      result[p][c] = []
      for s in ['double', 'mixed', 'single']:
        if self.lib.lammps_config_accelerator(p.encode(),c.encode(),s.encode()):
          result[p][c].append(s)
    return result

  # -------------------------------------------------------------------------

  @property
  def has_gpu_device(self):
    """ Availability of GPU package compatible device

    This is a wrapper around the :cpp:func:`lammps_has_gpu_device`
    function of the C library interface.

    :return: True if a GPU package compatible device is present, otherwise False
    :rtype: bool
    """
    return self.lib.lammps_has_gpu_device() != 0

  # -------------------------------------------------------------------------

  def get_gpu_device_info(self):
    """Return a string with detailed information about any devices that are
    usable by the GPU package.

    This is a wrapper around the :cpp:func:`lammps_get_gpu_device_info`
    function of the C-library interface.

    :return: GPU device info string
    :rtype:  string
    """

    sb = create_string_buffer(8192)
    self.lib.lammps_get_gpu_device_info(sb,8192)
    return sb.value.decode()

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
      with ExceptionCheck(self):
        nstyles = self.lib.lammps_style_count(self.lmp, category.encode())
      sb = create_string_buffer(100)
      for idx in range(nstyles):
        with ExceptionCheck(self):
          self.lib.lammps_style_name(self.lmp, category.encode(), idx, sb, 100)
        self._available_styles[category].append(sb.value.decode())
    return self._available_styles[category]

  # -------------------------------------------------------------------------

  def has_id(self, category, name):
    """Returns whether a given ID name is available in a given category

    .. versionadded:: 9Oct2020

    This is a wrapper around the function :cpp:func:`lammps_has_id`
    of the library interface.

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

    .. versionadded:: 9Oct2020

    This is a wrapper around the functions :cpp:func:`lammps_id_count()`
    and :cpp:func:`lammps_id_name()` of the library interface.

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

  def available_plugins(self, category):
    """Returns a list of plugins available for a given category

    .. versionadded:: 10Mar2021

    This is a wrapper around the functions :cpp:func:`lammps_plugin_count()`
    and :cpp:func:`lammps_plugin_name()` of the library interface.

    :return: list of style/name pairs of loaded plugins
    :rtype:  list
    """

    available_plugins = []
    num = self.lib.lammps_plugin_count(self.lmp)
    sty = create_string_buffer(100)
    nam = create_string_buffer(100)
    for idx in range(num):
      self.lib.lammps_plugin_name(idx, sty, nam, 100)
      available_plugins.append([sty.value.decode(), nam.value.decode()])
    return available_plugins

  # -------------------------------------------------------------------------

  def set_fix_external_callback(self, fix_id, callback, caller=None):
    """Set the callback function for a fix external instance with a given fix ID.

    Optionally also set a reference to the calling object.

    This is a wrapper around the :cpp:func:`lammps_set_fix_external_callback` function
    of the C-library interface.  However this is set up to call a Python function with
    the following arguments.

    .. code-block: python

       def func(object, ntimestep, nlocal, tag, x, f):

    - object is the value of the "caller" argument
    - ntimestep is the current timestep
    - nlocal is the number of local atoms on the current MPI process
    - tag is a 1d NumPy array of integers representing the atom IDs of the local atoms
    - x is a 2d NumPy array of doubles of the coordinates of the local atoms
    - f is a 2d NumPy array of doubles of the forces on the local atoms that will be added

    .. versionchanged:: 28Jul2021

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param callback: Python function that will be called from fix external
    :type: function
    :param caller: reference to some object passed to the callback function
    :type: object, optional
    """
    import numpy as np

    def callback_wrapper(caller, ntimestep, nlocal, tag_ptr, x_ptr, fext_ptr):
      tag = self.numpy.iarray(self.c_tagint, tag_ptr, nlocal, 1)
      x   = self.numpy.darray(x_ptr, nlocal, 3)
      f   = self.numpy.darray(fext_ptr, nlocal, 3)
      callback(caller, ntimestep, nlocal, tag, x, f)

    cFunc   = self.FIX_EXTERNAL_CALLBACK_FUNC(callback_wrapper)
    cCaller = caller

    self.callback[fix_id] = { 'function': cFunc, 'caller': caller }
    with ExceptionCheck(self):
      self.lib.lammps_set_fix_external_callback(self.lmp, fix_id.encode(), cFunc, cCaller)

  # -------------------------------------------------------------------------

  def fix_external_get_force(self, fix_id):
    """Get access to the array with per-atom forces of a fix external instance with a given fix ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_get_force` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :return: requested data
    :rtype: ctypes.POINTER(ctypes.POINTER(ctypes.double))
    """

    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_get_force(self.lmp, fix_id.encode())

  # -------------------------------------------------------------------------

  def fix_external_set_energy_global(self, fix_id, eng):
    """Set the global energy contribution for a fix external instance with the given ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_energy_global` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param eng:     potential energy value to be added by fix external
    :type: float
    """

    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_energy_global(self.lmp, fix_id.encode(), eng)

  # -------------------------------------------------------------------------

  def fix_external_set_virial_global(self, fix_id, virial):
    """Set the global virial contribution for a fix external instance with the given ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_virial_global` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param eng:     list of 6 floating point numbers with the virial to be added by fix external
    :type: float
    """

    cvirial = (6*c_double)(*virial)
    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_virial_global(self.lmp, fix_id.encode(), cvirial)

  # -------------------------------------------------------------------------

  def fix_external_set_energy_peratom(self, fix_id, eatom):
    """Set the per-atom energy contribution for a fix external instance with the given ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_energy_peratom` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param eatom:   list of potential energy values for local atoms to be added by fix external
    :type: float
    """

    nlocal = self.extract_setting('nlocal')
    if len(eatom) < nlocal:
      raise Exception('per-atom energy list length must be at least nlocal')
    ceatom = (nlocal*c_double)(*eatom)
    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_energy_peratom(self.lmp, fix_id.encode(), ceatom)

  # -------------------------------------------------------------------------

  def fix_external_set_virial_peratom(self, fix_id, vatom):
    """Set the per-atom virial contribution for a fix external instance with the given ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_virial_peratom` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param vatom:     list of natoms lists with 6 floating point numbers to be added by fix external
    :type: float
    """

    # copy virial data to C compatible buffer
    nlocal = self.extract_setting('nlocal')
    if len(vatom) < nlocal:
      raise Exception('per-atom virial first dimension must be at least nlocal')
    if len(vatom[0]) != 6:
      raise Exception('per-atom virial second dimension must be 6')
    vbuf = (c_double * 6)
    vptr = POINTER(c_double)
    c_virial = (vptr * nlocal)()
    for i in range(nlocal):
        c_virial[i] = vbuf()
        for j in range(6):
          c_virial[i][j] = vatom[i][j]

    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_virial_peratom(self.lmp, fix_id.encode(), c_virial)

  # -------------------------------------------------------------------------
  def fix_external_set_vector_length(self, fix_id, length):
    """Set the vector length for a global vector stored with fix external for analysis

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_vector_length` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param length:  length of the global vector
    :type: int
    """

    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_vector_length(self.lmp, fix_id.encode(), length)

  # -------------------------------------------------------------------------
  def fix_external_set_vector(self, fix_id, idx, val):
    """Store a global vector value for a fix external instance with the given ID.

    .. versionadded:: 28Jul2021

    This is a wrapper around the :cpp:func:`lammps_fix_external_set_vector` function
    of the C-library interface.

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param idx:     1-based index of the value in the global vector
    :type: int
    :param val:     value to be stored in the global vector
    :type: float
    """

    with ExceptionCheck(self):
      return self.lib.lammps_fix_external_set_vector(self.lmp, fix_id.encode(), idx, val)

  # -------------------------------------------------------------------------

  def get_neighlist(self, idx):
    """Returns an instance of :class:`NeighList` which wraps access to the neighbor list with the given index

    See :py:meth:`lammps.numpy.get_neighlist() <lammps.numpy_wrapper.numpy_wrapper.get_neighlist()>` if you want to use
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

  def find_pair_neighlist(self, style, exact=True, nsub=0, reqid=0):
    """Find neighbor list index of pair style neighbor list

    Search for a neighbor list requested by a pair style instance that
    matches "style".  If exact is True, the pair style name must match
    exactly. If exact is False, the pair style name is matched against
    "style" as regular expression or sub-string. If the pair style is a
    hybrid pair style, the style is instead matched against the hybrid
    sub-styles. If the same pair style is used as sub-style multiple
    types, you must set nsub to a value n > 0 which indicates the nth
    instance of that sub-style to be used (same as for the pair_coeff
    command). The default value of 0 will fail to match in that case.

    Once the pair style instance has been identified, it may have
    requested multiple neighbor lists. Those are uniquely identified by
    a request ID > 0 as set by the pair style. Otherwise the request
    ID is 0.

    :param style: name of pair style that should be searched for
    :type  style: string
    :param exact: controls whether style should match exactly or only must be contained in pair style name, defaults to True
    :type  exact: bool, optional
    :param nsub:  match nsub-th hybrid sub-style, defaults to 0
    :type  nsub:  int, optional
    :param reqid: list request id, > 0 in case there are more than one, defaults to 0
    :type  reqid:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int

    """
    style = style.encode()
    exact = int(exact)
    idx = self.lib.lammps_find_pair_neighlist(self.lmp, style, exact, nsub, reqid)
    return idx

  # -------------------------------------------------------------------------

  def find_fix_neighlist(self, fixid, reqid=0):
    """Find neighbor list index of fix neighbor list

    The fix instance requesting the neighbor list is uniquely identified
    by the fix ID.  In case the fix has requested multiple neighbor
    lists, those are uniquely identified by a request ID > 0 as set by
    the fix.  Otherwise the request ID is 0 (the default).

    :param fixid: name of fix
    :type  fixid: string
    :param reqid:   id of neighbor list request, in case there are more than one request, defaults to 0
    :type  reqid:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int

    """
    fixid = fixid.encode()
    idx = self.lib.lammps_find_fix_neighlist(self.lmp, fixid, reqid)
    return idx

  # -------------------------------------------------------------------------

  def find_compute_neighlist(self, computeid, reqid=0):
    """Find neighbor list index of compute neighbor list

    The compute instance requesting the neighbor list is uniquely
    identified by the compute ID.  In case the compute has requested
    multiple neighbor lists, those are uniquely identified by a request
    ID > 0 as set by the compute.  Otherwise the request ID is 0 (the
    default).

    :param computeid: name of compute
    :type  computeid: string
    :param reqid:   index of neighbor list request, in case there are more than one request, defaults to 0
    :type  reqid:   int, optional
    :return: neighbor list index if found, otherwise -1
    :rtype:  int

    """
    computeid = computeid.encode()
    idx = self.lib.lammps_find_compute_neighlist(self.lmp, computeid, reqid)
    return idx
