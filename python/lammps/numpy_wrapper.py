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

################################################################################
# NumPy additions
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

import warnings
from ctypes import POINTER, c_void_p, c_char_p, c_double, c_int, c_int32, c_int64, cast


from .constants import *                # lgtm [py/polluting-import]
from .data import NeighList


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
        if name in ("x", "v", "f", "x0","omega", "angmom", "torque", "csforce", "vforce", "vest"):
          dim = 3
        elif name == "smd_data_9":
          dim = 9
        elif name == "smd_stress":
          dim = 6
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

  def extract_compute(self, cid, cstyle, ctype):
    """Retrieve data from a LAMMPS compute

    This is a wrapper around the
    :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>` method.
    It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

    :param cid: compute ID
    :type cid:  string
    :param cstyle: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
    :type cstyle:  int
    :param ctype: type of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
    :type ctype:  int
    :return: requested data either as float, as NumPy array with direct access to C data, or None
    :rtype: float, numpy.array, or NoneType
    """
    value = self.lmp.extract_compute(cid, cstyle, ctype)

    if cstyle == LMP_STYLE_GLOBAL:
      if ctype == LMP_TYPE_VECTOR:
        nrows = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_VECTOR)
        return self.darray(value, nrows)
      elif ctype == LMP_TYPE_ARRAY:
        nrows = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_ROWS)
        ncols = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_COLS)
        return self.darray(value, nrows, ncols)
    elif cstyle == LMP_STYLE_LOCAL:
      nrows = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_ROWS)
      ncols = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_COLS)
      if ncols == 0:
        return self.darray(value, nrows)
      else:
        return self.darray(value, nrows, ncols)
    elif cstyle == LMP_STYLE_ATOM:
      if ctype == LMP_TYPE_VECTOR:
        nlocal = self.lmp.extract_global("nlocal")
        return self.darray(value, nlocal)
      elif ctype == LMP_TYPE_ARRAY:
        nlocal = self.lmp.extract_global("nlocal")
        ncols = self.lmp.extract_compute(cid, cstyle, LMP_SIZE_COLS)
        return self.darray(value, nlocal, ncols)
    return value

  # -------------------------------------------------------------------------

  def extract_fix(self, fid, fstyle, ftype, nrow=0, ncol=0):
    """Retrieve data from a LAMMPS fix

    This is a wrapper around the :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>` method.
    It behaves the same as the original method, but returns NumPy arrays
    instead of ``ctypes`` pointers.

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
    :return: requested data
    :rtype: integer or double value, pointer to 1d or 2d double array  or None

    """
    value = self.lmp.extract_fix(fid, fstyle, ftype, nrow, ncol)
    if fstyle == LMP_STYLE_ATOM:
      if ftype == LMP_TYPE_VECTOR:
        nlocal = self.lmp.extract_global("nlocal")
        return self.darray(value, nlocal)
      elif ftype == LMP_TYPE_ARRAY:
        nlocal = self.lmp.extract_global("nlocal")
        ncols = self.lmp.extract_fix(fid, fstyle, LMP_SIZE_COLS, 0, 0)
        return self.darray(value, nlocal, ncols)
    elif fstyle == LMP_STYLE_LOCAL:
      if ftype == LMP_TYPE_VECTOR:
        nrows = self.lmp.extract_fix(fid, fstyle, LMP_SIZE_ROWS, 0, 0)
        return self.darray(value, nrows)
      elif ftype == LMP_TYPE_ARRAY:
        nrows = self.lmp.extract_fix(fid, fstyle, LMP_SIZE_ROWS, 0, 0)
        ncols = self.lmp.extract_fix(fid, fstyle, LMP_SIZE_COLS, 0, 0)
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

  def gather_bonds(self):
    """Retrieve global list of bonds as NumPy array

    This is a wrapper around :py:meth:`lammps.gather_bonds() <lammps.lammps.gather_bonds()>`
    It behaves the same as the original method, but returns a NumPy array instead
    of a ``ctypes`` list.

    .. versionadded:: 28Jul2021

    :return: the requested data as a 2d-integer numpy array
    :rtype: numpy.array(nbonds,3)
    """
    import numpy as np
    nbonds, value = self.lmp.gather_bonds()
    return np.ctypeslib.as_array(value).reshape(nbonds,3)

    # -------------------------------------------------------------------------

  def gather_angles(self):
    """ Retrieve global list of angles as NumPy array

    This is a wrapper around :py:meth:`lammps.gather_angles() <lammps.lammps.gather_angles()>`
    It behaves the same as the original method, but returns a NumPy array instead
    of a ``ctypes`` list.

    .. versionadded:: TBD

    :return: the requested data as a 2d-integer numpy array
    :rtype: numpy.array(nangles,4)
    """
    import numpy as np
    nangles, value = self.lmp.gather_angles()
    return np.ctypeslib.as_array(value).reshape(nangles,4)

    # -------------------------------------------------------------------------

  def gather_dihedrals(self):
    """ Retrieve global list of dihedrals as NumPy array

    This is a wrapper around :py:meth:`lammps.gather_dihedrals() <lammps.lammps.gather_dihedrals()>`
    It behaves the same as the original method, but returns a NumPy array instead
    of a ``ctypes`` list.

    .. versionadded:: TBD

    :return: the requested data as a 2d-integer numpy array
    :rtype: numpy.array(ndihedrals,5)
    """
    import numpy as np
    ndihedrals, value = self.lmp.gather_dihedrals()
    return np.ctypeslib.as_array(value).reshape(ndihedrals,5)

    # -------------------------------------------------------------------------

  def gather_impropers(self):
    """ Retrieve global list of impropers as NumPy array

    This is a wrapper around :py:meth:`lammps.gather_impropers() <lammps.lammps.gather_impropers()>`
    It behaves the same as the original method, but returns a NumPy array instead
    of a ``ctypes`` list.

    .. versionadded:: TBD

    :return: the requested data as a 2d-integer numpy array
    :rtype: numpy.array(nimpropers,5)
    """
    import numpy as np
    nimpropers, value = self.lmp.gather_impropers()
    return np.ctypeslib.as_array(value).reshape(nimpropers,5)

    # -------------------------------------------------------------------------

  def fix_external_get_force(self, fix_id):
    """Get access to the array with per-atom forces of a fix external instance with a given fix ID.

    This function is a wrapper around the
    :py:meth:`lammps.fix_external_get_force() <lammps.lammps.fix_external_get_force()>`
    method.  It behaves the same as the original method, but returns a NumPy array instead
    of a ``ctypes`` pointer.

    .. versionchanged:: 28Jul2021

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :return: requested data
    :rtype: numpy.array
    """
    import numpy as np
    nlocal = self.lmp.extract_setting('nlocal')
    value = self.lmp.fix_external_get_force(fix_id)
    return self.darray(value,nlocal,3)

    # -------------------------------------------------------------------------

  def fix_external_set_energy_peratom(self, fix_id, eatom):
    """Set the per-atom energy contribution for a fix external instance with the given ID.

    This function is an alternative to
    :py:meth:`lammps.fix_external_set_energy_peratom() <lammps.lammps.fix_external_set_energy_peratom()>`
    method.  It behaves the same as the original method, but accepts a NumPy array
    instead of a list as argument.

    .. versionadded:: 28Jul2021

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param eatom:   per-atom potential energy
    :type: numpy.array
    """
    import numpy as np
    nlocal = self.lmp.extract_setting('nlocal')
    if len(eatom) < nlocal:
      raise Exception('per-atom energy dimension must be at least nlocal')

    c_double_p = POINTER(c_double)
    value = eatom.astype(np.double)
    return self.lmp.lib.lammps_fix_external_set_energy_peratom(self.lmp.lmp, fix_id.encode(),
                                                               value.ctypes.data_as(c_double_p))

    # -------------------------------------------------------------------------

  def fix_external_set_virial_peratom(self, fix_id, vatom):
    """Set the per-atom virial contribution for a fix external instance with the given ID.

    This function is an alternative to
    :py:meth:`lammps.fix_external_set_virial_peratom() <lammps.lammps.fix_external_set_virial_peratom()>`
    method.  It behaves the same as the original method, but accepts a NumPy array
    instead of a list as argument.

    .. versionadded:: 28Jul2021

    :param fix_id:  Fix-ID of a fix external instance
    :type: string
    :param eatom:   per-atom potential energy
    :type: numpy.array
    """
    import numpy as np
    nlocal = self.lmp.extract_setting('nlocal')
    if len(vatom) < nlocal:
      raise Exception('per-atom virial first dimension must be at least nlocal')
    if len(vatom[0]) != 6:
      raise Exception('per-atom virial second dimension must be 6')

    c_double_pp = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C')

    # recast numpy array to be compatible with library interface
    value = (vatom.__array_interface__['data'][0]
                   + np.arange(vatom.shape[0])*vatom.strides[0]).astype(np.uintp)

    # change prototype to our custom type
    self.lmp.lib.lammps_fix_external_set_virial_peratom.argtypes = [ c_void_p, c_char_p, c_double_pp ]

    self.lmp.lib.lammps_fix_external_set_virial_peratom(self.lmp.lmp, fix_id.encode(), value)

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
    if raw_ptr is None:
      return None

    import numpy as np
    np_int_type = self._ctype_to_numpy_int(c_int_type)

    if dim == 1:
      ptr = cast(raw_ptr, POINTER(c_int_type * nelem))
    else:
      ptr = cast(raw_ptr[0], POINTER(c_int_type * nelem * dim))

    a = np.frombuffer(ptr.contents, dtype=np_int_type)

    if dim > 1:
      a.shape = (nelem, dim)
    else:
      a.shape = (nelem)
    return a

  # -------------------------------------------------------------------------

  def darray(self, raw_ptr, nelem, dim=1):
    if raw_ptr is None:
      return None

    import numpy as np

    if dim == 1:
      ptr = cast(raw_ptr, POINTER(c_double * nelem))
    else:
      ptr = cast(raw_ptr[0], POINTER(c_double * nelem * dim))

    a = np.frombuffer(ptr.contents)

    if dim > 1:
      a.shape = (nelem, dim)
    else:
      a.shape = (nelem)
    return a

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
        Access a specific neighbor list entry. "element" must be a number from 0 to the size-1 of the list

        :return: tuple with atom local index, numpy array of neighbor local atom indices
        :rtype:  (int, numpy.array)
        """
        iatom, neighbors = self.lmp.numpy.get_neighlist_element_neighbors(self.idx, element)
        return iatom, neighbors

    def find(self, iatom):
        """
        Find the neighbor list for a specific (local) atom iatom.
        If there is no list for iatom, None is returned.

        :return: numpy array of neighbor local atom indices
        :rtype:  numpy.array or None
        """
        inum = self.size
        for ii in range(inum):
          idx, neighbors = self.get(ii)
          if idx == iatom:
            return neighbors
        return None
