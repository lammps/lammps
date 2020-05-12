The ``lammps`` Python module
****************************

The LAMMPS Python interface is implemented in the ``lammps.py`` file in
the ``python`` folder.  After compilation of LAMMPS, in can be installed
into a Python system folder or a user folder with ``make
install-python``.

There are multiple Python interface classes implemented:

- the :py:class:`lammps <lammps.lammps>` class. This is a wrapper around
  the C-library interface and its member functions try to replicate the
  :doc:`C-library API <pg_library>` closely.  This is the most
  feature-complete python API.
- the :py:class:`PyLammps <lammps.PyLammps>` class. This is a more high-level
  and more Python style class implemented on top of the
  :py:class:`lammps <lammps.lammps>` class.
- the :py:class:`IPyLammps <lammps.IPyLammps>` class is derived from
  :py:class:`PyLammps <lammps.PyLammps>` and adds embedded graphics
  features to conveniently include LAMMPS into `Jupyter
  <https://jupyter.org/>`_ notebooks.

----------

The ``lammps`` class API
************************

The :py:class:`lammps <lammps.lammps>` class is the core of the LAMMPS
Python interfaces.  It is a wrapper around the :doc:`LAMMPS C library
API <pg_library>` using the `Python ctypes module
<https://docs.python.org/3/library/ctypes.html>`_ and a shared library
compiled from the LAMMPS sources code.  The individual methods in this
class try to closely follow the corresponding C functions.  The handle
argument that needs to be passed to the C functions is stored internally
in the class and automatically added when calling the C library
functions. Below is a detailed documentation of the API.

.. autoclass:: lammps.lammps
   :members:

----------

The ``PyLammps`` class API
**************************

.. autoclass:: lammps.PyLammps
   :members:

----------

The ``IPyLammps`` class API
***************************

.. autoclass:: lammps.IPyLammps
   :members:

----------

Additional components of the ``lammps`` module
**********************************************

The :py:mod:`lammps` module additionally contains several constants
and the :py:class:`NeighList <lammps.NeighList>` class:

.. _py_data_constants:
.. py:data:: LAMMPS_INT, LAMMPS_DOUBLE, LAMMPS_BIGINT, LAMMPS_TAGINT, LAMMPS_STRING
   :type: int

   Constants in the :py:mod:`lammps` module to indicate how to
   cast data when the C library function returns a void pointer.
   Used in :py:func:`lammps.extract_global`.

.. _py_style_constants:
.. py:data:: LMP_STYLE_GLOBAL, LMP_STYLE_ATOM, LMP_STYLE_LOCAL
   :type: int

   Constants in the :py:mod:`lammps` module to select what style of data
   to request from computes or fixes. See :cpp:enum:`_LMP_STYLE_CONST`
   for the equivalent constants in the C library interface. Used in
   :py:func:`lammps.extract_compute` and :py:func:`lammps.extract_fix`.

.. _py_type_constants:
.. py:data:: LMP_TYPE_SCALAR, LMP_TYLE_VECTOR, LMP_TYPE_ARRAY, LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS
   :type: int

   Constants in the :py:mod:`lammps` module to select what type of data
   to request  from computes  or fixes.  See :cpp:enum:`_LMP_TYPE_CONST`
   for the equivalent constants in the C library interface. Used in
   :py:func:`lammps.extract_compute` and :py:func:`lammps.extract_fix`.

.. _py_var_constants:
.. py:data:: LMP_VAR_EQUAL, LMP_VAR_ATOM
   :type: int

   Constants in the :py:mod:`lammps` module to select what style of
   variable to query when calling :py:func:`lammps.extract_variable`.

.. autoclass:: lammps.NeighList
   :members:
   :no-undoc-members:

