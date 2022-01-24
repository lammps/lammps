The ``lammps`` Python module
****************************

.. py:module:: lammps

The LAMMPS Python interface is implemented as a module called :py:mod:`lammps`
which is defined in the ``lammps`` package in the ``python`` folder of the
LAMMPS source code distribution.  After compilation of LAMMPS, the module can
be installed into a Python system folder or a user folder with ``make
install-python``.  Components of the module can then loaded into a Python
session with the ``import`` command.

There are multiple Python interface classes in the :py:mod:`lammps` module:

- the :py:class:`lammps <lammps.lammps>` class. This is a wrapper around
  the C-library interface and its member functions try to replicate the
  :ref:`C-library API <lammps_c_api>` closely.  This is the most
  feature-complete Python API.
- the :py:class:`PyLammps <lammps.PyLammps>` class. This is a more high-level
  and more Python style class implemented on top of the
  :py:class:`lammps <lammps.lammps>` class.
- the :py:class:`IPyLammps <lammps.IPyLammps>` class is derived from
  :py:class:`PyLammps <lammps.PyLammps>` and adds embedded graphics
  features to conveniently include LAMMPS into `Jupyter
  <https://jupyter.org/>`_ notebooks.

.. _mpi4py_url: https://mpi4py.readthedocs.io

.. admonition:: Version check
   :class: note

   The :py:mod:`lammps` module stores the version number of the LAMMPS
   version it is installed from.  When initializing the
   :py:class:`lammps <lammps.lammps>` class, this version is checked to
   be the same as the result from :py:func:`lammps.version`, the version
   of the LAMMPS shared library that the module interfaces to.  If the
   they are not the same an AttributeError exception is raised since a
   mismatch of versions (e.g.  due to incorrect use of the
   ``LD_LIBRARY_PATH`` or ``PYTHONPATH`` environment variables can lead
   to crashes or data corruption and otherwise incorrect behavior.

.. automodule:: lammps
   :members:
   :noindex:

----------

The ``lammps`` class API
========================

The :py:class:`lammps <lammps.lammps>` class is the core of the LAMMPS
Python interfaces.  It is a wrapper around the :ref:`LAMMPS C library
API <lammps_c_api>` using the `Python ctypes module
<https://docs.python.org/3/library/ctypes.html>`_ and a shared library
compiled from the LAMMPS sources code.  The individual methods in this
class try to closely follow the corresponding C functions.  The handle
argument that needs to be passed to the C functions is stored internally
in the class and automatically added when calling the C library
functions. Below is a detailed documentation of the API.

.. autoclass:: lammps.lammps
   :members:

.. autoclass:: lammps.numpy_wrapper::numpy_wrapper
   :members:

----------

The ``PyLammps`` class API
==========================

The :py:class:`PyLammps <lammps.PyLammps>` class is a wrapper that creates a
simpler, more "Pythonic" interface to common LAMMPS functionality. LAMMPS
data structures are exposed through objects and properties. This makes Python
scripts shorter and more concise. See the :doc:`PyLammps Tutorial
<Howto_pylammps>` for an introduction on how to use this interface.

.. autoclass:: lammps.PyLammps
   :members:

.. autoclass:: lammps.AtomList
   :members:

.. autoclass:: lammps.Atom
   :members:

.. autoclass:: lammps.Atom2D
   :members:

----------

The ``IPyLammps`` class API
===========================

The :py:class:`IPyLammps <lammps.PyLammps>` class is an extension of
:py:class:`PyLammps <lammps.PyLammps>`, adding additional functions to
quickly display visualizations such as images and videos inside of IPython.
See the :doc:`PyLammps Tutorial <Howto_pylammps>` for examples.

.. autoclass:: lammps.IPyLammps
   :members:

----------

Additional components of the ``lammps`` module
==============================================

The :py:mod:`lammps` module additionally contains several constants
and the :py:class:`NeighList <lammps.NeighList>` class:

.. _py_datatype_constants:

Data Types
----------

.. py:data:: LAMMPS_INT, LAMMPS_INT_2D, LAMMPS_DOUBLE, LAMMPS_DOUBLE_2D, LAMMPS_INT64, LAMMPS_INT64_2D, LAMMPS_STRING
   :type: int

   Constants in the :py:mod:`lammps` module to indicate how to
   cast data when the C library function returns a void pointer.
   Used in :py:func:`lammps.extract_global` and :py:func:`lammps.extract_atom`.
   See :cpp:enum:`_LMP_DATATYPE_CONST` for the equivalent constants in the
   C library interface.

.. _py_style_constants:

Style Constants
---------------

.. py:data:: LMP_STYLE_GLOBAL, LMP_STYLE_ATOM, LMP_STYLE_LOCAL
   :type: int

   Constants in the :py:mod:`lammps` module to select what style of data
   to request from computes or fixes. See :cpp:enum:`_LMP_STYLE_CONST`
   for the equivalent constants in the C library interface. Used in
   :py:func:`lammps.extract_compute`, :py:func:`lammps.extract_fix`, and their NumPy variants
   :py:func:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute>` and
   :py:func:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix>`.

.. _py_type_constants:

Type Constants
--------------

.. py:data:: LMP_TYPE_SCALAR, LMP_TYPE_VECTOR, LMP_TYPE_ARRAY, LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS
   :type: int

   Constants in the :py:mod:`lammps` module to select what type of data
   to request  from computes  or fixes.  See :cpp:enum:`_LMP_TYPE_CONST`
   for the equivalent constants in the C library interface. Used in
   :py:func:`lammps.extract_compute`, :py:func:`lammps.extract_fix`, and their NumPy variants
   :py:func:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute>` and
   :py:func:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix>`.

.. _py_vartype_constants:

Variable Type Constants
------------------------

.. py:data:: LMP_VAR_EQUAL, LMP_VAR_ATOM
   :type: int

   Constants in the :py:mod:`lammps` module to select what type of
   variable to query when calling :py:func:`lammps.extract_variable`. See also: :doc:`variable command <variable>`.

Classes representing internal objects
-------------------------------------

.. autoclass:: lammps.NeighList
   :members:
   :no-undoc-members:

.. autoclass:: lammps.numpy_wrapper::NumPyNeighList
   :members:
   :no-undoc-members:
