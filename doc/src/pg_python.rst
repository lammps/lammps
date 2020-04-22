LAMMPS Python APIs
******************

The LAMMPS Python interface requires a shared LAMMPS library or to be
called from inside LAMMPS via the :ref:`PYTHON <python>` package.  It
uses the `ctypes <https://docs.python.org/3/library/ctypes.html>`_ Python
module to dynamically import the compiled functions via the
:doc:`C-library API <pg_library>` into the Python interpreter.  This is
implemented in the ``lammps.py`` file in the ``python`` folder.

There are multiple Python interface classes implemented:

- the :py:class:`lammps <lammps.lammps> class`. This is a wrapper around
  the C-library interface and its member functions try to replicate the
  :doc:`C-library API <pg_library>` closely.  This is the most
  feature-complete python API.
- the :py:class:`PyLammps <lammps.PyLammps>`. This is a more high-level
  and more Python style class implemented on top of the
  :py:class:`lammps <lammps.lammps>` class.
- the :py:class:`IPyLammps <lammps.IPyLammps>` class is derived from
  :py:class:`PyLammps <lammps.PyLammps>` and adds embedded graphics
  features to conveniently embed LAMMPS into `Jupyter <https://jupyter.org/>`_
  notebooks.

The lammps module
=================

The :py:mod:`lammps` module contains various functions, constants,
and the :py:class:`lammps <lammps.lammps>`,
:py:class:`PyLammps <lammps.PyLammps>`,
and :py:class:`IPyLammps <lammps.IPyLammps>` classes that form
the LAMMPS Python API.

Here is a summary of the constants and their uses.

.. py:module:: lammps
   :synopsis: Python wrapper for the LAMMPS library API via ctypes

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

----------

The lammps class
================

The :py:class:`lammps <lammps.lammps>` class is a wrapper around
the :doc:`LAMMPS C library API <pg_library>` using the Python ``ctypes`` module
and a shared library compiled from the LAMMPS sources code.
The individual methods in this class try to closely follow the
corresponding C functions.  The handle argument that needs to be
passed to the C functions is stored internally in the class and
automatically added when calling the C library functions.

----------

.. autoclass:: lammps.lammps
   :members:

----------

The NeighList class
===================

.. autoclass:: lammps.NeighList
   :members:
   :no-undoc-members:

----------

The PyLammps class
==================

.. autoclass:: lammps.PyLammps
   :members:

----------

The IPyLammps class
===================

.. autoclass:: lammps.IPyLammps
   :members:
