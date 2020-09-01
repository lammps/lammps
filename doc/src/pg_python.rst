The ``lammps`` Python module
****************************

.. py:module:: lammps

The LAMMPS Python interface is implemented as a module called
:py:mod:`lammps` in the ``lammps.py`` file in the ``python`` folder of
the LAMMPS source code distribution.  After compilation of LAMMPS, the
module can be installed into a Python system folder or a user folder
with ``make install-python``.  Components of the module can then loaded
into a Python session with the ``import`` command.

There are multiple Python interface classes in the :py:mod:`lammps` module:

- the :py:class:`lammps <lammps.lammps>` class. This is a wrapper around
  the C-library interface and its member functions try to replicate the
  :doc:`C-library API <pg_library>` closely.  This is the most
  feature-complete Python API.
- the :py:class:`PyLammps <lammps.PyLammps>` class. This is a more high-level
  and more Python style class implemented on top of the
  :py:class:`lammps <lammps.lammps>` class.
- the :py:class:`IPyLammps <lammps.IPyLammps>` class is derived from
  :py:class:`PyLammps <lammps.PyLammps>` and adds embedded graphics
  features to conveniently include LAMMPS into `Jupyter
  <https://jupyter.org/>`_ notebooks.

.. _mpi4py_url: https://mpi4py.readthedocs.io

----------

Creating or deleting a LAMMPS object
************************************

With the Python interface the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructor for the
:py:func:`lammps <lammps.lammps>` class.  Internally it will call either
:cpp:func:`lammps_open` or :cpp:func:`lammps_open_no_mpi` from the C
library API to create the class instance.

All arguments are optional.  The *name* argument is to allow loading a
LAMMPS shared library that is named ``liblammps_machine.so`` instead of
the default name of ``liblammps.so``.  In most cases the latter will be
installed or used.  The *ptr* argument is for use of the
:py:mod:`lammps` module from inside a LAMMPS instance, e.g. with the
:doc:`python <python>` command, where a pointer to the already existing
:cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class instance can be passed
to the Python class and used instead of creating a new instance.  The
*comm* argument may be used in combination with the `mpi4py <mpi4py_url_>`_
module to pass an MPI communicator to LAMMPS and thus it is possible
to run the Python module like the library interface on a subset of the
MPI ranks after splitting the communicator. Here is a simple example:

.. code-block:: python

   from lammps import lammps

   # NOTE: argv[0] is set by the Python module
   args = ["-log", "none"]
   # create LAMMPS instance
   lmp = lammps(cmdargs=args)
   # get and print numerical version code
   print("LAMMPS Version: ", lmp.version())
   # explicitly close and delete LAMMPS instance (optional)
   lmp.close()

Same as with the :doc:`C library API <pg_lib_create>` this will use the
``MPI_COMM_WORLD`` communicator for the MPI library that LAMMPS was
compiled with.  The :py:func:`lmp.close() <lammps.lammps.close>` call is
optional since the LAMMPS class instance will also be deleted
automatically during the :py:class:`lammps <lammps.lammps>` class
destructor.

Executing LAMMPS commands
*************************

Once an instance of the :py:class:`lammps <lammps.lammps>` class is
created, there are multiple ways to "feed" it commands. In a way that is
not very different from running a LAMMPS input script, except that
Python has many more facilities for structured programming than the
LAMMPS input script syntax.  Furthermore it is possible to "compute"
what the next LAMMPS command should be. Same as in the equivalent `C
library functions <pg_lib_execute>`, commands can be read from a file, a
single string, a list of strings and a block of commands in a single
multi-line string. They are processed under the same boundary conditions
as the C library counterparts.  The example below demonstrates the use
of :py:func:`lammps.file`, :py:func:`lammps.command`,
:py:func:`lammps.commands_list`, and :py:func:`lammps.commands_string`:

.. code-block:: python

   from lammps import lammps

   lmp = lammps()
   # read commands from file 'in.melt'
   lmp.file('in.melt')
   # issue a single command
   lmp.command('variable zpos index 1.0')
   # create 10 groups with 10 atoms each
   cmds = ["group g{} id {}:{}".format(i,10*i+1,10*(i+1)) for i in range(10)]
   lmp.commands_list(cmds)
   # run commands from a multi-line string
   block = """
   clear
   region  box block 0 2 0 2 0 2
   create_box 1 box
   create_atoms 1 single 1.0 1.0 ${zpos}
   """
   lmp.commands_string(block)

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

