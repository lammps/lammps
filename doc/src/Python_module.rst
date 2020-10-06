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

----------

Setting up a Python virtual environment
***************************************

LAMMPS and its Python module can be installed together into a Python virtual
environment. This lets you isolate your customized Python environment from
your user or system installation. The following is a minimal working example:

.. code-block:: bash

   # create and change into build directory
   mkdir build
   cd build

   # create virtual environment
   virtualenv myenv

   # Add venv lib folder to LD_LIBRARY_PATH when activating it
   echo 'export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH' >> myenv/bin/activate

   # Add LAMMPS_POTENTIALS path when activating venv
   echo 'export LAMMPS_POTENTIALS=$VIRTUAL_ENV/share/lammps/potentials' >> myenv/bin/activate

   # activate environment
   source myenv/bin/activate

   # configure LAMMPS compilation
   # compiles as shared library with PYTHON package and C++ exceptions
   # and installs into myvenv
   (myenv)$ cmake -C ../cmake/presets/minimal.cmake    \
                  -D BUILD_SHARED_LIBS=on              \
                  -D PKG_PYTHON=on                     \
                  -D LAMMPS_EXCEPTIONS=on              \
                  -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
                  ../cmake

    # compile LAMMPS
    (myenv)$ cmake --build . --parallel

    # install LAMMPS into myvenv
    (myenv)$ cmake --install .

Creating or deleting a LAMMPS object
************************************

With the Python interface the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructors for the
:py:meth:`lammps <lammps.lammps.__init__()>`, :py:meth:`PyLammps <lammps.PyLammps.__init__()>`,
and :py:meth:`PyLammps <lammps.IPyLammps.__init__()>` classes.
Internally it will call either :cpp:func:`lammps_open` or :cpp:func:`lammps_open_no_mpi` from the C
library API to create the class instance.

All arguments are optional.  The *name* argument allows loading a
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
MPI ranks after splitting the communicator.


Here are simple examples using all three Python interfaces:

.. tabs::

   .. tab:: lammps API

      .. code-block:: python

         from lammps import lammps

         # NOTE: argv[0] is set by the lammps class constructor
         args = ["-log", "none"]
         # create LAMMPS instance
         lmp = lammps(cmdargs=args)
         # get and print numerical version code
         print("LAMMPS Version: ", lmp.version())
         # explicitly close and delete LAMMPS instance (optional)
         lmp.close()

   .. tab:: PyLammps API

      The :py:class:`PyLammps` class is a wrapper around the
      :py:class:`lammps` class and all of its lower level functions.
      By default, it will create a new instance of :py:class:`lammps` passing
      along all arguments to the constructor of :py:class:`lammps`.

      .. code-block:: python

         from lammps import PyLammps

         # NOTE: argv[0] is set by the lammps class constructor
         args = ["-log", "none"]
         # create LAMMPS instance
         L = PyLammps(cmdargs=args)
         # get and print numerical version code
         print("LAMMPS Version: ", L.version())
         # explicitly close and delete LAMMPS instance (optional)
         L.close()

      :py:class:`PyLammps` objects can also be created on top of an existing :py:class:`lammps` object:

      .. code-block:: Python

         from lammps import lammps, PyLammps
         ...
         # create LAMMPS instance
         lmp = lammps(cmdargs=args)
         # create PyLammps instance using previously created LAMMPS instance
         L = PyLammps(ptr=lmp)

      This is useful if you have to create the :py:class:`lammps <lammps.lammps>`
      instance is a specific way, but want to take advantage of the
      :py:class:`PyLammps <lammps.PyLammps>` interface.

   .. tab:: IPyLammps API

      The :py:class:`IPyLammps` class is an extension of the
      :py:class:`PyLammps` class. It has the same construction behavior. By
      default, it will create a new instance of :py:class:`lammps` passing
      along all arguments to the constructor of :py:class:`lammps`.

      .. code-block:: python

         from lammps import IPyLammps

         # NOTE: argv[0] is set by the lammps class constructor
         args = ["-log", "none"]
         # create LAMMPS instance
         L = IPyLammps(cmdargs=args)
         # get and print numerical version code
         print("LAMMPS Version: ", L.version())
         # explicitly close and delete LAMMPS instance (optional)
         L.close()

      You can also initialize IPyLammps on top of an existing :py:class:`lammps` or :py:class:`PyLammps` object:

      .. code-block:: Python

         from lammps import lammps, IPyLammps
         ...
         # create LAMMPS instance
         lmp = lammps(cmdargs=args)
         # create PyLammps instance using previously created LAMMPS instance
         L = PyLammps(ptr=lmp)

      This is useful if you have to create the :py:class:`lammps <lammps.lammps>`
      instance is a specific way, but want to take advantage of the
      :py:class:`IPyLammps <lammps.IPyLammps>` interface.

In all of the above cases, same as with the :ref:`C library API <lammps_c_api>`, this will use the
``MPI_COMM_WORLD`` communicator for the MPI library that LAMMPS was
compiled with.  The :py:func:`lmp.close() <lammps.lammps.close>` call is
optional since the LAMMPS class instance will also be deleted
automatically during the :py:class:`lammps <lammps.lammps>` class
destructor.

Executing LAMMPS commands
*************************

Once an instance of the :py:class:`lammps`, :py:class:`PyLammps`, or
:py:class:`IPyLammps` class is created, there are multiple ways to "feed" it
commands. In a way that is not very different from running a LAMMPS input
script, except that Python has many more facilities for structured
programming than the LAMMPS input script syntax. Furthermore it is possible
to "compute" what the next LAMMPS command should be.

.. tabs::

   .. tab:: lammps API

      Same as in the equivalent
      :doc:`C library functions <Library_execute>`, commands can be read from a file, a
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

   .. tab:: PyLammps/IPyLammps API

      Unlike the lammps API, the PyLammps/IPyLammps APIs allow running LAMMPS
      commands by calling equivalent member functions.

      For instance, the following LAMMPS command

      .. code-block:: LAMMPS

         region box block 0 10 0 5 -0.5 0.5

      can be executed using the following Python code if *L* is a :py:class:`lammps` instance:

      .. code-block:: Python

         L.command("region box block 0 10 0 5 -0.5 0.5")

      With the PyLammps interface, any LAMMPS command can be split up into arbitrary parts.
      These parts are then passed to a member function with the name of the command.
      For the ``region`` command that means the :code:`region` method can be called.
      The arguments of the command can be passed as one string, or
      individually.

      .. code-block:: Python

         L.region("box block", 0, 10, 0, 5, -0.5, 0.5)

      In this example all parameters except the first are Python floating-point literals. The
      PyLammps interface takes the entire parameter list and transparently
      merges it to a single command string.

      The benefit of this approach is avoiding redundant command calls and easier
      parameterization. In the original interface parameterization this needed to be done
      manually by creating formatted strings.

      .. code-block:: Python

         L.command("region box block %f %f %f %f %f %f" % (xlo, xhi, ylo, yhi, zlo, zhi))

      In contrast, methods of PyLammps accept parameters directly and will convert
      them automatically to a final command string.

      .. code-block:: Python

         L.region("box block", xlo, xhi, ylo, yhi, zlo, zhi)

      Using these facilities, the example shown for the lammps API can be rewritten as follows:

      .. code-block:: python

         from lammps import PyLammps
         L = PyLammps()
         # read commands from file 'in.melt'
         L.file('in.melt')
         # issue a single command
         L.variable('zpos', 'index', 1.0)
         # create 10 groups with 10 atoms each
         for i in range(10):
            L.group(f"g{i}", "id", f"{10*i+1}:{10*(i+1)}")

         L.clear()
         L.region("box block", 0, 2, 0, 2, 0, 2)
         L.create_box(1, "box")
         L.create_atoms(1, "single", 1.0, 1.0, "${zpos}")

----------

The ``lammps`` class API
************************

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

----------

The ``PyLammps`` class API
**************************

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
***************************

The :py:class:`IPyLammps <lammps.PyLammps>` class is an extension of
:py:class:`PyLammps <lammps.PyLammps>`, adding additional functions to
quickly display visualizations such as images and videos inside of IPython.
See the :doc:`PyLammps Tutorial <Howto_pylammps>` for examples.

.. autoclass:: lammps.IPyLammps
   :members:

----------

Additional components of the ``lammps`` module
**********************************************

The :py:mod:`lammps` module additionally contains several constants
and the :py:class:`NeighList <lammps.NeighList>` class:

.. _py_data_constants:

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
   :py:func:`lammps.extract_compute` and :py:func:`lammps.extract_fix`.

.. _py_type_constants:

Type Constants
--------------

.. py:data:: LMP_TYPE_SCALAR, LMP_TYLE_VECTOR, LMP_TYPE_ARRAY, LMP_SIZE_VECTOR, LMP_SIZE_ROWS, LMP_SIZE_COLS
   :type: int

   Constants in the :py:mod:`lammps` module to select what type of data
   to request  from computes  or fixes.  See :cpp:enum:`_LMP_TYPE_CONST`
   for the equivalent constants in the C library interface. Used in
   :py:func:`lammps.extract_compute` and :py:func:`lammps.extract_fix`.

.. _py_var_constants:

Variable Style Constants
------------------------

.. py:data:: LMP_VAR_EQUAL, LMP_VAR_ATOM
   :type: int

   Constants in the :py:mod:`lammps` module to select what style of
   variable to query when calling :py:func:`lammps.extract_variable`.

Classes representing internal objects
-------------------------------------

.. autoclass:: lammps.NeighList
   :members:
   :no-undoc-members:


LAMMPS error handling in Python
*******************************

Compiling the shared library with :ref:`C++ exception support <exceptions>` provides a better error
handling experience. Without exceptions the LAMMPS code will terminate the
current Python process with an error message.  C++ exceptions allow capturing
them on the C++ side and rethrowing them on the Python side. This way
LAMMPS errors can be handled through the Python exception handling mechanism.

.. warning::

   Capturing a LAMMPS exception in Python can still mean that the
   current LAMMPS process is in an illegal state and must be terminated. It is
   advised to save your data and terminate the Python instance as quickly as
   possible.
