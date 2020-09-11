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

Same as with the :ref:`C library API <lammps_c_api>` this will use the
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
simpler, more "Pythonic" interface to common LAMMPS functionality. LAMMPS data
structures are exposed through objects and properties. This makes Python scripts
shorter and more concise.

Creating a new instance of PyLammps
-----------------------------------

To create a PyLammps object you need to first import the class from the lammps
module. By using the default constructor, a new :py:class:`lammps <lammps.lammps>` instance is created.

.. code-block:: Python

   from lammps import PyLammps
   L = PyLammps()

You can also initialize PyLammps on top of this existing :py:class:`lammps <lammps.lammps>` object:

.. code-block:: Python

   from lammps import lammps, PyLammps
   lmp = lammps()
   L = PyLammps(ptr=lmp)

This is useful if you have create the :py:class:`lammps <lammps.lammps>`
instance is a specific way, but want to take advantage of the
:py:class:`PyLammps <lammps.PyLammps>` interface.

Commands
--------

Sending a LAMMPS command with the existing library interfaces is done using
the command method of the lammps object instance.

For instance, let's take the following LAMMPS command:

.. code-block:: LAMMPS

   region box block 0 10 0 5 -0.5 0.5

In the original interface this command can be executed with the following
Python code if *L* was a lammps instance:

.. code-block:: Python

   L.command("region box block 0 10 0 5 -0.5 0.5")

With the PyLammps interface, any command can be split up into arbitrary parts
separated by white-space, passed as individual arguments to a region method.

.. code-block:: Python

   L.region("box block", 0, 10, 0, 5, -0.5, 0.5)

Note that each parameter is set as Python literal floating-point number. In the
PyLammps interface, each command takes an arbitrary parameter list and transparently
merges it to a single command string, separating individual parameters by white-space.

The benefit of this approach is avoiding redundant command calls and easier
parameterization. In the original interface parameterization needed to be done
manually by creating formatted strings.

.. code-block:: Python

   L.command("region box block %f %f %f %f %f %f" % (xlo, xhi, ylo, yhi, zlo, zhi))

In contrast, methods of PyLammps accept parameters directly and will convert
them automatically to a final command string.

.. code-block:: Python

   L.region("box block", xlo, xhi, ylo, yhi, zlo, zhi)

System state
------------

In addition to dispatching commands directly through the PyLammps object, it
also provides several properties which allow you to query the system state.

:py:attr:`lammps.PyLammps.system`
   Is a dictionary describing the system such as the bounding box or number of atoms

L.system.xlo, L.system.xhi
   bounding box limits along x-axis

L.system.ylo, L.system.yhi
   bounding box limits along y-axis

L.system.zlo, L.system.zhi
   bounding box limits along z-axis

L.communication
   configuration of communication subsystem, such as the number of threads or processors

L.communication.nthreads
   number of threads used by each LAMMPS process

L.communication.nprocs
   number of MPI processes used by LAMMPS

L.fixes
   List of fixes in the current system

L.computes
   List of active computes in the current system

L.dump
   List of active dumps in the current system

L.groups
   List of groups present in the current system

Working with LAMMPS variables
-----------------------------

LAMMPS variables can be both defined and accessed via the PyLammps interface.

To define a variable you can use the :doc:`variable <variable>` command:

.. code-block:: Python

   L.variable("a index 2")

A dictionary of all variables is returned by L.variables

you can access an individual variable by retrieving a variable object from the
L.variables dictionary by name

.. code-block:: Python

   a = L.variables['a']

The variable value can then be easily read and written by accessing the value
property of this object.

.. code-block:: Python

   print(a.value)
   a.value = 4

Retrieving the value of an arbitrary LAMMPS expressions
-------------------------------------------------------

LAMMPS expressions can be immediately evaluated by using the eval method. The
passed string parameter can be any expression containing global thermo values,
variables, compute or fix data.

.. code-block:: Python

   result = L.eval("ke") # kinetic energy
   result = L.eval("pe") # potential energy

   result = L.eval("v_t/2.0")

Accessing atom data
-------------------

All atoms in the current simulation can be accessed by using the L.atoms list.
Each element of this list is an object which exposes its properties (id, type,
position, velocity, force, etc.).

.. code-block:: Python

   # access first atom
   L.atoms[0].id
   L.atoms[0].type

   # access second atom
   L.atoms[1].position
   L.atoms[1].velocity
   L.atoms[1].force

Some properties can also be used to set:

.. code-block:: Python

   # set position in 2D simulation
   L.atoms[0].position = (1.0, 0.0)

   # set position in 3D simulation
   L.atoms[0].position = (1.0, 0.0, 1.)

Evaluating thermo data
----------------------

Each simulation run usually produces thermo output based on system state,
computes, fixes or variables. The trajectories of these values can be queried
after a run via the L.runs list. This list contains a growing list of run data.
The first element is the output of the first run, the second element that of
the second run.

.. code-block:: Python

   L.run(1000)
   L.runs[0] # data of first 1000 time steps

   L.run(1000)
   L.runs[1] # data of second 1000 time steps

Each run contains a dictionary of all trajectories. Each trajectory is
accessible through its thermo name:

.. code-block:: Python

   L.runs[0].thermo.Step # list of time steps in first run
   L.runs[0].thermo.Ke   # list of kinetic energy values in first run

Together with matplotlib plotting data out of LAMMPS becomes simple:

.. code-block:: Python

   import matplotlib.plot as plt
   steps = L.runs[0].thermo.Step
   ke    = L.runs[0].thermo.Ke
   plt.plot(steps, ke)

Error handling with PyLammps
----------------------------

Compiling the shared library with C++ exception support provides a better error
handling experience.  Without exceptions the LAMMPS code will terminate the
current Python process with an error message.  C++ exceptions allow capturing
them on the C++ side and rethrowing them on the Python side. This way you
can handle LAMMPS errors through the Python exception handling mechanism.

.. warning::

   Capturing a LAMMPS exception in Python can still mean that the
   current LAMMPS process is in an illegal state and must be terminated. It is
   advised to save your data and terminate the Python instance as quickly as
   possible.

.. autoclass:: lammps.PyLammps
   :members:

.. autoclass:: lammps.AtomList
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

