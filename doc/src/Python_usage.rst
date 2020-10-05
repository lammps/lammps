.. _mpi4py_url: https://mpi4py.readthedocs.io/

.. _python_create_lammps:

Creating or deleting a LAMMPS object
************************************

With the Python interface the creation of a :cpp:class:`LAMMPS
<LAMMPS_NS::LAMMPS>` instance is included in the constructors for the
:py:class:`lammps <lammps.lammps>`, :py:class:`PyLammps <lammps.PyLammps>`,
and :py:class:`IPyLammps <lammps.IPyLammps>` classes.
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

      The :py:class:`PyLammps <lammps.PyLammps>` class is a wrapper around the
      :py:class:`lammps <lammps.lammps>` class and all of its lower level functions.
      By default, it will create a new instance of :py:class:`lammps <lammps.lammps>` passing
      along all arguments to the constructor of :py:class:`lammps <lammps.lammps>`.

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

      :py:class:`PyLammps <lammps.PyLammps>` objects can also be created on top of an existing
      :py:class:`lammps <lammps.lammps>` object:

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

      The :py:class:`IPyLammps <lammps.IPyLammps>` class is an extension of the
      :py:class:`PyLammps <lammps.PyLammps>` class. It has the same construction behavior. By
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
compiled with.

The :py:func:`lmp.close() <lammps.lammps.close()>` call is
optional since the LAMMPS class instance will also be deleted
automatically during the :py:class:`lammps <lammps.lammps>` class
destructor.

Note that you can create multiple LAMMPS objects in your Python
script, and coordinate and run multiple simulations, e.g.

.. code-block:: Python

   from lammps import lammps
   lmp1 = lammps()
   lmp2 = lammps()
   lmp1.file("in.file1")
   lmp2.file("in.file2")

Executing LAMMPS commands
*************************

Once an instance of the :py:class:`lammps <lammps.lammps>`,
:py:class:`PyLammps <lammps.PyLammps>`, or
:py:class:`IPyLammps <lammps.IPyLammps>` class is created, there are
multiple ways to "feed" it commands. In a way that is not very different from
running a LAMMPS input script, except that Python has many more facilities
for structured programming than the LAMMPS input script syntax. Furthermore
it is possible to "compute" what the next LAMMPS command should be.

.. tabs::

   .. tab:: lammps API

      Same as in the equivalent
      :doc:`C library functions <Library_execute>`, commands can be read from a file, a
      single string, a list of strings and a block of commands in a single
      multi-line string. They are processed under the same boundary conditions
      as the C library counterparts.  The example below demonstrates the use
      of :py:func:`lammps.file()`, :py:func:`lammps.command()`,
      :py:func:`lammps.commands_list()`, and :py:func:`lammps.commands_string()`:

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
      commands by calling equivalent member functions of :py:class:`PyLammps <lammps.PyLammps>`
      and :py:class:`IPyLammps <lammps.IPyLammps>` instances.

      For instance, the following LAMMPS command

      .. code-block:: LAMMPS

         region box block 0 10 0 5 -0.5 0.5

      can be executed using with the lammps AI with the following Python code if *L* is an
      instance of :py:class:`lammps <lammps.lammps>`:

      .. code-block:: Python

         L.command("region box block 0 10 0 5 -0.5 0.5")

      With the PyLammps interface, any LAMMPS command can be split up into arbitrary parts.
      These parts are then passed to a member function with the name of the command.
      For the ``region`` command that means the :code:`region()` method can be called.
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


Retrieving or setting LAMMPS system properties
**********************************************

Similar to what is described in :doc:`Library_properties`, the instances of
:py:class:`lammps <lammps.lammps>`, :py:class:`PyLammps <lammps.PyLammps>`, or
:py:class:`IPyLammps <lammps.IPyLammps>` can be used to extract different kinds
of information about the active LAMMPS instance and also to modify some of it. The
main difference between the interfaces is how the information is exposed.

While the :py:class:`lammps <lammps.lammps>` is just a thin layer that wraps C API calls,
:py:class:`PyLammps <lammps.PyLammps>` and :py:class:`IPyLammps <lammps.IPyLammps>` expose
information as objects and properties.

In some cases the data returned is a direct reference to the original data
inside LAMMPS cast to ``ctypes`` pointers. Where possible, the wrappers will
determine the ``ctypes`` data type and cast pointers accordingly. If
``numpy`` is installed arrays can also be extracted as numpy arrays, which
will access the C arrays directly and have the correct dimensions to protect
against invalid accesses.

.. warning::

   When accessing per-atom data,
   please note that this data is the per-processor local data and indexed
   accordingly. These arrays can change sizes and order at every neighbor list
   rebuild and atom sort event as atoms are migrating between sub-domains.

.. tabs::

   .. tab:: lammps API

      .. code-block:: python

         from lammps import lammps

         lmp = lammps()
         lmp.file("in.sysinit")

         natoms = lmp.get_natoms()
         print(f"running simulation with {natoms} atoms")

         lmp.command("run 1000 post no");

         for i in range(10):
            lmp.command("run 100 pre no post no")
            pe = lmp.get_thermo("pe")
            ke = lmp.get_thermo("ke")
            print(f"PE = {pe}\nKE = {ke}")

         lmp.close()

      **Methods**:

      * :py:meth:`version() <lammps.lammps.version()>`: return the numerical version id, e.g. LAMMPS 2 Sep 2015 -> 20150902
      * :py:meth:`get_thermo() <lammps.lammps.get_thermo()>`: return current value of a thermo keyword
      * :py:meth:`get_natoms() <lammps.lammps.get_natoms()>`: total # of atoms as int
      * :py:meth:`reset_box() <lammps.lammps.reset_box()>`: reset the simulation box size
      * :py:meth:`extract_setting() <lammps.lammps.extract_setting()>`: return a global setting
      * :py:meth:`extract_global() <lammps.lammps.extract_global()>`: extract a global quantity
      * :py:meth:`extract_atom() <lammps.lammps.extract_atom()>`: extract a per-atom quantity
      * :py:meth:`extract_box() <lammps.lammps.extract_box()>`: extract box info
      * :py:meth:`create_atoms() <lammps.lammps.create_atoms()>`: create N atoms with IDs, types, x, v, and image flags

      **Numpy Methods**:

      * :py:meth:`numpy.extract_atom() <lammps.numpy_wrapper.extract_atom()>`: extract a per-atom quantity as numpy array

   .. tab:: PyLammps/IPyLammps API

      In addition to the functions provided by :py:class:`lammps <lammps.lammps>`, :py:class:`PyLammps <lammps.PyLammps>` objects
      have several properties which allow you to query the system state:

      L.system
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

      **Retrieving the value of an arbitrary LAMMPS expressions**

      LAMMPS expressions can be immediately evaluated by using the ``eval`` method. The
      passed string parameter can be any expression containing global :doc:`thermo` values,
      variables, compute or fix data (see :doc:`Howto_output`):


      .. code-block:: Python

         result = L.eval("ke") # kinetic energy
         result = L.eval("pe") # potential energy

         result = L.eval("v_t/2.0")

      **Example**

      .. code-block:: python

         from lammps import PyLammps

         L = PyLammps()
         L.file("in.sysinit")

         print(f"running simulation with {L.system.natoms} atoms")

         L.run(1000, "post no");

         for i in range(10):
            L.run(100, "pre no post no")
            pe = L.eval("pe")
            ke = L.eval("ke")
            print(f"PE = {pe}\nKE = {ke}")



Retrieving or setting properties of LAMMPS objects
**************************************************

This section documents accessing or modifying data from objects like
computes, fixes, or variables in LAMMPS using the :py:mod:`lammps` module.

.. tabs::

   .. tab:: lammps API

      For :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>` and
      :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>`, the global, per-atom,
      or local data calculated by the compute or fix can be accessed. What is returned
      depends on whether the compute or fix calculates a scalar or vector or array.
      For a scalar, a single double value is returned.  If the compute or fix calculates
      a vector or array, a pointer to the internal LAMMPS data is returned, which you can
      use via normal Python subscripting.

      The one exception is that for a fix that calculates a
      global vector or array, a single double value from the vector or array
      is returned, indexed by I (vector) or I and J (array).  I,J are
      zero-based indices.
      See the :doc:`Howto output <Howto_output>` doc page for a discussion of
      global, per-atom, and local data, and of scalar, vector, and array
      data types.  See the doc pages for individual :doc:`computes <compute>`
      and :doc:`fixes <fix>` for a description of what they calculate and
      store.

      For :py:meth:`lammps.extract_variable() <lammps.lammps.extract_variable()>`,
      an :doc:`equal-style or atom-style variable <variable>` is evaluated and
      its result returned.

      For equal-style variables a single ``c_double`` value is returned and the
      group argument is ignored.  For atom-style variables, a vector of
      ``c_double`` is returned, one value per atom, which you can use via normal
      Python subscripting. The values will be zero for atoms not in the
      specified group.

      :py:meth:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.extract_compute()>`,
      :py:meth:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.extract_fix()>`, and
      :py:meth:`lammps.numpy.extract_variable() <lammps.numpy_wrapper.extract_variable()>` are
      equivlanent NumPy implementations that return NumPy arrays instead of ``ctypes`` pointers.

      The :py:meth:`lammps.set_variable() <lammps.lammps.set_variable()>` method sets an
      existing string-style variable to a new string value, so that subsequent LAMMPS
      commands can access the variable.

      **Methods**:

      * :py:meth:`lammps.extract_compute() <lammps.lammps.extract_compute()>`: extract value(s) from a compute
      * :py:meth:`lammps.extract_fix() <lammps.lammps.extract_fix()>`: extract value(s) from a fix
      * :py:meth:`lammps.extract_variable() <lammps.lammps.extract_variable()>`: extract value(s) from a variable
      * :py:meth:`lammps.set_variable() <lammps.lammps.set_variable()>`: set existing named string-style variable to value

      **NumPy Methods**:

      * :py:meth:`lammps.numpy.extract_compute() <lammps.numpy_wrapper.extract_compute()>`: extract value(s) from a compute, return arrays as numpy arrays
      * :py:meth:`lammps.numpy.extract_fix() <lammps.numpy_wrapper.extract_fix()>`: extract value(s) from a fix, return arrays as numpy arrays
      * :py:meth:`lammps.numpy.extract_variable() <lammps.numpy_wrapper.extract_variable()>`: extract value(s) from a variable, return arrays as numpy arrays


   .. tab:: PyLammps/IPyLammps API

      PyLammps and IPyLammps classes currently do not add any additional ways of
      retrieving information out of computes and fixes. This information can still be accessed by using the lammps API:

      .. code-block:: python

         L.lmp.extract_compute(...)
         L.lmp.extract_fix(...)
         # OR
         L.lmp.numpy.extract_compute(...)
         L.lmp.numpy.extract_fix(...)

      LAMMPS variables can be both defined and accessed via the :py:class:`PyLammps <lammps.PyLammps>` interface.

      To define a variable you can use the :doc:`variable <variable>` command:

      .. code-block:: Python

         L.variable("a index 2")

      A dictionary of all variables is returned by the :py:attr:`PyLammps.variables <lammps.PyLammps.variables>` property:

      you can access an individual variable by retrieving a variable object from the
      ``L.variables`` dictionary by name

      .. code-block:: Python

         a = L.variables['a']

      The variable value can then be easily read and written by accessing the value
      property of this object.

      .. code-block:: Python

         print(a.value)
         a.value = 4



Gather and Scatter Data between MPI processors
**********************************************

.. code-block:: Python

   data = lmp.gather_atoms(name,type,count)  # return per-atom property of all atoms gathered into data, ordered by atom ID
                                             # name = "x", "charge", "type", etc
   data = lmp.gather_atoms_concat(name,type,count)  # ditto, but concatenated atom values from each proc (unordered)
   data = lmp.gather_atoms_subset(name,type,count,ndata,ids)  # ditto, but for subset of Ndata atoms with IDs

   lmp.scatter_atoms(name,type,count,data)   # scatter per-atom property to all atoms from data, ordered by atom ID
                                             # name = "x", "charge", "type", etc
                                             # count = # of per-atom values, 1 or 3, etc

   lmp.scatter_atoms_subset(name,type,count,ndata,ids,data)  # ditto, but for subset of Ndata atoms with IDs


The gather methods collect peratom info of the requested type (atom
coords, atom types, forces, etc) from all processors, and returns the
same vector of values to each calling processor.  The scatter
functions do the inverse.  They distribute a vector of peratom values,
passed by all calling processors, to individual atoms, which may be
owned by different processors.

Note that the data returned by the gather methods,
e.g. gather_atoms("x"), is different from the data structure returned
by extract_atom("x") in four ways.  (1) Gather_atoms() returns a
vector which you index as x[i]; extract_atom() returns an array
which you index as x[i][j].  (2) Gather_atoms() orders the atoms
by atom ID while extract_atom() does not.  (3) Gather_atoms() returns
a list of all atoms in the simulation; extract_atoms() returns just
the atoms local to each processor.  (4) Finally, the gather_atoms()
data structure is a copy of the atom coords stored internally in
LAMMPS, whereas extract_atom() returns an array that effectively
points directly to the internal data.  This means you can change
values inside LAMMPS from Python by assigning a new values to the
extract_atom() array.  To do this with the gather_atoms() vector, you
need to change values in the vector, then invoke the scatter_atoms()
method.

For the scatter methods, the array of coordinates passed to must be a
ctypes vector of ints or doubles, allocated and initialized something
like this:

.. code-block:: Python

   from ctypes import c_double
   natoms = lmp.get_natoms()
   n3 = 3*natoms
   x = (n3*c_double)()
   x[0] = x coord of atom with ID 1
   x[1] = y coord of atom with ID 1
   x[2] = z coord of atom with ID 1
   x[3] = x coord of atom with ID 2
   ...
   x[n3-1] = z coord of atom with ID natoms
   lmp.scatter_atoms("x",1,3,x)

Alternatively, you can just change values in the vector returned by
the gather methods, since they are also ctypes vectors.

