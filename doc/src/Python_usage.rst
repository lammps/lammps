.. _mpi4py_url: https://mpi4py.readthedocs.io/

.. _python_create_lammps:

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

.. code-block:: Python

   lmp = lammps()           # create a LAMMPS object using the default liblammps.so library
                            # 4 optional args are allowed: name, cmdargs, ptr, comm
   lmp = lammps(ptr=lmpptr) # use lmpptr as previously created LAMMPS object
   lmp = lammps(comm=split) # create a LAMMPS object with a custom communicator, requires mpi4py 2.0.0 or later
   lmp = lammps(name="g++")   # create a LAMMPS object using the liblammps_g++.so library
   lmp = lammps(name="g++",cmdargs=list)    # add LAMMPS command-line args, e.g. list = ["-echo","screen"]

   lmp.close()              # destroy a LAMMPS object

The lines

.. code-block:: Python

   from lammps import lammps
   lmp = lammps()

create an instance of LAMMPS, wrapped in a Python class by the lammps
Python module, and return an instance of the Python class as lmp.  It
is used to make all subsequent calls to the LAMMPS library.

Additional arguments to lammps() can be used to tell Python the name
of the shared library to load or to pass arguments to the LAMMPS
instance, the same as if LAMMPS were launched from a command-line
prompt.

If the ptr argument is set like this:

.. code-block:: Python

   lmp = lammps(ptr=lmpptr)

then lmpptr must be an argument passed to Python via the LAMMPS
:doc:`python <python>` command, when it is used to define a Python
function that is invoked by the LAMMPS input script.  This mode of
calling Python from LAMMPS is described in the :doc:`Python call <Python_call>` doc page.  The variable lmpptr refers to the
instance of LAMMPS that called the embedded Python interpreter.  Using
it as an argument to lammps() allows the returned Python class
instance "lmp" to make calls to that instance of LAMMPS.  See the
:doc:`python <python>` command doc page for examples using this syntax.

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


.. code-block:: Python

   lmp.file(file)           # run an entire input script, file = "in.lj"
   lmp.command(cmd)         # invoke a single LAMMPS command, cmd = "run 100"
   lmp.commands_list(cmdlist)     # invoke commands in cmdlist = **"run 10", "run 20"**
   lmp.commands_string(multicmd)  # invoke commands in multicmd = "run 10\nrun 20"

The file(), command(), commands_list(), commands_string() methods
allow an input script, a single command, or multiple commands to be
invoked.


Retrieving or setting LAMMPS system properties
**********************************************

.. code-block:: Python

   version = lmp.version()  # return the numerical version id, e.g. LAMMPS 2 Sep 2015 -> 20150902

   size = lmp.extract_setting(name)     # return data type info

   xlo = lmp.extract_global(name,type)  # extract a global quantity
                                        # name = "boxxlo", "nlocal", etc
                                        # type = 0 = int
                                        #        1 = double

   boxlo,boxhi,xy,yz,xz,periodicity,box_change = lmp.extract_box()  # extract box info


   value = lmp.get_thermo(name)              # return current value of a thermo keyword
   natoms = lmp.get_natoms()                 # total # of atoms as int

   lmp.reset_box(boxlo,boxhi,xy,yz,xz)       # reset the simulation box size

   lmp.create_atoms(n,ids,types,x,v,image,shrinkexceed)   # create N atoms with IDs, types, x, v, and image flags


The :py:meth:`get_thermo() <lammps.lammps.get_thermo()>` method returns the current value of a thermo
keyword as a float.

The get_natoms() method returns the total number of atoms in the
simulation, as an int.


The reset_box() method resets the size and shape of the simulation
box, e.g. as part of restoring a previously extracted and saved state
of a simulation.


The extract_setting(), extract_global(), extract_box(),
extract_atom(), extract_compute(), extract_fix(), and
extract_variable() methods return values or pointers to data
structures internal to LAMMPS.

For extract_global() see the src/library.cpp file for the list of
valid names.  New names could easily be added.  A double or integer is
returned.  You need to specify the appropriate data type via the type
argument.

Retrieving or setting properties of LAMMPS objects
**************************************************

.. code-block:: Python

   coords = lmp.extract_atom(name,type)      # extract a per-atom quantity
                                             # name = "x", "type", etc
                                             # type = 0 = vector of ints
                                             #        1 = array of ints
                                             #        2 = vector of doubles
                                             #        3 = array of doubles

   eng = lmp.extract_compute(id,style,type)  # extract value(s) from a compute
   v3 = lmp.extract_fix(id,style,type,i,j)   # extract value(s) from a fix
                                             # id = ID of compute or fix
                                             # style = 0 = global data
                                             #         1 = per-atom data
                                             #         2 = local data
                                             # type = 0 = scalar
                                             #        1 = vector
                                             #        2 = array
                                             # i,j = indices of value in global vector or array

   var = lmp.extract_variable(name,group,flag)  # extract value(s) from a variable
                                                # name = name of variable
                                                # group = group ID (ignored for equal-style variables)
                                                # flag = 0 = equal-style variable
                                                #        1 = atom-style variable

   flag = lmp.set_variable(name,value)       # set existing named string-style variable to value, flag = 0 if successful

For extract_atom(), a pointer to internal LAMMPS atom-based data is
returned, which you can use via normal Python subscripting.  See the
extract() method in the src/atom.cpp file for a list of valid names.
Again, new names could easily be added if the property you want is not
listed.  A pointer to a vector of doubles or integers, or a pointer to
an array of doubles (double \*\*) or integers (int \*\*) is returned.  You
need to specify the appropriate data type via the type argument.

For extract_compute() and extract_fix(), the global, per-atom, or
local data calculated by the compute or fix can be accessed.  What is
returned depends on whether the compute or fix calculates a scalar or
vector or array.  For a scalar, a single double value is returned.  If
the compute or fix calculates a vector or array, a pointer to the
internal LAMMPS data is returned, which you can use via normal Python
subscripting.  The one exception is that for a fix that calculates a
global vector or array, a single double value from the vector or array
is returned, indexed by I (vector) or I and J (array).  I,J are
zero-based indices.  The I,J arguments can be left out if not needed.
See the :doc:`Howto output <Howto_output>` doc page for a discussion of
global, per-atom, and local data, and of scalar, vector, and array
data types.  See the doc pages for individual :doc:`computes <compute>`
and :doc:`fixes <fix>` for a description of what they calculate and
store.

For extract_variable(), an :doc:`equal-style or atom-style variable <variable>` is evaluated and its result returned.

For equal-style variables a single double value is returned and the
group argument is ignored.  For atom-style variables, a vector of
doubles is returned, one value per atom, which you can use via normal
Python subscripting. The values will be zero for atoms not in the
specified group.

The set_variable() method sets an existing string-style variable to a
new string value, so that subsequent LAMMPS commands can access the
variable.


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

   from ctypes import \*
   natoms = lmp.get_natoms()
   n3 = 3\*natoms
   x = (n3\*c_double)()
   x[0] = x coord of atom with ID 1
   x[1] = y coord of atom with ID 1
   x[2] = z coord of atom with ID 1
   x[3] = x coord of atom with ID 2
   ...
   x[n3-1] = z coord of atom with ID natoms
   lmp.scatter_atoms("x",1,3,x)

Alternatively, you can just change values in the vector returned by
the gather methods, since they are also ctypes vectors.

