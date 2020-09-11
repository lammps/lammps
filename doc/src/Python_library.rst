Python library interface
========================

As described previously, the Python interface to LAMMPS consists of a
Python "lammps" module, the source code for which is in
python/lammps.py, which creates a "lammps" object, with a set of
methods that can be invoked on that object.  The sample Python code
below assumes you have first imported the "lammps" module in your
Python script, as follows:

.. code-block:: Python

   from lammps import lammps

These are the methods defined by the lammps module.  If you look at
the files src/library.cpp and src/library.h you will see they
correspond one-to-one with calls you can make to the LAMMPS library
from a C++ or C or Fortran program, and which are described on the
:doc:`Howto library <Howto_library>` doc page.

The python/examples directory has Python scripts which show how Python
can run LAMMPS, grab data, change it, and put it back into LAMMPS.

.. code-block:: Python

   lmp = lammps()           # create a LAMMPS object using the default liblammps.so library
                            # 4 optional args are allowed: name, cmdargs, ptr, comm
   lmp = lammps(ptr=lmpptr) # use lmpptr as previously created LAMMPS object
   lmp = lammps(comm=split) # create a LAMMPS object with a custom communicator, requires mpi4py 2.0.0 or later
   lmp = lammps(name="g++")   # create a LAMMPS object using the liblammps_g++.so library
   lmp = lammps(name="g++",cmdargs=list)    # add LAMMPS command-line args, e.g. list = ["-echo","screen"]

   lmp.close()              # destroy a LAMMPS object

   version = lmp.version()  # return the numerical version id, e.g. LAMMPS 2 Sep 2015 -> 20150902

   lmp.file(file)           # run an entire input script, file = "in.lj"
   lmp.command(cmd)         # invoke a single LAMMPS command, cmd = "run 100"
   lmp.commands_list(cmdlist)     # invoke commands in cmdlist = **"run 10", "run 20"**
   lmp.commands_string(multicmd)  # invoke commands in multicmd = "run 10\nrun 20"

   size = lmp.extract_setting(name)     # return data type info

   xlo = lmp.extract_global(name,type)  # extract a global quantity
                                        # name = "boxxlo", "nlocal", etc
                                        # type = 0 = int
                                        #        1 = double

   boxlo,boxhi,xy,yz,xz,periodicity,box_change = lmp.extract_box()  # extract box info

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

   value = lmp.get_thermo(name)              # return current value of a thermo keyword
   natoms = lmp.get_natoms()                 # total # of atoms as int

   flag = lmp.set_variable(name,value)       # set existing named string-style variable to value, flag = 0 if successful
   lmp.reset_box(boxlo,boxhi,xy,yz,xz)       # reset the simulation box size

   data = lmp.gather_atoms(name,type,count)  # return per-atom property of all atoms gathered into data, ordered by atom ID
                                             # name = "x", "charge", "type", etc
   data = lmp.gather_atoms_concat(name,type,count)  # ditto, but concatenated atom values from each proc (unordered)
   data = lmp.gather_atoms_subset(name,type,count,ndata,ids)  # ditto, but for subset of Ndata atoms with IDs

   lmp.scatter_atoms(name,type,count,data)   # scatter per-atom property to all atoms from data, ordered by atom ID
                                             # name = "x", "charge", "type", etc
                                             # count = # of per-atom values, 1 or 3, etc

   lmp.scatter_atoms_subset(name,type,count,ndata,ids,data)  # ditto, but for subset of Ndata atoms with IDs

   lmp.create_atoms(n,ids,types,x,v,image,shrinkexceed)   # create N atoms with IDs, types, x, v, and image flags

----------

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

The file(), command(), commands_list(), commands_string() methods
allow an input script, a single command, or multiple commands to be
invoked.

The extract_setting(), extract_global(), extract_box(),
extract_atom(), extract_compute(), extract_fix(), and
extract_variable() methods return values or pointers to data
structures internal to LAMMPS.

For extract_global() see the src/library.cpp file for the list of
valid names.  New names could easily be added.  A double or integer is
returned.  You need to specify the appropriate data type via the type
argument.

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

The get_thermo() method returns the current value of a thermo
keyword as a float.

The get_natoms() method returns the total number of atoms in the
simulation, as an int.

The set_variable() method sets an existing string-style variable to a
new string value, so that subsequent LAMMPS commands can access the
variable.

The reset_box() method resets the size and shape of the simulation
box, e.g. as part of restoring a previously extracted and saved state
of a simulation.

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

----------

As noted above, these Python class methods correspond one-to-one with
the functions in the LAMMPS library interface in src/library.cpp and
library.h.  This means you can extend the Python wrapper via the
following steps:

* Add a new interface function to src/library.cpp and
  src/library.h.
* Rebuild LAMMPS as a shared library.
* Add a wrapper method to python/lammps.py for this interface
  function.
* You should now be able to invoke the new interface function from a
  Python script.

