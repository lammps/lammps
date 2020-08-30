LAMMPS Developer Guide
**********************

This section describes the internal structure and basic algorithms
of the LAMMPS code. This is a work in progress and additional
information will be added incrementally depending on availability
of time and requests from the LAMMPS user community.


LAMMPS source files
===================

The source files of the LAMMPS code are distributed across two
directories of the distribution.  The core of the code is located in the
``src`` folder and its sub-directories. Almost all of those are C++ files
(implementation files have a ``.cpp`` extension and and headers a
``.h``).  A sizable number of these files are in the ``src`` directory
itself, but there are plenty of :doc:`packages <Packages>`, which can be
included or excluded when LAMMPS is built.  See the :doc:`Include
packages in build <Build_package>` section of the manual for more
information about that part of the build process.  LAMMPS currently
supports building with :doc:`conventional makefiles <Build_make>` and
through :doc:`CMake <Build_cmake>` which differ in how packages are
enabled or disabled for a LAMMPS binary.  The source files for each
package are in all-uppercase sub-directories of the ``src`` folder, for
example ``src/MOLECULE`` or ``src/USER-MISC``.  The ``src/STUBS``
sub-directory is not a package but contains a dummy MPI library, that is
used when building a serial version of the code. the ``src/MAKE``
directory contains makefiles with settings and flags for a variety of
configuration and machines for the build process with traditional
makefiles.

The ``lib`` directory contains the source code for several supporting
libraries or files with configuration settings to use globally installed
libraries, that are required by some of the optional packages.
Each sub-directory, like ``lib/poems`` or ``lib/gpu``, contains the
source files, some of which are in different languages such as Fortran
or CUDA. These libraries are linked to during a LAMMPS build, if the
corresponding package is installed.

LAMMPS C++ source files almost always come in pairs, such as
``src/run.cpp`` and ``src/run.h``.  The pair of files defines a C++
class, for example the :cpp:class:`LAMMPS_NS::Run` class which contains
the code invoked by the :doc:`run <run>` command in a LAMMPS input script.
As this example illustrates, source file and class names often have a
one-to-one correspondence with a command used in a LAMMPS input script.
Some source files and classes do not have a corresponding input script
command, e.g. ``src/force.cpp`` and the :cpp:class:`LAMMPS_NS::Force`
class.  They are discussed in the next section.

LAMMPS class topology
=====================

Though LAMMPS has a lot of source files and classes, its class topology
is relative flat, as outlined in the :ref:`class-topology` figure.  Each
name refers to a class and has a pair of associated source files in the
``src`` folder, for example the class :cpp:class:`LAMMPS_NS::Memory`
corresponds to the files ``memory.cpp`` and ``memory.h``, or the class
:cpp:class:`LAMMPS_NS::AtomVec` corresponds to the files
``atom_vec.cpp`` and ``atom_vec.h``.  Full lines in the figure represent
compositing: that is the class to the left holds a pointer to an
instance of the class to the right.  Dashed lines instead represent
inheritance: the class to the right is derived from the class on the
left. Classes with a red boundary are not instantiated directly, but
they represent the base classes for "styles".  Those "styles" make up
the bulk of the LAMMPS code and only a few typical examples are included
in the figure for demonstration purposes.

.. _class-topology:
.. figure:: JPG/lammps-classes.png

   LAMMPS class topology

   This figure shows some of the relations of the base classes of the
   LAMMPS simulation package.  Full lines indicate that a class holds an
   instance of the class it is pointing to; dashed lines point to
   derived classes that are given as examples of what classes may be
   instantiated during a LAMMPS run based on the input commands and
   accessed through the API define by their respective base classes.  At
   the core is the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>` class, which
   holds pointers to class instances with specific purposes.  Those may
   hold instances of other classes, sometimes directly, or only
   temporarily, sometimes as derived classes or derived classes or
   derived classes, which may also hold instances of other classes.

The :cpp:class:`LAMMPS_NS::LAMMPS` class is the topmost class and
represents what is referred to an "instance" of LAMMPS.  It is a
composite holding references to instances of other core classes
providing the core functionality of the MD engine in LAMMPS and through
them abstractions of the required operations.  The constructor of the
LAMMPS class will instantiate those instances, process the command line
flags, initialize MPI (if not already done) and set up file pointers for
input and output. The destructor will shut everything down and free all
associated memory.  Thus code for the standalone LAMMPS executable in
``main.cpp`` simply initializes MPI, instantiates a single instance of
LAMMPS, and passes it the command line flags and input script. It
deletes the LAMMPS instance after the method reading the input returns
and shuts down the MPI environment before it exits the executable.

The :cpp:class:`LAMMPS_NS::Pointers` is not shown in the
:ref:`class-topology` figure, it holds references to members of the
`LAMMPS_NS::LAMMPS`, so that all classes derived from
:cpp:class:`LAMMPS_NS::Pointers` have direct access to those reference.
From the class topology all classes with blue boundary are referenced in
this class and all classes in the second and third columns, that are not
listed as derived classes are instead derived from
:cpp:class:`LAMMPS_NS::Pointers`.

Since all storage is encapsulated, the LAMMPS class can also be
instantiated multiple times by a calling code, and that can be either
simultaneously or consecutively.  When running in parallel with MPI,
care has to be taken, that suitable communicators are used to not
create conflicts between different instances.

The LAMMPS class currently holds instances of 19 classes representing
different core functionalities There are a handful of virtual parent
classes in LAMMPS that define what LAMMPS calls ``styles``.  They are
shaded red in the :ref:`class-topology` figure.  Each of these are
parents of a number of child classes that implement the interface
defined by the parent class.  There are two main categories of these
``styles``: some may only have one instance active at a time (e.g. atom,
pair, bond, angle, dihedral, improper, kspace, comm) and there is a
dedicated pointer variable in the composite class that manages them.
Setups that require a mix of different such styles have to use a
*hybrid* class that manages and forwards calls to the corresponding
sub-styles for the designated subset of atoms or data. or the composite
class may have lists of class instances, e.g. Modify handles lists of
compute and fix styles, while Output handles dumps class instances.

The exception to this scheme are the ``command`` style classes. These
implement specific commands that can be invoked before, after, or between
runs or are commands which launch a simulation.  For these an instance
of the class is created, its command() method called and then, after
completion, the class instance deleted.  Examples for this are the
create_box, create_atoms, minimize, run, or velocity command styles.

For all those ``styles`` certain naming conventions are employed: for
the fix nve command the class is called FixNVE and the files are
``fix_nve.h`` and ``fix_nve.cpp``. Similarly for fix ave/time we have
FixAveTime and ``fix_ave_time.h`` and ``fix_ave_time.cpp``. Style names
are lower case and without spaces or special characters. A suffix or
multiple appended with a forward slash '/' denotes a variant of the
corresponding class without the suffix. To connect the style name and
the class name, LAMMPS uses macros like the following ATOM\_CLASS,
PAIR\_CLASS, BOND\_CLASS, REGION\_CLASS, FIX\_CLASS, COMPUTE\_CLASS,
or DUMP\_CLASS in the corresponding header file.  During compilation
files with the pattern ``style_name.h`` are created that contain include
statements including all headers of all styles of a given type that
are currently active (or "installed).


More details on individual classes in the :ref:`class-topology` are as
follows:

- The Memory class handles allocation of all large vectors and arrays.

- The Error class prints all error and warning messages.

- The Universe class sets up partitions of processors so that multiple
  simulations can be run, each on a subset of the processors allocated
  for a run, e.g. by the mpirun command.

- The Input class reads and processes input input strings and files,
  stores variables, and invokes :doc:`commands <Commands_all>`.

- As discussed above, command style classes are directly derived from
  the Pointers class. They provide input script commands that perform
  one-time operations before/after/between simulations or which invoke a
  simulation.  They are instantiated from within the Input class,
  invoked, then immediately destructed.

- The Finish class is instantiated to print statistics to the screen
  after a simulation is performed, by commands like run and minimize.

- The Special class walks the bond topology of a molecular system to
  find first, second, third neighbors of each atom.  It is invoked by
  several commands, like :doc:`read_data <read_data>`,
  :doc:`read_restart <read_restart>`, or :doc:`replicate <replicate>`.

- The Atom class stores per-atom properties associated with atom styles.
  More precisely, they are allocated and managed by a class derived from
  the AtomVec class, and the Atom class simply stores pointers to them.
  The classes derived from AtomVec represent the different atom styles
  and they are instantiated through the :doc:`atom_style <atom_style>`
  command.

- The Update class holds instances of an integrator and a minimizer
  class.  The Integrate class is a parent style for the Verlet and
  r-RESPA time integrators, as defined by the :doc:`run_style
  <run_style>` command.  The Min class is a parent style for various
  energy minimizers.

- The Neighbor class builds and stores neighbor lists.  The NeighList
  class stores a single list (for all atoms).  A NeighRequest class
  instance is created by pair, fix, or compute styles when they need a
  particular kind of neighbor list and use the NeighRequest properties
  to select the neighbor list settings for the given request. There can
  be multiple instances of the NeighRequest class and the Neighbor class
  will try to optimize how they are computed by creating copies or
  sub-lists where possible.

- The Comm class performs inter-processor communication, typically of
  ghost atom information.  This usually involves MPI message exchanges
  with 6 neighboring processors in the 3d logical grid of processors
  mapped to the simulation box. There are two :doc:`communication styles
  <comm_style>` enabling different ways to do the domain decomposition.
  Sometimes the Irregular class is used, when atoms may migrate to
  arbitrary processors.

- The Domain class stores the simulation box geometry, as well as
  geometric Regions and any user definition of a Lattice.  The latter
  are defined by the :doc:`region <region>` and :doc:`lattice <lattice>`
  commands in an input script.

- The Force class computes various forces between atoms.  The Pair
  parent class is for non-bonded or pair-wise forces, which in LAMMPS
  also includes many-body forces such as the Tersoff 3-body potential if
  those are computed by walking pairwise neighbor lists.  The Bond,
  Angle, Dihedral, Improper parent classes are styles for bonded
  interactions within a static molecular topology.  The KSpace parent
  class is for computing long-range Coulombic interactions.  One of its
  child classes, PPPM, uses the FFT3D and Remap classes to redistribute
  and communicate grid-based information across the parallel processors.

- The Modify class stores lists of class instances derived from the
  :doc:`Fix <fix>` and :doc:`Compute <compute>` base classes.

- The Group class manipulates groups that atoms are assigned to via the
  :doc:`group <group>` command.  It also has functions to compute
  various attributes of groups of atoms.

- The Output class is used to generate 3 kinds of output from a LAMMPS
  simulation: thermodynamic information printed to the screen and log
  file, dump file snapshots, and restart files.  These correspond to the
  :doc:`Thermo <thermo_style>`, :doc:`Dump <dump>`, and
  :doc:`WriteRestart <write_restart>` classes respectively.  The Dump
  class is a base class with several derived classes implementing
  various dump style variants.

- The Timer class logs timing information, output at the end
  of a run.

.. TODO section on "Spatial decomposition and parallel operations"
..       diagram of 3d processor grid, brick vs. tiled. local vs. ghost
..       atoms, 6-way communication with pack/unpack functions,
..       PBC as part of the communication

.. TODO section on "Fixes, Computes, and Variables"
..      how and when data is computed and provided and how it is
..      referenced. flags in Fix/Compute/Variable classes tell
..      style and amount of available data.


How a timestep works
====================

The first and most fundamental operation within LAMMPS to understand is
how a timestep is structured.  Timestepping is performed by calling
methods of the Integrate class instance within the Update class.  Since
Integrate is a base class, it will point to an instance of a derived
class corresponding to what is selected by the :doc:`run_style
<run_style>` input script command.

In this section, the timestep implemented by the Verlet class is
described.  A similar timestep protocol is implemented by the Respa
class, for the r-RESPA hierarchical timestepping method.

The Min base class performs energy minimization, so does not perform a
literal timestep.  But it has logic similar to what is described here,
to compute forces and invoke fixes at each iteration of a minimization.
Differences between time integration and minimization are highlighted at
the end of this section.

The Verlet class is encoded in the ``src/verlet.cpp`` and ``verlet.h``
files.  It implements the velocity-Verlet timestepping algorithm.  The
workhorse method is ``Verlet::run()``, but first we highlight several
other methods in the class.

- The ``init()`` method is called at the beginning of each dynamics
  run.  It simply sets some internal flags, based on user settings in
  other parts of the code.

- The ``setup()`` or ``setup_minimal()`` methods are also called before
  each run.  The velocity-Verlet method requires current forces be
  calculated before the first timestep, so these routines compute
  forces due to all atomic interactions, using the same logic that
  appears in the timestepping described next.  A few fixes are also
  invoked, using the mechanism described in the next section.  Various
  counters are also initialized before the run begins.  The
  ``setup_minimal()`` method is a variant that has a flag for performing
  less setup.  This is used when runs are continued and information
  from the previous run is still valid.  For example, if repeated
  short LAMMPS runs are being invoked, interleaved by other commands,
  via the *pre no* and *every* options of the run command, the
  ``setup_minimal()`` method is used.

- The ``force_clear()`` method initializes force and other arrays to
  zero before each timestep, so that forces (torques, etc) can be
  accumulated.

Now for the ``Verlet::run()`` method.  Its basic structure in hi-level pseudo
code is shown below.  In the actual code in ``src/verlet.cpp`` some of
these operations are conditionally invoked.

.. code-block:: python

   loop over N timesteps:
     if timeout condition: break
     ev_set()

     fix->initial_integrate()
     fix->post_integrate()

     nflag = neighbor->decide()
     if nflag:
       fix->pre_exchange()
       domain->pbc()
       domain->reset_box()
       comm->setup()
       neighbor->setup_bins()
       comm->exchange()
       comm->borders()
       fix->pre_neighbor()
       neighbor->build()
       fix->post_neighbor()
     else:
       comm->forward_comm()

     force_clear()
     fix->pre_force()

     pair->compute()
     bond->compute()
     angle->compute()
     dihedral->compute()
     improper->compute()
     kspace->compute()

     fix->pre_reverse()
     comm->reverse_comm()

     fix->post_force()
     fix->final_integrate()
     fix->end_of_step()

     if any output on this step:
       output->write()

   # after loop
   fix->post_run()


The ``ev_set()`` method (in the parent Integrate class), sets two flags
(*eflag* and *vflag*) for energy and virial computation.  Each flag
encodes whether global and/or per-atom energy and virial should be
calculated on this timestep, because some fix or variable or output will
need it.  These flags are passed to the various methods that compute
particle interactions, so that they either compute and tally the
corresponding data or can skip the extra calculations if the energy and
virial are not needed.  See the comments for the ``Integrate::ev_set()``
method which document the flag values.

At various points of the timestep, fixes are invoked,
e.g. ``fix->initial_integrate()``.  In the code, this is actually done
via the Modify class which stores all the Fix objects and lists of which
should be invoked at what point in the timestep.  Fixes are the LAMMPS
mechanism for tailoring the operations of a timestep for a particular
simulation.  As described elsewhere, each fix has one or more methods,
each of which is invoked at a specific stage of the timestep, as show in
the timestep pseudo-code.  All the active fixes defined in an input
script, that are flagged to have an ``initial_integrate()`` method are
invoked at the beginning of each timestep.  Examples are :doc:`fix nve
<fix_nve>` or :doc:`fix nvt or fix npt <fix_nh>` which perform the
start-of-timestep velocity-Verlet integration operations to update
velocities by a half-step, and coordinates by a full step.  The
``post_integrate()`` method is next for operations that need to happen
immediately after those updates.  Only a few fixes use this, e.g. to
reflect particles off box boundaries in the :doc:`FixWallReflect class
<fix_wall_reflect>`.

The ``decide()`` method in the Neighbor class determines whether
neighbor lists need to be rebuilt on the current timestep (conditions
can be changed using the :doc:`neigh_modify every/delay/check
<neigh_modify>` command.  If not, coordinates of ghost atoms are
acquired by each processor via the ``forward_comm()`` method of the Comm
class.  If neighbor lists need to be built, several operations within
the inner if clause of the pseudo-code are first invoked.  The
``pre_exchange()`` method of any defined fixes is invoked first.
Typically this inserts or deletes particles from the system.

Periodic boundary conditions are then applied by the Domain class via
its ``pbc()`` method to remap particles that have moved outside the
simulation box back into the box.  Note that this is not done every
timestep, but only when neighbor lists are rebuilt.  This is so that
each processor's sub-domain will have consistent (nearby) atom
coordinates for its owned and ghost atoms.  It is also why dumped atom
coordinates may be slightly outside the simulation box if not dumped
on a step where the neighbor lists are rebuilt.

The box boundaries are then reset (if needed) via the ``reset_box()``
method of the Domain class, e.g. if box boundaries are shrink-wrapped to
current particle coordinates.  A change in the box size or shape
requires internal information for communicating ghost atoms (Comm class)
and neighbor list bins (Neighbor class) be updated.  The ``setup()``
method of the Comm class and ``setup_bins()`` method of the Neighbor
class perform the update.

The code is now ready to migrate atoms that have left a processor's
geometric sub-domain to new processors.  The ``exchange()`` method of
the Comm class performs this operation.  The ``borders()`` method of the
Comm class then identifies ghost atoms surrounding each processor's
sub-domain and communicates ghost atom information to neighboring
processors.  It does this by looping over all the atoms owned by a
processor to make lists of those to send to each neighbor processor.  On
subsequent timesteps, the lists are used by the ``Comm::forward_comm()``
method.

Fixes with a ``pre_neighbor()`` method are then called.  These typically
re-build some data structure stored by the fix that depends on the
current atoms owned by each processor.

Now that each processor has a current list of its owned and ghost
atoms, LAMMPS is ready to rebuild neighbor lists via the ``build()``
method of the Neighbor class.  This is typically done by binning all
owned and ghost atoms, and scanning a stencil of bins around each
owned atom's bin to make a Verlet list of neighboring atoms within the
force cutoff plus neighbor skin distance.

In the next portion of the timestep, all interaction forces between
particles are computed, after zeroing the per-atom force vector via the
``force_clear()`` method.  If the newton flag is set to *on* by the
newton command, forces are added to both owned and ghost atoms, otherwise
only to owned (aka local) atoms.

Pairwise forces are calculated first, which enables the global virial
(if requested) to be calculated cheaply (at O(N) cost instead of O(N**2)
at the end of the ``Pair::compute()`` method), by a dot product of atom
coordinates and forces.  By including owned and ghost atoms in the dot
product, the effect of periodic boundary conditions is correctly
accounted for.  Molecular topology interactions (bonds, angles,
dihedrals, impropers) are calculated next (if supported by the current
atom style).  The final contribution is from long-range Coulombic
interactions, invoked by the KSpace class.

The ``pre_reverse()`` method in fixes is used for operations that have to
be done *before* the upcoming reverse communication (e.g. to perform
additional data transfers or reductions for data computed during the
force computation and stored with ghost atoms).

If the newton flag is on, forces on ghost atoms are communicated and
summed back to their corresponding owned atoms.  The ``reverse_comm()``
method of the Comm class performs this operation, which is essentially
the inverse operation of sending copies of owned atom coordinates to
other processor's ghost atoms.

At this point in the timestep, the total force on each (local) atom is
known.  Additional force constraints (external forces, SHAKE, etc) are
applied by Fixes that have a ``post_force()`` method.  The second half
of the velocity-Verlet integration, ``final_integrate()`` is then
performed (another half-step update of the velocities) via fixes like
nve, nvt, npt.

At the end of the timestep, fixes that contain an ``end_of_step()``
method are invoked.  These typically perform a diagnostic calculation,
e.g. the ave/time and ave/spatial fixes.  The final operation of the
timestep is to perform any requested output, via the ``write()`` method
of the Output class.  There are 3 kinds of LAMMPS output: thermodynamic
output to the screen and log file, snapshots of atom data to a dump
file, and restart files.  See the :doc:`thermo_style <thermo_style>`,
:doc:`dump <dump>`, and :doc:`restart <restart>` commands for more
details.

The the flow of control during energy minimization iterations is
similar to that of a molecular dynamics timestep.  Forces are computed,
neighbor lists are built as needed, atoms migrate to new processors, and
atom coordinates and forces are communicated to neighboring processors.
The only difference is what Fix class operations are invoked when.  Only
a subset of LAMMPS fixes are useful during energy minimization, as
explained in their individual doc pages.  The relevant Fix class methods
are ``min_pre_exchange()``, ``min_pre_force()``, and ``min_post_force()``.
Each fix is invoked at the appropriate place within the minimization
iteration.  For example, the ``min_post_force()`` method is analogous to
the ``post_force()`` method for dynamics; it is used to alter or constrain
forces on each atom, which affects the minimization procedure.

After all iterations are completed there is a ``cleanup`` step which
calls the ``post_run()`` method of fixes to perform operations only required
at the end of a calculations (like freeing temporary storage or creating
final outputs).

Writing LAMMPS styles
=====================

The :doc:`Modify` section of the manual gives an overview of how LAMMPS can
be extended by writing new classes that derive from existing
parent classes in LAMMPS.  Here, some specific coding
details are provided for writing code for LAMMPS.

Writing a new fix style
^^^^^^^^^^^^^^^^^^^^^^^

Writing fixes is a flexible way of extending LAMMPS.  Users can
implement many things using fixes:

- changing particles attributes (positions, velocities, forces, etc.). Examples: FixNVE, FixFreeze.
- reading/writing data. Example: FixRestart.
- adding or modifying properties due to geometry. Example: FixWall.
- interacting with other subsystems or external code: Examples: FixTTM, FixExternal, FixLATTE
- saving information for analysis or future use (previous positions,
  for instance). Examples: Fix AveTime, FixStoreState.


All fixes are derived from the Fix base class and must have a
constructor with the signature: ``FixPrintVel(class LAMMPS *, int, char **)``.

Every fix must be registered in LAMMPS by writing the following lines
of code in the header before include guards:

.. code-block:: c

   #ifdef FIX_CLASS
   FixStyle(print/vel,FixPrintVel)
   #else
   /* the definition of the FixPrintVel class comes here */
   ...
   #endif

Where ``print/vel`` is the style name of your fix in the input script and
``FixPrintVel`` is the name of the class. The header file would be called
``fix_print_vel.h`` and the implementation file ``fix_print_vel.cpp``.
These conventions allow LAMMPS to automatically integrate it into the
executable when compiling and associate your new fix class with the designated
keyword when it parses the input script.

Let's write a simple fix which will print the average velocity at the end
of each timestep. First of all, implement a constructor:

.. code-block:: C++

   FixPrintVel::FixPrintVel(LAMMPS *lmp, int narg, char **arg)
   : Fix(lmp, narg, arg)
   {
     if (narg < 4)
       error->all(FLERR,"Illegal fix print/vel command");

     nevery = force->inumeric(FLERR,arg[3]);
     if (nevery <= 0)
       error->all(FLERR,"Illegal fix print/vel command");
   }

In the constructor you should parse your fix arguments which are
specified in the script. All fixes have pretty the same syntax:
``fix <fix-ID> <fix group> <fix name> <fix arguments ...>``. The
first 3 parameters are parsed by Fix base class constructor, while
``<fix arguments>`` should be parsed by you. In our case, we need to
specify how often we want to print an average velocity. For instance,
once in 50 timesteps: ``fix 1 print/vel 50``. There is a special variable
in the Fix class called ``nevery`` which specifies how often the method
``end_of_step()`` is called. Thus all we need to do is just set it up.

The next method we need to implement is ``setmask()``:

.. code-block:: C++

   int FixPrintVel::setmask()
   {
     int mask = 0;
     mask |= FixConst::END_OF_STEP;
     return mask;
   }

Here the user specifies which methods of your fix should be called
during execution. The constant ``END_OF_STEP`` corresponds to the
``end_of_step()`` method. The most important available methods that
are called during a timestep and the order in which they are called
are shown in the previous section.

.. code-block:: C++

   void FixPrintVel::end_of_step()
   {
     // for add3, scale3
     using namespace MathExtra;

     double** v = atom->v;
     int nlocal = atom->nlocal;
     double localAvgVel[4]; // 4th element for particles count
     memset(localAvgVel, 0, 4 * sizeof(double));
     for (int particleInd = 0; particleInd < nlocal; ++particleInd) {
       add3(localAvgVel, v[particleInd], localAvgVel);
     }
     localAvgVel[3] = nlocal;
     double globalAvgVel[4];
     memset(globalAvgVel, 0, 4 * sizeof(double));
     MPI_Allreduce(localAvgVel, globalAvgVel, 4, MPI_DOUBLE, MPI_SUM, world);
     scale3(1.0 / globalAvgVel[3], globalAvgVel);
     if ((comm->me == 0) && screen) {
       fmt::print(screen,"{}, {}, {}\n",
                  globalAvgVel[0], globalAvgVel[1], globalAvgVel[2]);
     }
   }

In the code above, we use MathExtra routines defined in
``math_extra.h``.  There are bunch of math functions to work with
arrays of doubles as with math vectors.  It is also important to note
that LAMMPS code should always assume to be run in parallel and that
atom data is thus distributed across the MPI ranks.  Thus you can
only process data from local atoms directly and need to use MPI library
calls to combine or exchange data.  For serial execution, LAMMPS
comes bundled with the MPI STUBS library that contains the MPI library
function calls in dummy versions that only work for a single MPI rank.

In this code we use an instance of Atom class. This object is stored
in the Pointers class (see ``pointers.h``) which is the base class of
the Fix base class. This object contains references to various class
instances (the original instances are created and held by the LAMMPS
class) with all global information about the simulation system.
Data from the Pointers class is available to all classes inherited from
it using protected inheritance. Hence when you write you own class,
which is going to use LAMMPS data, don't forget to inherit from Pointers
or pass an Pointer to it to all functions that need access. When writing
fixes we inherit from class Fix which is inherited from Pointers so
there is no need to inherit from it directly.

The code above computes average velocity for all particles in the
simulation.  Yet you have one unused parameter in fix call from the
script: ``group_name``.  This parameter specifies the group of atoms
used in the fix. So we should compute average for all particles in the
simulation only if ``group_name == "all"``, but it can be any group.
The group membership information of an atom is contained in the *mask*
property of and atom and the bit corresponding to a given group is
stored in the groupbit variable which is defined in Fix base class:

.. code-block:: C++

   for (int i = 0; i < nlocal; ++i) {
     if (atom->mask[i] & groupbit) {
     // Do all job here
     }
   }

Class Atom encapsulates atoms positions, velocities, forces, etc. User
can access them using particle index. Note, that particle indexes are
usually changed every few timesteps because of neighbor list rebuilds
and spatial sorting (to improve cache efficiency).

Let us consider another Fix example: We want to have a fix which stores
atoms position from previous time step in your fix. The local atoms
indexes may not be valid on the next iteration. In order to handle
this situation there are several methods which should be implemented:

- ``double memory_usage()``: return how much memory the fix uses (optional)
- ``void grow_arrays(int)``: do reallocation of the per particle arrays in your fix
- ``void copy_arrays(int i, int j, int delflag)``: copy i-th per-particle
  information to j-th. Used when atom sorting is performed. if delflag is set
  and atom j owns a body, move the body information to atom i.
- ``void set_arrays(int i)``: sets i-th particle related information to zero

Note, that if your class implements these methods, it must call add calls of
add_callback and delete_callback to constructor and destructor. Since we want
to store positions of atoms from previous timestep, we need to add
``double** xold`` to the header file. Than add allocation code
to the constructor:

.. code-block:: C++

   FixSavePos::FixSavePos(LAMMPS *lmp, int narg, char **arg), xold(nullptr)
   {
   //...
     memory->create(xold, atom->nmax, 3, "FixSavePos:x");
     atom->add_callback(0);
   }

   FixSavePos::~FixSavePos() {
     atom->delete_callback(id, 0);
     memory->destroy(xold);
   }

Implement the aforementioned methods:

.. code-block:: C++

   double FixSavePos::memory_usage()
   {
     int nmax = atom->nmax;
     double bytes = 0.0;
     bytes += nmax * 3 * sizeof(double);
     return bytes;
   }

   void FixSavePos::grow_arrays(int nmax)
   {
     memory->grow(xold, nmax, 3, "FixSavePos:xold");
   }

   void FixSavePos::copy_arrays(int i, int j, int delflag)
   {
     memcpy(xold[j], xold[i], sizeof(double) * 3);
   }

   void FixSavePos::set_arrays(int i)
   {
     memset(xold[i], 0, sizeof(double) * 3);
   }

   int FixSavePos::pack_exchange(int i, double *buf)
   {
     int m = 0;
     buf[m++] = xold[i][0];
     buf[m++] = xold[i][1];
     buf[m++] = xold[i][2];

     return m;
   }

   int FixSavePos::unpack_exchange(int nlocal, double *buf)
   {
     int m = 0;
     xold[nlocal][0] = buf[m++];
     xold[nlocal][1] = buf[m++];
     xold[nlocal][2] = buf[m++];

     return m;
   }

Now, a little bit about memory allocation. We use the Memory class which
is just a bunch of template functions for allocating 1D and 2D
arrays. So you need to add include ``memory.h`` to have access to them.

Finally, if you need to write/read some global information used in
your fix to the restart file, you might do it by setting flag
``restart_global = 1`` in the constructor and implementing methods void
``write_restart(FILE *fp)`` and ``void restart(char *buf)``.
If, in addition, you want to write the per-atom property to restart
files additional settings and functions are needed:

- a fix flag indicating this needs to be set ``restart_peratom = 1;``
- ``atom->add_callback()`` and ``atom->delete_callback()`` must be called
  a second time with the final argument set to 1 instead of 0 (indicating
  restart processing instead of per-atom data memory management).
- the functions ``void pack_restart(int i, double *buf)`` and
  ``void unpack_restart(int nlocal, int nth)`` need to be implemented

---------------------------

LAMMPS utility functions
========================

The ``utils`` sub-namespace inside the ``LAMMPS_NS`` namespace provides
a collection of convenience functions and utilities that perform common
tasks that are required repeatedly throughout the LAMMPS code like
reading or writing to files with error checking or translation of
strings into specific types of numbers with checking for validity.  This
reduces redundant implementations and encourages consistent behavior.

I/O functions with validity check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are wrappers around the corresponding C library calls like
``fgets()`` or ``fread()``.  They will check if there were errors
on reading or an unexpected end-of-file state was reached.  In that
case, the functions will stop the calculation with an error message,
indicating the name of the problematic file, if possible.

.. doxygenfunction:: sfgets
   :project: progguide

.. doxygenfunction:: sfread
   :project: progguide

String to number conversions with validity check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are functions equivalent to those in the ``Force`` class that
were implemented with the aim to replace and supersede those.  Unlike
the versions in ``Force``, these can be used in cases where only
a single MPI rank is trying to convert strings to numbers, as you
can select through an argument, whether ``Error->all()`` or ``Error->one()``
should be called on improper strings.

These functions are preferred over C library calls like ``atoi()`` or
``atof()`` since they check if the **entire** provided string is a valid
(floating-point or integer) number, and will error out instead of silently
return the result of a partial conversion or zero in cases where the
string is not a valid number.  This behavior allows to more easily detect
typos or issues when processing input files.

.. doxygenfunction:: numeric
   :project: progguide

.. doxygenfunction:: inumeric
   :project: progguide

.. doxygenfunction:: bnumeric
   :project: progguide

.. doxygenfunction:: tnumeric
   :project: progguide


String processing functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following are functions to help with processing strings
and parsing files or arguments.

.. doxygenfunction:: trim
   :project: progguide

.. doxygenfunction:: trim_comment
   :project: progguide

.. doxygenfunction:: count_words(const char *text)
   :project: progguide

.. doxygenfunction:: count_words(const std::string &text)
   :project: progguide

.. doxygenfunction:: count_words(const std::string &text, const std::string &separators)
   :project: progguide

.. doxygenfunction:: trim_and_count_words
   :project: progguide

.. doxygenfunction:: split_words
   :project: progguide

.. doxygenfunction:: strmatch
   :project: progguide

.. doxygenfunction:: is_integer
   :project: progguide

.. doxygenfunction:: is_double
   :project: progguide

Filename and path functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: guesspath
   :project: progguide

.. doxygenfunction:: path_basename
   :project: progguide

.. doxygenfunction:: path_join
   :project: progguide

.. doxygenfunction:: file_is_readable
   :project: progguide

Potential file functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: get_potential_file_path
   :project: progguide

.. doxygenfunction:: get_potential_date
   :project: progguide

.. doxygenfunction:: get_potential_units
   :project: progguide

.. doxygenfunction:: get_supported_conversions
   :project: progguide

.. doxygenfunction:: get_conversion_factor
   :project: progguide

Convenience functions
^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: logmesg
   :project: progguide

.. doxygenfunction:: getsyserror
   :project: progguide

.. doxygenfunction:: check_packages_for_style
   :project: progguide

.. doxygenfunction:: timespec2seconds
   :project: progguide

---------------------------

Tokenizer classes
=================

The purpose of the tokenizer classes is to simplify the recurring task
of breaking lines of text down into words and/or numbers.
Traditionally, LAMMPS code would be using the ``strtok()`` function from
the C library for that purpose, but that function has two significant
disadvantages: 1) it cannot be used concurrently from different LAMMPS
instances since it stores its status in a global variable and 2) it
modifies the string that it is processing.  These classes were
implemented to avoid both of these issues and also to reduce the amount
of code that needs to be written.

The basic procedure is to create an instance of the class with the
string to be processed as an argument and then do a loop until all
available tokens are read.  The constructor has a default set of
separator characters, but that can be overridden. The default separators
are all "whitespace" characters, i.e. the space character, the tabulator
character, the carriage return character, the linefeed character, and
the form feed character.

.. doxygenclass:: LAMMPS_NS::Tokenizer
   :project: progguide

.. doxygenclass:: LAMMPS_NS::ValueTokenizer
   :project: progguide
