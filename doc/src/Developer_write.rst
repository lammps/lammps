Writing new styles
------------------

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
specified in the script. All fixes have pretty much the same syntax:
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

