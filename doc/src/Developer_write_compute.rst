Writing new compute styles
^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute styles are used to calculate global or per-atom scalar, vector,
or array quantities or local or per-chunk or grid based properties from
the current state of the simulation.  Examples include temperatures,
pressures, kinetic energies, and user-defined quantities.  All compute
styles are derived from the ``Compute`` base class, which is defined in
``src/compute.h``.  An overview of the available methods and flags is
given on the page for :doc:`modifying or extending compute styles
<Modify_compute>`.

Computes are executed *on demand* and thus require a "consumer" like
:doc:`fix ave/time <fix_ave_time>`, :doc:`compute reduce
<compute_reduce>`, or a :doc:`custom thermo style <thermo_style>` or
:doc:`custom dump style <dump>`.  This is different from
:doc:`Developer_write_fix` which are executed at :doc:`pre-determined
steps during a run <Developer_flow>` and either at every steps or with a
set frequency.  Some compute instances are also created and used
internally, e.g. for :doc:`thermo output <thermo>`.  It should be noted
that some computes do not actually perform computations (e.g.  computes
for potential energy or virial) but rather suitable flags are set and
then the corresponding data is collected by the force computation as
needed and then the compute merely makes that pre-collected data
available.  Such computes may not be executed between runs but only
during a run (or minimization).

In general, new compute styles should be added to the
:ref:`EXTRA-COMPUTE package <PKG-EXTRA-COMPUTE>`.  It is usually not
necessary to create accelerated compute styles except for the
:ref:`KOKKOS <PKG-KOKKOS>` package where there is a significant benefit
to avoid time consuming data transfers between GPU and host.  If you
feel that your contribution should be added to a different package,
please consult with the :doc:`LAMMPS developers <Intro_authors>` first.
The contributed code needs to support the :doc:`traditional GNU make
build process <Build_make>` **and** the :doc:`CMake build process
<Build_cmake>`.

----

Case 1: A global scalar compute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we describe how to implement a compute that produces a
single global scalar value.  As a concrete example we use :doc:`compute
temp <compute_temp>`, which computes the kinetic temperature of a group
of atoms.  The full implementation can be found in
``src/compute_temp.cpp`` and ``src/compute_temp.h``.

Header file
"""""""""""

Every compute style must be registered in LAMMPS by including the
following lines after the copyright block and before the include guards
for the class definition:

.. code-block:: c++

   #ifdef COMPUTE_CLASS
   // clang-format off
   ComputeStyle(temp,ComputeTemp);
   // clang-format on
   #else

The block between ``#ifdef COMPUTE_CLASS`` and ``#else`` will be
included by the ``Modify`` class in ``modify.cpp`` to build a map of
factory functions.  The map connects the style name ``temp`` with the
class name ``ComputeTemp``.  During compilation, LAMMPS constructs a
file ``style_compute.h`` that ``#include``\s all installed compute style
headers.  The ``// clang-format`` comments prevent the ``clang-format``
tool from incorrectly inserting spaces inside the macro arguments.

The class definition:

.. code-block:: c++

   #ifndef LMP_COMPUTE_TEMP_H
   #define LMP_COMPUTE_TEMP_H

   #include "compute.h"

   namespace LAMMPS_NS {

   class ComputeTemp : public Compute {
    public:
     ComputeTemp(class LAMMPS *, int, char **);
     ~ComputeTemp() override;
     void init() override {}
     void setup() override;
     double compute_scalar() override;
     void compute_vector() override;

    protected:
     double tfactor;

     virtual void dof_compute();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

The class derives from ``Compute`` and overrides the pure virtual method
``init()`` (here with an empty inline body since there is no
initialization needed for this example) as well as several optional
methods.  All overriding methods carry the ``override`` keyword.

Implementation file
"""""""""""""""""""

Constructor
'''''''''''

The constructor calls the base class constructor and sets the flags
that tell LAMMPS what kind of output this compute produces:

.. code-block:: c++

   ComputeTemp::ComputeTemp(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
   {
     if (narg != 3) error->all(FLERR, "Illegal compute temp command");

     scalar_flag = vector_flag = 1;
     size_vector = 6;
     extscalar = 0;
     extvector = 1;
     tempflag = 1;

     vector = new double[size_vector];
   }

The flags have the following meanings:

- ``scalar_flag = 1``: this compute provides a global scalar via
  ``compute_scalar()``.
- ``vector_flag = 1``: this compute also provides a global vector via
  ``compute_vector()``.
- ``size_vector = 6``: the vector has 6 components (the six components
  of the kinetic energy tensor).
- ``extscalar = 0``: the scalar is *intensive* (does not scale with the
  number of atoms since temperature is an intensive property).
- ``extvector = 1``: each vector component is *extensive* (since kinetic
  energy scales with the number of atoms).
- ``tempflag = 1``: marks this compute as providing a temperature, which
  is used by thermostat fixes to find a compatible temperature compute.

The ``vector`` array is allocated in the constructor; the base class
``Compute`` provides a pointer ``vector`` that derived classes should
point to their allocated storage.

The destructor frees the ``vector`` array if ``copymode`` is false. The
``copymode`` flag is set to ``0`` by default and set to ``1`` by the
``ComputeTempKokkos`` class.  This class is derived from ``ComputeTemp``
and replaces the storage assigned to ``vector`` with a Kokkos specific
storage object to facilitate convenient data transfer between host and
accelerator devices:

.. code-block:: c++

   ComputeTemp::~ComputeTemp()
   {
     if (!copymode) delete[] vector;
   }

The init() method (required)
'''''''''''''''''''''''''''''

The ``init()`` method must be overridden in every compute style; it is
pure virtual in the base class.  It is called once before each run or
minimization to perform one-time initialization (e.g., checking that all
required fixes or computes exist, requesting neighbor lists, etc.).
Here, ``init()`` has nothing to do and is given an empty body inline in
the header.

The setup() method (optional)
''''''''''''''''''''''''''''''

The optional ``setup()`` method is called at the start of each run,
after ``init()``.  For ``ComputeTemp``, it computes the number of degrees
of freedom:

.. code-block:: c++

   void ComputeTemp::setup()
   {
     dynamic = 0;
     if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
     dof_compute();
   }

The ``dynamic`` flag controls whether the number of degrees of freedom
is re-computed every timestep (needed when atoms enter or leave the
group or the simulation).

The compute_scalar() method (optional)
'''''''''''''''''''''''''''''''''''''''

The ``compute_scalar()`` method computes the temperature and returns it
as a double:

.. code-block:: c++

   double ComputeTemp::compute_scalar()
   {
     invoked_scalar = update->ntimestep;

     double **v = atom->v;
     double *mass = atom->mass;
     double *rmass = atom->rmass;
     int *type = atom->type;
     int *mask = atom->mask;
     int nlocal = atom->nlocal;

     double t = 0.0;

     if (rmass) {
       for (int i = 0; i < nlocal; i++)
         if (mask[i] & groupbit)
           t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
     } else {
       for (int i = 0; i < nlocal; i++)
         if (mask[i] & groupbit)
           t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
     }

     MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
     if (dynamic) dof_compute();
     if (dof < 0.0 && natoms_temp > 0.0)
       error->all(FLERR, "Temperature compute degrees of freedom < 0");
     scalar *= tfactor;
     return scalar;
   }

There are several important points:

- ``invoked_scalar`` is set to the current timestep number to allow
  other computes and fixes to check whether the result is fresh and
  avoid redundant re-computation.
- The group membership is checked via ``mask[i] & groupbit``.  The
  ``groupbit`` variable is defined in the base class ``Compute`` and
  corresponds to the group specified in the compute command.
- Atom data may be distributed across MPI ranks.  Only *local* atoms
  (``i < atom->nlocal``) are iterated.  ``MPI_Allreduce`` is used to
  sum the result across all ranks.
- The result is stored in ``scalar`` (declared in the base class) and
  also returned.
- The branch on ``rmass`` handles the case where atoms have individual
  masses (e.g. for granular simulations) versus per-type masses.

The compute_vector() method (optional)
'''''''''''''''''''''''''''''''''''''''

The ``compute_vector()`` method computes the kinetic energy tensor and
stores it in the pre-allocated ``vector`` array:

.. code-block:: c++

   void ComputeTemp::compute_vector()
   {
     invoked_vector = update->ntimestep;

     double **v = atom->v;
     double *mass = atom->mass;
     double *rmass = atom->rmass;
     int *type = atom->type;
     int *mask = atom->mask;
     int nlocal = atom->nlocal;

     double massone, t[6];
     for (int i = 0; i < 6; i++) t[i] = 0.0;

     for (int i = 0; i < nlocal; i++)
       if (mask[i] & groupbit) {
         if (rmass)
           massone = rmass[i];
         else
           massone = mass[type[i]];
         t[0] += massone * v[i][0] * v[i][0];
         t[1] += massone * v[i][1] * v[i][1];
         t[2] += massone * v[i][2] * v[i][2];
         t[3] += massone * v[i][0] * v[i][1];
         t[4] += massone * v[i][0] * v[i][2];
         t[5] += massone * v[i][1] * v[i][2];
       }

     MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
     for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
   }

The result must be stored in the ``vector`` array (declared in the base
class) and ``invoked_vector`` must be set to mark the result as
fresh.

----

Case 2: A per-atom compute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we describe a compute that produces a per-atom quantity.
As an example we use :doc:`compute ke/atom <compute_ke_atom>`, which
computes the kinetic energy of each atom.  The full implementation can be
found in ``src/compute_ke_atom.cpp`` and ``src/compute_ke_atom.h``.

Header file
"""""""""""

The registration block uses ``ComputeStyle(ke/atom,ComputeKEAtom)``.
The class definition:

.. code-block:: c++

   class ComputeKEAtom : public Compute {
    public:
     ComputeKEAtom(class LAMMPS *, int, char **);
     ~ComputeKEAtom() override;
     void init() override;
     void compute_peratom() override;
     double memory_usage() override;

    private:
     int nmax;
     double *ke;
   };

The compute overrides ``init()``, ``compute_peratom()``, and
``memory_usage()``.  The per-atom data array ``ke`` and its current
allocated size ``nmax`` are private members.

Constructor
'''''''''''

.. code-block:: c++

   ComputeKEAtom::ComputeKEAtom(LAMMPS *lmp, int narg, char **arg) :
       Compute(lmp, narg, arg), ke(nullptr)
   {
     if (narg != 3) error->all(FLERR, "Illegal compute ke/atom command");

     peratom_flag = 1;
     size_peratom_cols = 0;

     nmax = 0;
   }

The flags have the following meanings:

- ``peratom_flag = 1``: marks this compute as providing per-atom data
  via ``compute_peratom()``.
- ``size_peratom_cols = 0``: the per-atom data is a vector (one value
  per atom); set to a non-zero value for arrays with multiple values per
  atom.

The destructor must free ``ke``:

.. code-block:: c++

   ComputeKEAtom::~ComputeKEAtom()
   {
     memory->destroy(ke);
   }

The init() method (required)
'''''''''''''''''''''''''''''

.. code-block:: c++

   void ComputeKEAtom::init()
   {
     if (modify->get_compute_by_style(style).size() > 1)
       if (comm->me == 0) error->warning(FLERR, "More than one compute {}", style);
   }

This ``init()`` implementation emits a warning if more than one compute
of this style exists, since having duplicates wastes computation.  Most
compute styles implement at least a minimal ``init()`` that checks for
required fixes or computes and requests neighbor lists when needed.

The compute_peratom() method (optional)
''''''''''''''''''''''''''''''''''''''''

.. code-block:: c++

   void ComputeKEAtom::compute_peratom()
   {
     invoked_peratom = update->ntimestep;

     // grow ke array if necessary

     if (atom->nmax > nmax) {
       memory->destroy(ke);
       nmax = atom->nmax;
       memory->create(ke, nmax, "ke/atom:ke");
       vector_atom = ke;
     }

     // compute kinetic energy for each atom in group

     double mvv2e = force->mvv2e;
     double **v = atom->v;
     double *mass = atom->mass;
     double *rmass = atom->rmass;
     int *mask = atom->mask;
     int *type = atom->type;
     int nlocal = atom->nlocal;

     if (rmass)
       for (int i = 0; i < nlocal; i++) {
         if (mask[i] & groupbit) {
           ke[i] = 0.5 * mvv2e * rmass[i] *
               (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
         } else
           ke[i] = 0.0;
       }
     else
       for (int i = 0; i < nlocal; i++) {
         if (mask[i] & groupbit) {
           ke[i] = 0.5 * mvv2e * mass[type[i]] *
               (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
         } else
           ke[i] = 0.0;
       }
   }

There are several important points:

- ``invoked_peratom`` is set to the current timestep to mark the result
  as fresh.
- The ``ke`` array is grown if the number of owned+ghost atoms
  (``atom->nmax``) has increased since the last allocation.  The base
  class pointer ``vector_atom`` must be updated whenever ``ke`` is
  reallocated.  Setting ``vector_atom = ke`` after each reallocation
  ensures that callers always see the current storage.
- Atoms not in the group have their per-atom value set to 0.0.
- Only *local* atoms (``i < nlocal``) are processed; ghost atoms are
  not included.

The memory_usage() method (optional)
'''''''''''''''''''''''''''''''''''''

The ``memory_usage()`` method returns an estimate of the memory usage
in bytes.  This is used by the :doc:`info command <info>` and for
diagnostic output:

.. code-block:: c++

   double ComputeKEAtom::memory_usage()
   {
     double bytes = (double) nmax * sizeof(double);
     return bytes;
   }

----

Notes on flags and output types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Compute`` base class provides several flags and member variables
that control how the compute's output is used elsewhere in LAMMPS.
The most important ones for new compute styles are summarized below.

+-------------------------------+------------------------------------------------------+
| Flag                          | Meaning                                              |
+===============================+======================================================+
| ``scalar_flag``               | 1 if compute_scalar() is implemented                 |
+-------------------------------+------------------------------------------------------+
| ``vector_flag``               | 1 if compute_vector() is implemented                 |
+-------------------------------+------------------------------------------------------+
| ``array_flag``                | 1 if compute_array() is implemented                  |
+-------------------------------+------------------------------------------------------+
| ``peratom_flag``              | 1 if compute_peratom() is implemented                |
+-------------------------------+------------------------------------------------------+
| ``local_flag``                | 1 if compute_local() is implemented                  |
+-------------------------------+------------------------------------------------------+
| ``size_vector``               | length of the global vector                          |
+-------------------------------+------------------------------------------------------+
| ``size_array_rows``           | number of rows of the global array                   |
+-------------------------------+------------------------------------------------------+
| ``size_array_cols``           | number of columns of the global array                |
+-------------------------------+------------------------------------------------------+
| ``size_peratom_cols``         | number of columns in the per-atom array (0=vector)   |
+-------------------------------+------------------------------------------------------+
| ``extscalar``                 | 0 if scalar is intensive, 1 if extensive             |
+-------------------------------+------------------------------------------------------+
| ``extvector``                 | 0 if each vector component is intensive, 1 extensive |
+-------------------------------+------------------------------------------------------+
| ``extarray``                  | 0 if each array column is intensive, 1 extensive     |
+-------------------------------+------------------------------------------------------+
| ``tempflag``                  | 1 if compute produces a temperature                  |
+-------------------------------+------------------------------------------------------+
| ``pressflag``                 | 1 if compute produces a pressure                     |
+-------------------------------+------------------------------------------------------+
| ``vector``                    | pointer to the global vector output array            |
+-------------------------------+------------------------------------------------------+
| ``array``                     | pointer to the global 2D array output                |
+-------------------------------+------------------------------------------------------+
| ``vector_atom``               | pointer to the per-atom vector output array          |
+-------------------------------+------------------------------------------------------+
| ``array_atom``                | pointer to the per-atom 2D array output              |
+-------------------------------+------------------------------------------------------+

For the output arrays (``vector``, ``array``, ``vector_atom``,
``array_atom``), the derived class is responsible for allocating storage
and pointing the base class pointer to it.  The ``scalar`` variable is
declared in the base class and need not be separately allocated.

The ``invoked_scalar``, ``invoked_vector``, ``invoked_array``, and
``invoked_peratom`` variables should be set to ``update->ntimestep``
at the start of the corresponding compute method to allow other
routines to detect whether the output is current.
