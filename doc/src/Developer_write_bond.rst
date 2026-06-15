Writing new bond, angle, dihedral, and improper styles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bond, angle, dihedral, and improper styles are used to compute bonded
interactions between groups of two, three, or four atoms whose topology
is described in a data file or set by other commands.  They are derived
from the ``Bond``, ``Angle``, ``Dihedral``, and ``Improper`` base
classes, respectively.  As shown on the corresponding pages :doc:`Bond,
angle, dihedral, improper styles <Modify_bond>`, the classes that
compute these molecular interactions share a common design.

In general, new styles for bonded interactions should be added to the
:ref:`EXTRA-MOLECULE package <PKG-EXTRA-MOLECULE>`.  If you feel that
your contribution should be added to a different package (e.g. the
:ref:`MOLECULE package <PKG-MOLECULE>`), please consult with the
:doc:`LAMMPS developers <Intro_authors>` first.

Before implementing an accelerated style the corresponding plain CPU
version should be implemented and properly tested.  The accelerated
version should then be derived from the plain CPU version so that only
the code relevant for acceleration is re-implemented.  Those should then
be added to the corresponding accelerator package (:ref:`KOKKOS
<PKG-KOKKOS>`, :ref:`OPENMP <PKG-OPENMP>`, :ref:`INTEL <PKG-INTEL>`).

The contributed code needs to support the :doc:`traditional GNU make
build process <Build_make>` **and** the :doc:`CMake build process
<Build_cmake>`.  This is usually automatic unless the new style uses
some external library, which is uncommon for bonded interactions.

It is highly recommended also to include one or more example input decks
and :doc:`force-style unit tests <Developer_unittest>`.  The tests will
be particularly useful to test consistency between the ``compute()`` and
``single()`` methods, plain and accelerated versions, and to test
correct restarting and writing and reading of data files.

----

Case 1: Implementing a bond style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we describe the procedure of implementing a new bond
style.  As a concrete example we use the harmonic bond style
:doc:`bond_style harmonic <bond_harmonic>`, the simplest non-trivial
bond style.  The implementation can be found in the files
``src/MOLECULE/bond_harmonic.cpp`` and ``src/MOLECULE/bond_harmonic.h``
of the LAMMPS source code.  It implements the potential:

.. math::

   E = K (r - r_0)^2

with the spring constant :math:`K` and the equilibrium distance
:math:`r_0` as per-type coefficients.

Header file
"""""""""""

The header file ``bond_harmonic.h`` starts with the standard LAMMPS
copyright and license block (see, for example, the discussion in
:doc:`Writing a new pair style <Developer_write_pair>` for details).

After the copyright block, every bond style must be registered in LAMMPS
by including the following lines before the include guards:

.. code-block:: c++

   #ifdef BOND_CLASS
   // clang-format off
   BondStyle(harmonic,BondHarmonic);
   // clang-format on
   #else

This block between ``#ifdef BOND_CLASS`` and ``#else`` will be included
by the ``Force`` class in ``force.cpp`` to build a map of factory
functions for bond styles.  The map connects the name of the bond style,
"harmonic", with the name of the class, ``BondHarmonic``.  During
compilation, LAMMPS generates a file ``style_bond.h`` that contains
``#include`` statements for all installed bond styles.  Before including
``style_bond.h`` into ``force.cpp``, the ``BOND_CLASS`` define is set
and the ``BondStyle(name,class)`` macro defined.  The ``//
clang-format`` comments prevent ``clang-format`` from reformatting the
macro argument in a way that would break it.

Analogously, an angle style would use ``#ifdef ANGLE_CLASS`` and
``AngleStyle(name,class)``, a dihedral style would use ``#ifdef
DIHEDRAL_CLASS`` and ``DihedralStyle(name,class)``, and an improper
style would use ``#ifdef IMPROPER_CLASS`` and
``ImproperStyle(name,class)``.

The class definition itself follows after the include guard:

.. code-block:: c++

   #ifndef LMP_BOND_HARMONIC_H
   #define LMP_BOND_HARMONIC_H

   #include "bond.h"

   namespace LAMMPS_NS {

   class BondHarmonic : public Bond {
    public:
     BondHarmonic(class LAMMPS *);
     ~BondHarmonic() override;
     void compute(int, int) override;
     void coeff(int, char **) override;
     double equilibrium_distance(int) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_data(FILE *) override;
     double single(int, double, int, int, double &) override;
     void born_matrix(int, double, int, int, double &, double &) override;
     void *extract(const char *, int &) override;

    protected:
     double *k, *r0;

     virtual void allocate();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

The class is derived from ``Bond`` and overrides all *pure* virtual
methods as require by the C++ standard (``compute()``, ``coeff()``,
``equilibrium_distance()``, ``write_restart()``, ``read_restart()``, and
``single()``) and several optional ones.  All overriding methods use the
``override`` keyword to allow the compiler to detect mismatches.  The
per-type coefficient arrays ``k`` and ``r0`` are declared as protected
member variables.  The ``allocate()`` method is declared as ``virtual``
so that derived classes can override it if needed.

.. note::

   Ideally, there should be no additional ``#include`` statements
   outside of ``#include "bond.h"``.  Where possible, forward
   declarations should be used and the proper include statements only
   added to the implementation file (discussed below).  Exceptions are
   headers for C++ container classes.  Please see
   :doc:`Modify_requirements` and :doc:`Modify_style` for more
   information and this and related issues.

Implementation file
"""""""""""""""""""

The implementation file ``bond_harmonic.cpp`` starts with the same
copyright block and includes, followed by ``using namespace
LAMMPS_NS;``.

Constructor and destructor
''''''''''''''''''''''''''

The constructor is minimal; it only calls the base class constructor
and sets the ``born_matrix_enable`` flag to enable the optional
``born_matrix()`` method:

.. code-block:: c++

   BondHarmonic::BondHarmonic(LAMMPS *_lmp) : Bond(_lmp), k(nullptr), r0(nullptr)
   {
     born_matrix_enable = 1;
   }

It is recommended to always initialize all pointers to ``nullptr`` in
the initializer list so that they can be safely deleted in the
destructor.

The destructor frees the per-type coefficient arrays, but only if
``allocated`` is true and ``copymode`` is false.  The ``copymode``
flag is set by Kokkos-accelerated styles that copy the base class
rather than owning the storage:

.. code-block:: c++

   BondHarmonic::~BondHarmonic()
   {
     if (allocated && !copymode) {
       memory->destroy(setflag);
       memory->destroy(k);
       memory->destroy(r0);
     }
   }

The allocate() method
'''''''''''''''''''''

The ``allocate()`` helper allocates all per-type arrays.  Types are
one-indexed, so arrays of size ``nbondtypes + 1`` are created. The
``setflag`` array, which tracks whether coefficients have been set for
each type, is also allocated here and initialized to zero:

.. code-block:: c++

   void BondHarmonic::allocate()
   {
     allocated = 1;
     const int np1 = atom->nbondtypes + 1;

     memory->create(k, np1, "bond:k");
     memory->create(r0, np1, "bond:r0");

     memory->create(setflag, np1, "bond:setflag");
     for (int i = 1; i < np1; i++) setflag[i] = 0;
   }

The ``setflag`` array is declared in the ``Bond`` base class.  It is
used by ``Bond::init()`` to check that coefficients have been set for
all bond types.

The coeff() method (required)
'''''''''''''''''''''''''''''

The ``coeff()`` method processes the arguments of the
:doc:`bond_coeff <bond_coeff>` command.  The first argument is always
the type range (e.g. ``1*3``, handled by ``utils::bounds()``);
subsequent arguments are the potential coefficients.

.. code-block:: c++

   void BondHarmonic::coeff(int narg, char **arg)
   {
     if (narg != 3) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
     if (!allocated) allocate();

     int ilo, ihi;
     utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

     double k_one = utils::numeric(FLERR, arg[1], false, lmp);
     double r0_one = utils::numeric(FLERR, arg[2], false, lmp);

     int count = 0;
     for (int i = ilo; i <= ihi; i++) {
       k[i] = k_one;
       r0[i] = r0_one;
       setflag[i] = 1;
       count++;
     }

     if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
   }

The ``utils::bounds()`` function converts the type string into integer
lower and upper bounds ``ilo`` and ``ihi``.  After setting the
coefficients for each type in the range, ``setflag[i]`` is set to 1.
The final check guards against a range that matches no types.

The compute() method (required)
''''''''''''''''''''''''''''''''

The ``compute()`` method iterates over all bonds in the local bond
list (``neighbor->bondlist``) and accumulates forces and energies.
The two arguments ``eflag`` and ``vflag`` control whether energies and
virial contributions are computed.  The actual energy and virial
accumulation is handled by calling ``ev_tally()`` at the end of each
iteration step.

.. code-block:: c++

   void BondHarmonic::compute(int eflag, int vflag)
   {
     int i1, i2, n, type;
     double delx, dely, delz, ebond, fbond;
     double rsq, r, dr, rk;

     ebond = 0.0;
     ev_init(eflag, vflag);

     double **x = atom->x;
     double **f = atom->f;
     int **bondlist = neighbor->bondlist;
     int nbondlist = neighbor->nbondlist;
     int nlocal = atom->nlocal;
     int newton_bond = force->newton_bond;

     for (n = 0; n < nbondlist; n++) {
       i1 = bondlist[n][0];
       i2 = bondlist[n][1];
       type = bondlist[n][2];

       delx = x[i1][0] - x[i2][0];
       dely = x[i1][1] - x[i2][1];
       delz = x[i1][2] - x[i2][2];

       rsq = delx * delx + dely * dely + delz * delz;
       r = sqrt(rsq);
       dr = r - r0[type];
       rk = k[type] * dr;

       // force & energy

       if (r > 0.0)
         fbond = -2.0 * rk / r;
       else
         fbond = 0.0;

       if (eflag) ebond = rk * dr;

       // apply force to each of 2 atoms

       if (newton_bond || i1 < nlocal) {
         f[i1][0] += delx * fbond;
         f[i1][1] += dely * fbond;
         f[i1][2] += delz * fbond;
       }

       if (newton_bond || i2 < nlocal) {
         f[i2][0] -= delx * fbond;
         f[i2][1] -= dely * fbond;
         f[i2][2] -= delz * fbond;
       }

       if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
     }
   }

Before the loop, ``ev_init(eflag, vflag)`` initializes the energy and
virial accumulators and sets several flags (including ``evflag``) that
control which quantities are computed.  The force scalar ``fbond``
is defined such that the force on atom ``i1`` is
``fbond * (x[i1] - x[i2])`` and the force on atom ``i2`` is
``-fbond * (x[i1] - x[i2])``.  This sign convention matches the
force as the negative gradient of the energy.

The Newton's third law check (``newton_bond || i < nlocal``) is
important for correctness in parallel runs: when ``newton_bond`` is
on, forces are computed for both atoms even when they are on different
MPI ranks, relying on a subsequent reverse communication to add up the
contributions; when ``newton_bond`` is off, only forces on local atoms
(``i < nlocal``) are accumulated.

The equilibrium_distance() method (required for bond styles)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The ``equilibrium_distance()`` method returns the equilibrium bond
length for the given type.  This is used, for instance, by the :doc:`fix
shake <fix_shake>` and :doc:`fix rattle <fix_shake>` commands:

.. code-block:: c++

   double BondHarmonic::equilibrium_distance(int i)
   {
     return r0[i];
   }

The single() method (required for bond styles)
''''''''''''''''''''''''''''''''''''''''''''''

The ``single()`` method computes the force divided by *r* and the
potential energy for a single pair of atoms.  It is called by
:doc:`compute bond/local <compute_bond_local>` and related commands:

.. code-block:: c++

   double BondHarmonic::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
   {
     double r = sqrt(rsq);
     double dr = r - r0[type];
     double rk = k[type] * dr;
     fforce = 0;
     if (r > 0.0) fforce = -2.0 * rk / r;
     return rk * dr;
   }

The return value is the potential energy.  ``fforce`` is set to ``-dE/dr
/ r`` (i.e. the force divided by the interatomic distance) so that it
can be easily converted into the x-, y-, and z-direction force
components by multiplying with :math:`\Delta x`, :math:`\Delta y`, and
:math:`\Delta z`, respectively.

Restart and data file methods (required)
''''''''''''''''''''''''''''''''''''''''

The ``write_restart()`` method is called by proc 0 to write the
per-type coefficients to the binary restart file.  The matching
``read_restart()`` method is called when reading a restart file; it
allocates storage, reads coefficients on proc 0, and broadcasts them
to all ranks.

.. code-block:: c++

   void BondHarmonic::write_restart(FILE *fp)
   {
     fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
     fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
   }

   void BondHarmonic::read_restart(FILE *fp)
   {
     allocate();

     if (comm->me == 0) {
       utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
       utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
     }
     MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
     MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

     for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
   }

   void BondHarmonic::write_data(FILE *fp)
   {
     for (int i = 1; i <= atom->nbondtypes; i++) fprintf(fp, "%d %g %g\n", i, k[i], r0[i]);
   }

Note that index 0 of each array is unused (types are one-indexed), so
``&k[1]`` is the start of the actual data.  Using ``utils::sfread()``
(instead of ``fread()``) provides error checking.

----

Case 2: Implementing an angle style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Angle styles are similar to bond styles but involve three atoms.  The
implementation follows the same pattern as for bond styles, so for
details not mentioned here, you can refer to the corresponding section
documenting the implementation of bond styles.  As an example we use
:doc:`angle_style harmonic <angle_harmonic>` which implements the
potential:

.. math::

   E = K (\theta - \theta_0)^2

with spring constant :math:`K` and equilibrium angle :math:`\theta_0`.

The header file uses ``#ifdef ANGLE_CLASS`` and ``AngleStyle(name,class)``
for the registration block, and the class is derived from ``Angle``:

.. code-block:: c++

   #ifdef ANGLE_CLASS
   // clang-format off
   AngleStyle(harmonic,AngleHarmonic);
   // clang-format on
   #else

   #ifndef LMP_ANGLE_HARMONIC_H
   #define LMP_ANGLE_HARMONIC_H

   #include "angle.h"

   namespace LAMMPS_NS {

   class AngleHarmonic : public Angle {
    public:
     AngleHarmonic(class LAMMPS *);
     ~AngleHarmonic() override;
     void compute(int, int) override;
     void coeff(int, char **) override;
     double equilibrium_angle(int) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_data(FILE *) override;
     double single(int, int, int, int) override;
     void *extract(const char *, int &) override;

    protected:
     double *k, *theta0;

     virtual void allocate();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

Key differences from a bond style:

- The registration macro is ``AngleStyle`` and the include guard
  ``#ifdef ANGLE_CLASS``.
- The base class is ``Angle`` (``#include "angle.h"``).
- The ``equilibrium_distance()`` method is replaced by
  ``equilibrium_angle()`` (used by :doc:`fix shake <fix_shake>`).
- The ``single()`` method signature takes four atom indices instead
  of a distance squared plus two atom indices.
- The ``compute()`` loop iterates over ``neighbor->anglelist`` and
  ``neighbor->nanglelist`` involving three-atom groups.
- The types array is ``atom->nangletypes``.

The ``compute()`` method for angle styles uses three atom indices
(``i1``, ``i2``, ``i3``) from ``neighbor->anglelist[n][0..3]``, where
the type is at index 3.  The equilibrium angle ``theta0`` is stored
in radians (converted from degrees in ``coeff()``).

----

Case 3: Implementing a dihedral style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dihedral styles compute interactions among four atoms forming a
torsion angle.  Harmonic dihedrals implement the potential:

.. math::

   E = K \left[1 + d \cos(n\phi) \right]

where :math:`K` is the force constant, :math:`d` is the sign (+1 or
-1), and :math:`n` is the multiplicity (periodicity).  The implementation
can be found in ``src/MOLECULE/dihedral_harmonic.cpp`` and
``src/MOLECULE/dihedral_harmonic.h``.

The header file uses ``#ifdef DIHEDRAL_CLASS`` and
``DihedralStyle(name,class)``, and the class is derived from
``Dihedral``:

.. code-block:: c++

   #ifdef DIHEDRAL_CLASS
   // clang-format off
   DihedralStyle(harmonic,DihedralHarmonic);
   // clang-format on
   #else

   #ifndef LMP_DIHEDRAL_HARMONIC_H
   #define LMP_DIHEDRAL_HARMONIC_H

   #include "dihedral.h"

   namespace LAMMPS_NS {

   class DihedralHarmonic : public Dihedral {
    public:
     DihedralHarmonic(class LAMMPS *);
     ~DihedralHarmonic() override;
     void compute(int, int) override;
     void coeff(int, char **) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_data(FILE *) override;

    protected:
     double *k, *cos_shift, *sin_shift;
     int *sign, *multiplicity;

     virtual void allocate();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

Key differences from bond and angle styles:

- The base class is ``Dihedral`` (``#include "dihedral.h"``).
- There is no ``equilibrium_distance()`` or ``equilibrium_angle()``
  method; dihedral styles have no analogue required by other
  LAMMPS subsystems.
- There is no ``single()`` method in the base class for dihedrals.
- The ``compute()`` loop uses four atom indices from
  ``neighbor->dihedrallist[n][0..4]`` (type at index 4).
- The types array is ``atom->ndihedraltypes``.
- The ``writedata`` flag must be set to 1 in the constructor if
  ``write_data()`` is implemented (the default is 0 for dihedral styles).

In the constructor, the ``writedata`` flag is set:

.. code-block:: c++

   DihedralHarmonic::DihedralHarmonic(LAMMPS *_lmp) :
       Dihedral(_lmp), k(nullptr), cos_shift(nullptr), sin_shift(nullptr), sign(nullptr),
       multiplicity(nullptr)
   {
     writedata = 1;
   }

----

Case 4: Implementing an improper style
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Improper styles compute interactions involving four atoms where the
geometry represents an out-of-plane deformation rather than a
conventional torsion.  Harmonic impropers implement:

.. math::

   E = K (\chi - \chi_0)^2

where :math:`K` is the force constant and :math:`\chi_0` is the
equilibrium improper angle.

The header file uses ``#ifdef IMPROPER_CLASS`` and
``ImproperStyle(name,class)``, and the class is derived from ``Improper``:

.. code-block:: c++

   #ifdef IMPROPER_CLASS
   // clang-format off
   ImproperStyle(harmonic,ImproperHarmonic);
   // clang-format on
   #else

   #ifndef LMP_IMPROPER_HARMONIC_H
   #define LMP_IMPROPER_HARMONIC_H

   #include "improper.h"

   namespace LAMMPS_NS {

   class ImproperHarmonic : public Improper {
    public:
     ImproperHarmonic(class LAMMPS *);
     ~ImproperHarmonic() override;
     void compute(int, int) override;
     void coeff(int, char **) override;
     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_data(FILE *) override;
     void *extract(const char *, int &) override;

    protected:
     double *k, *chi;

     virtual void allocate();
   };

   }    // namespace LAMMPS_NS
   #endif
   #endif

Key differences from dihedral styles:

- The base class is ``Improper`` (``#include "improper.h"``).
- The ``compute()`` loop uses four atom indices from
  ``neighbor->improperlist[n][0..4]`` (type at index 4).
- The types array is ``atom->nimpropertypes``.
- A ``symmatoms`` array in the base class can be used to record
  which atom in the quadruplet is the central atom of symmetry.
  In ``ImproperHarmonic``, ``symmatoms[0] = 1`` indicates that the
  first atom in the quadruplet is the atom of symmetry.

Constructor setup:

.. code-block:: c++

   ImproperHarmonic::ImproperHarmonic(LAMMPS *_lmp) : Improper(_lmp), k(nullptr), chi(nullptr)
   {
     writedata = 1;
     // the first atom in the quadruplet is the atom of symmetry
     symmatoms[0] = 1;
   }
