Writing new pair styles
^^^^^^^^^^^^^^^^^^^^^^^

Pair styles are at the core of most simulations with LAMMPS, since they
are used to compute the forces (plus energy and virial contributions, if
needed) on atoms for pairs or small clusters of atoms within a given
cutoff.  This is often the dominant computation in LAMMPS, and sometimes
even the only one.  Pair styles can be grouped into multiple categories:

#. simple pairwise additive interactions of point particles
   (e.g. :doc:`Lennard-Jones <pair_lj>`, :doc:`Morse <pair_morse>`,
   :doc:`Buckingham <pair_buck>`)
#. pairwise additive interactions of point particles with added
   :doc:`Coulomb <pair_coul>` interactions or *only* the Coulomb interactions
#. manybody interactions of point particles (e.g. :doc:`EAM <pair_eam>`,
   :doc:`Tersoff <pair_tersoff>`)
#. complex interactions that include additional per-atom properties
   (e.g. Discrete Element Models (DEM), Peridynamics, Ellipsoids)
#. special purpose pair styles that may not even compute forces like
   :doc:`pair_style zero <pair_zero>` and :doc:`pair_style tracker
   <pair_tracker>`, or are a wrapper for multiple kinds of interactions
   like :doc:`pair_style hybrid <pair_hybrid>`, :doc:`pair_style list <pair_list>`,
   and :doc:`pair_style kim <pair_kim>`

In the text below, we will discuss aspects of implementing pair styles
in LAMMPS by looking at representative case studies.  The design of
LAMMPS allows developers to focus on the essentials, which is to compute
the forces (and energies or virial contributions), enter and manage the
global settings as well as the potential parameters, and the pair style
specific parts of reading and writing restart and data files.  Most of
the complex tasks like management of the atom positions, domain
decomposition and boundaries, or neighbor list creation are handled
transparently by other parts of the LAMMPS code.

As shown on the page for :doc:`writing or extending pair styles
<Modify_pair>`, in order to implement a new pair style, a new class must
be written that is either directly or indirectly derived from the
``Pair`` class.  If that class is directly derived from ``Pair``, there
are three methods that *must* be re-implemented, since they are "pure"
in the base class: ``Pair::compute()``, ``Pair::settings()``,
``Pair::coeff()``.  In addition, a custom constructor is needed.  All
other methods are optional and have default implementations in the base
class (most of which do nothing), but they may need to be overridden
depending on the requirements of the model.

We are looking at the following cases:

- `Case 1: a pairwise additive model`_
- `Case 2: a many-body potential`_
- `Case 3: a potential requiring communication`_
- `Case 4: potentials without a compute() function`_

----

Case 1: a pairwise additive model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, we will describe the procedure of adding a simple pair
style to LAMMPS: an empirical model that can be used to model liquid
mercury.  The pair style shall be called :doc:`bond/gauss
<pair_born_gauss>` and the complete implementation can be found in the
files ``src/EXTRA-PAIR/pair_born_gauss.cpp`` and
``src/EXTRA-PAIR/pair_born_gauss.h`` of the LAMMPS source code.

Model and general considerations
""""""""""""""""""""""""""""""""

The functional form of the model according to :ref:`(Bomont) <Bomont2>`
consists of a repulsive Born-Mayer exponential term and a temperature
dependent, attractive Gaussian term.

.. math::

   E = A_0 \exp \left( -\alpha r \right) - A_1 \exp\left[ -\beta \left(r - r_0 \right)^2 \right]

For the application to mercury, the following parameters are listed:

- :math:`A_0 = 8.2464 \times 10^{13} \; \textrm{eV}`
- :math:`\alpha = 12.48 \; \AA^{-1}`
- :math:`\beta = 0.44 \; \AA^{-2}`
- :math:`r_0 = 3.56 \; \AA`
- :math:`A_1` is temperature dependent and can be determined from
  :math:`A_1 = a_0 + a_1 T + a_2 T^2` with:

  - :math:`a_0 = 1.97475 \times 10^{-2} \; \textrm{eV}`
  - :math:`a_1 = 8.40841 \times 10^{-5} \; \textrm{eV/K}`
  - :math:`a_2 = -2.58717 \times 10^{-8} \; \textrm{eV/K}^{-2}`

With the optional cutoff, this means we have a total of 5 or 6
parameters for each pair of atom types. Additionally, we need to input a
default cutoff value as a global setting.

Because of the combination of Born-Mayer with a Gaussian, the pair style
shall be named "born/gauss" and thus the class name would be
``PairBornGauss`` and the source files ``pair_born_gauss.h`` and
``pair_born_gauss.cpp``.  Since this is a rather uncommon potential, it
shall be added to the :ref:`EXTRA-PAIR <PKG-EXTRA-PAIR>` package.

Header file
"""""""""""

The first segment of any LAMMPS source should be the copyright and
license statement.  Note the marker in the first line to indicate to
editors like emacs that this file is a C++ source, even though the .h
extension suggests a C source (this is a convention inherited from the
very beginning of the C++ version of LAMMPS).

.. code-block:: c++

   /* -*- c++ -*- ----------------------------------------------------------
      LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
      https://www.lammps.org/, Sandia National Laboratories
      LAMMPS development team: developers@lammps.org

      Copyright (2003) Sandia Corporation.  Under the terms of Contract
      DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
      certain rights in this software.  This software is distributed under
      the GNU General Public License.

      See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

Every pair style must be registered in LAMMPS by including the following
lines of code in the second part of the header after the copyright
message and before the include guards for the class definition:

.. code-block:: c++

   #ifdef PAIR_CLASS
   // clang-format off
   PairStyle(born/gauss,PairBornGauss);
   // clang-format on
   #else

   /* the definition of the PairBornGauss class (see below) is inserted here */

   #endif

This block of between ``#ifdef PAIR_CLASS`` and ``#else`` will be
included by the ``Force`` class in ``force.cpp`` to build a map of
"factory functions" that will create an instance of these classes and
return a pointer to it.  The map connects the name of the pair style,
"born/gauss", to the name of the class, ``PairBornGauss``.  During
compilation, LAMMPS constructs a file ``style_pair.h`` that contains
``#include`` statements for all "installed" pair styles.  Before
including ``style_pair.h`` into ``force.cpp``, the ``PAIR_CLASS`` define
is set and the ``PairStyle(name,class)`` macro defined.  The code of the
macro adds the installed pair styles to the "factory map" which enables
the :doc:`pair_style command <pair_style>` to create the pair style
instance.

The list of header files to include is automatically updated by the
build system if there are new files, so the presence of the new header
file in the ``src/EXTRA-PAIR`` folder and the enabling of the EXTRA-PAIR
package will trigger LAMMPS to include the new pair style when it is
(re-)compiled.  The "// clang-format" format comments are needed so that
running :ref:`clang-format <clang-format>` on the file will not insert
unwanted blanks between "born", "/", and "gauss" which would break the
``PairStyle`` macro.

The third part of the header file is the actual class definition of the
``PairBornGauss`` class.  This has the prototypes for all member
functions that will be implemented by this pair style.  This includes
:doc:`a few required and a number of optional functions <Modify_pair>`.
All functions that were labeled in the base class as "virtual" must be
given the "override" property, as it is done in the code shown below.

The "override" property helps to detect unexpected mismatches because
compilation will stop with an error in case the signature of a function
is changed in the base class without also changing it in all derived
classes.  For example, if this change added an optional argument with a
default value, then all existing source code *calling* the function
would not need changes and still compile, but the function in the
derived class would no longer override the one in the base class due to
the different number of arguments and the behavior of the pair style is
thus changed in an unintended way.  Using the "override" keyword
prevents such issues.

.. code-block:: c++

   #ifndef LMP_PAIR_BORN_GAUSS_H
   #define LMP_PAIR_BORN_GAUSS_H

   #include "pair.h"

   namespace LAMMPS_NS {

   class PairBornGauss : public Pair {
    public:
     PairBornGauss(class LAMMPS *);
     ~PairBornGauss() override;

     void compute(int, int) override;
     void settings(int, char **) override;
     void coeff(int, char **) override;
     double init_one(int, int) override;

     void write_restart(FILE *) override;
     void read_restart(FILE *) override;
     void write_restart_settings(FILE *) override;
     void read_restart_settings(FILE *) override;
     void write_data(FILE *) override;
     void write_data_all(FILE *) override;

     double single(int, int, int, int, double, double, double, double &) override;
     void *extract(const char *, int &) override;

Also, variables and arrays for storing global settings and potential
parameters are defined.  Since these are internal to the class, they are
placed after a "protected:" label.

.. code-block:: c++

    protected:
     double cut_global;
     double **cut;
     double **biga0, **alpha, **biga1, **beta, **r0;
     double **a0, **a1, **a2;
     double **offset;

     virtual void allocate();
   };
   }    // namespace LAMMPS_NS
   #endif

Implementation file
"""""""""""""""""""

We move on to the implementation of the ``PairBornGauss`` class in the
``pair_born_gauss.cpp`` file.  This file also starts with a LAMMPS
copyright and license header.  Below that notice is typically the space
where comments may be added with additional information about this
specific file, the author(s), affiliation(s), and email address(es).
This way the contributing author(s) can be easily contacted, when
there are questions about the implementation later.  Since the file(s)
may be around for a long time, it is beneficial to use some kind of
"permanent" email address, if possible.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
      https://www.lammps.org/, Sandia National Laboratories
      LAMMPS development team: developers@lammps.org

      Copyright (2003) Sandia Corporation.  Under the terms of Contract
      DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
      certain rights in this software.  This software is distributed under
      the GNU General Public License.

      See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

   // Contributing author: Axel Kohlmeyer, Temple University, akohlmey@gmail.com

   #include "pair_born_gauss.h"

   #include "atom.h"
   #include "comm.h"
   #include "error.h"
   #include "fix.h"
   #include "force.h"
   #include "memory.h"
   #include "neigh_list.h"

   #include <cmath>
   #include <cstring>

   using namespace LAMMPS_NS;

The second section of the implementation file has various include
statements.  The include file for the class header has to come first,
then a block of LAMMPS classes (sorted alphabetically) followed by a
block of system headers and others, if needed.  Note the standardized
C++ notation for headers of C-library functions (``cmath`` instead of
``math.h``).  The final statement of this segment imports the
``LAMMPS_NS::`` namespace globally for this file.  This way, all LAMMPS
specific functions and classes do not have to be prefixed with
``LAMMPS_NS::``.

Constructor and destructor (required)
"""""""""""""""""""""""""""""""""""""

The first two functions in the implementation source file are typically
the constructor and the destructor.

Pair styles are different from most classes in LAMMPS that define a
"style", as their constructor only uses the LAMMPS class instance
pointer as an argument, but **not** the command line arguments of the
:doc:`pair_style command <pair_style>`.  Instead, those arguments are
processed in the ``Pair::settings()`` function (or rather the version in
the derived class).  The constructor is the place where global defaults
are set and specifically flags are set indicating which optional
features of a pair style are available.

.. code-block:: c++

   /* ---------------------------------------------------------------------- */

   PairBornGauss::PairBornGauss(LAMMPS *lmp) : Pair(lmp)
   {
     writedata = 1;
   }

The `writedata = 1;` statement indicates that the pair style is capable
of writing the current pair coefficient parameters to data files.  That
is, the class implements specific versions for ``Pair::data()`` and
``Pair::data_all()``.  Other statements that could be added here would
be `single_enable = 1;` or `respa_enable = 0;` to indicate that the
``Pair::single()`` function is present and the
``Pair::compute_(inner|middle|outer)`` functions are not, but those are
also the default settings and already set in the base class.

In the destructor, we need to delete all memory that was allocated by the
pair style, usually to hold force field parameters that were entered
with the :doc:`pair_coeff command <pair_coeff>`.  Most of those array
pointers will need to be declared in the derived class header, but some
(e.g. setflag, cutsq) are already declared in the base class.

.. code-block:: c++

   PairBornGauss::~PairBornGauss()
   {
     if (allocated) {
       memory->destroy(setflag);
       memory->destroy(cutsq);
       memory->destroy(cut);
       memory->destroy(biga0);
       memory->destroy(alpha);
       memory->destroy(biga1);
       memory->destroy(beta);
       memory->destroy(r0);
       memory->destroy(offset);
     }
   }


Settings and coefficients (required)
""""""""""""""""""""""""""""""""""""

To enter the global pair style settings and the pair style parameters,
the functions ``Pair::settings()`` and ``Pair::coeff()`` need to be
re-implemented.  The arguments to the ``settings()`` function are the
arguments given to the :doc:`pair_style command <pair_style>`.
Normally, those would already be processed as part of the constructor,
but moving this to a separate function allows users to change global
settings like the default cutoff without having to reissue all
pair_coeff commands or re-read the ``Pair Coeffs`` sections from the
data file.  In the ``settings()`` function, also the arrays for storing
parameters, to define cutoffs, track with pairs of parameters have been
explicitly set are allocated and, if needed, initialized.  In this case,
the memory allocation and initialization is moved to a function
``allocate()``.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      allocate all arrays
   ------------------------------------------------------------------------- */

   void PairBornGauss::allocate()
   {
     allocated = 1;
     int np1 = atom->ntypes + 1;

     memory->create(setflag, np1, np1, "pair:setflag");
     for (int i = 1; i < np1; i++)
       for (int j = i; j < np1; j++) setflag[i][j] = 0;

     memory->create(cutsq, np1, np1, "pair:cutsq");
     memory->create(cut, np1, np1, "pair:cut");
     memory->create(biga0, np1, np1, "pair:biga0");
     memory->create(alpha, np1, np1, "pair:alpha");
     memory->create(biga1, np1, np1, "pair:biga1");
     memory->create(beta, np1, np1, "pair:beta");
     memory->create(r0, np1, np1, "pair:r0");
     memory->create(offset, np1, np1, "pair:offset");
   }

   /* ----------------------------------------------------------------------
      global settings
   ------------------------------------------------------------------------- */

   void PairBornGauss::settings(int narg, char **arg)
   {
     if (narg != 1) error->all(FLERR, "Pair style bond/gauss must have exactly one argument");
     cut_global = utils::numeric(FLERR, arg[0], false, lmp);

     // reset per-type pair cutoffs that have been explicitly set previously

     if (allocated) {
       for (int i = 1; i <= atom->ntypes; i++)
         for (int j = i; j <= atom->ntypes; j++)
           if (setflag[i][j]) cut[i][j] = cut_global;
     }
   }

The arguments to the ``coeff()`` function are the arguments to the
:doc:`pair_coeff command <pair_coeff>`.  The function is also called
when processing the ``Pair Coeffs`` or ``PairIJ Coeffs`` sections of
data files.  In the case of the ``Pair Coeffs`` section, there is only
one atom type per line and thus the first argument is duplicated.  Since
the atom type arguments of the :doc:`pair_coeff command <pair_coeff>`
may be a range (e.g. \*\ 3 for atom types 1, 2, and 3), the
corresponding arguments are passed to the :cpp:func:`utils::bounds()
<LAMMPS_NS::utils::bounds>` function which will then return the low
and high end of the range.  Note that the ``setflag`` array is set to 1
for all pairs of atom types processed by this call.  This information is
later used in the ``init_one()`` function to determine if any coefficients
are missing and, if supported by the potential, generate those missing
coefficients from the selected mixing rule.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      set coeffs for one or more type pairs
   ------------------------------------------------------------------------- */

   void PairBornGauss::coeff(int narg, char **arg)
   {
     if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
     if (!allocated) allocate();

     int ilo, ihi, jlo, jhi;
     utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
     utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

     double biga0_one = utils::numeric(FLERR, arg[2], false, lmp);
     double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
     double biga1_one = utils::numeric(FLERR, arg[4], false, lmp);
     double beta_one = utils::numeric(FLERR, arg[5], false, lmp);
     double r0_one = utils::numeric(FLERR, arg[6], false, lmp);
     double cut_one = cut_global;
     if (narg == 10) cut_one = utils::numeric(FLERR, arg[7], false, lmp);

     int count = 0;
     for (int i = ilo; i <= ihi; i++) {
       for (int j = MAX(jlo, i); j <= jhi; j++) {
         biga0[i][j] = biga0_one;
         alpha[i][j] = alpha_one;
         biga1[i][j] = biga1_one;
         beta[i][j] = beta_one;
         r0[i][j] = r0_one;
         cut[i][j] = cut_one;
         setflag[i][j] = 1;
         count++;
       }
     }

     if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
   }

Initialization
""""""""""""""

The ``init_one()`` function is called during the :doc:`"init" phase
<Developer_flow>` of a simulation.  This is where potential parameters
are checked for completeness, derived parameters computed (e.g. the
"offset" of the potential energy at the cutoff distance for use with the
:doc:`pair_modify shift yes <pair_modify>` command).  If a pair style
supports generating "mixed" parameters (i.e. where both atoms of a pair
have a different atom type) using a "mixing rule" from the parameters of
the type with itself, this is the place to compute and store those mixed
values.  The *born/gauss* pair style does not support mixing, so we only
check for completeness.  Another purpose of the ``init_one()`` function
is to symmetrize the potential parameter arrays.  The return value of
the function is the cutoff for the given pair of atom types.  This
information is used by the neighbor list code to determine the largest
cutoff and then build the neighbor lists accordingly.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      init for one type pair i,j and corresponding j,i
   ------------------------------------------------------------------------- */

   double PairBornGauss::init_one(int i, int j)
   {
     if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

     if (offset_flag) {
       double dr = cut[i][j] - r0[i][j];
       offset[i][j] =
           biga0[i][j] * exp(-alpha[i][j] * cut[i][j]) - biga1[i][j] * exp(-beta[i][j] * dr * dr);
     } else
       offset[i][j] = 0.0;

     biga0[j][i] = biga0[i][j];
     alpha[j][i] = alpha[i][j];
     biga1[j][i] = biga1[i][j];
     beta[j][i] = beta[i][j];
     r0[j][i] = r0[i][j];
     offset[j][i] = offset[i][j];

     return cut[i][j];
   }


Computing forces from the neighbor list (required)
""""""""""""""""""""""""""""""""""""""""""""""""""

The ``compute()`` function is the "workhorse" of a pair style.  This is
where we have the nested loops over all pairs of particles from the
neighbor list to compute forces and - if needed - energies and virials.

The first part is to define some variables for later use and store
cached copies of data or pointers that we need to access frequently.  Also,
this is a good place to call ``Pair::ev_init()``, which initializes
several flags derived from the `eflag` and `vflag` parameters signaling
whether the energy and virial need to be tallied and whether only globally
or also per-atom.

.. code-block:: c++

   /* ---------------------------------------------------------------------- */

   void PairBornGauss::compute(int eflag, int vflag)
   {
     int i, j, ii, jj, inum, jnum, itype, jtype;
     double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
     double rsq, r, dr, aexp, bexp, factor_lj;
     int *ilist, *jlist, *numneigh, **firstneigh;

     evdwl = 0.0;
     ev_init(eflag, vflag);

     double **x = atom->x;
     double **f = atom->f;
     int *type = atom->type;
     int nlocal = atom->nlocal;
     double *special_lj = force->special_lj;
     int newton_pair = force->newton_pair;

     inum = list->inum;
     ilist = list->ilist;
     numneigh = list->numneigh;
     firstneigh = list->firstneigh;

The outer loop (index *i*) is over local atoms of our sub-domain.
Typically, the value of `inum` (the number of neighbor lists) is the
same as the number of local atoms (= atoms *owned* by this sub-domain).
But when the pair style is used as a sub-style of a :doc:`hybrid pair
style <pair_hybrid>` or neighbor list entries are removed with
:doc:`neigh_modify exclude <neigh_modify>`, this number may be
smaller. The array ``list->ilist`` has the (local) indices of the atoms
for which neighbor lists have been created. Then ``list->numneigh`` is
an `inum` sized array with the number of entries of each list of
neighbors, and ``list->firstneigh`` is a list of pointers to those lists.

For efficiency reasons, cached copies of some properties of the outer
loop atoms are also initialized.

.. code-block:: c++

     // loop over neighbors of my atoms

     for (ii = 0; ii < inum; ii++) {
       i = ilist[ii];
       xtmp = x[i][0];
       ytmp = x[i][1];
       ztmp = x[i][2];
       itype = type[i];
       jlist = firstneigh[i];
       jnum = numneigh[i];

The inner loop (index *j*) processes the neighbor lists.  The neighbor
list code encodes in the upper 2 bits whether a pair is a regular pair
of neighbor (= 0) or a pair of 1-2 (= 1), 1-3 (= 2), or 1-4 (= 3)
:doc:`"special" neighbor <special_bonds>`.  The ``sbmask()`` inline
function extracts those bits and converts them into a number.  This
number is used to look up the corresponding scaling factor for the
non-bonded interaction from the ``force->special_lj`` array and stores
it in the `factor_lj` variable.  Due to the additional bits, the value
of *j* would be out of range when accessing data from per-atom arrays,
so we apply the NEIGHMASK constant with a bit-wise and operation to mask
them out.  This step *must* be done, even if a pair style does not use
special bond scaling of forces and energies to avoid segmentation faults.

With the corrected *j* index, it is now possible to compute the distance
of the pair.  For efficiency reasons, the square root is only taken
*after* the check for the cutoff (which has been stored as squared
cutoff by the ``Pair`` base class).  For some pair styles, like the 12-6
Lennard-Jones potential, computing the square root can be avoided
entirely.

.. code-block:: c++

       for (jj = 0; jj < jnum; jj++) {
         j = jlist[jj];
         factor_lj = special_lj[sbmask(j)];
         j &= NEIGHMASK;

         delx = xtmp - x[j][0];
         dely = ytmp - x[j][1];
         delz = ztmp - x[j][2];
         rsq = delx * delx + dely * dely + delz * delz;
         jtype = type[j];

The following block of code is the actual application of the model
potential to compute the force.  Note, that *fpair* is the pair-wise
force divided by the distance, as this simplifies the projection of the
x-, y-, and z-components of the force vector by simply multiplying with
the respective distances in those directions.

.. code-block:: c++

         if (rsq < cutsq[itype][jtype]) {
           r = sqrt(rsq);
           dr = r - r0[itype][jtype];
           aexp = biga0[itype][jtype] * exp(-alpha[itype][jtype] * r);
           bexp = biga1[itype][jtype] * exp(-beta[itype][jtype] * dr * dr);
           fpair = alpha[itype][jtype] * aexp;
           fpair -= 2.0 * beta[itype][jtype] * dr * bexp;
           fpair *= factor_lj / r;

In the next block, the force is added to the per-atom force arrays.  This
pair style uses a "half" neighbor list (each pair is listed only once)
so we take advantage of the fact that :math:`\vec{F}_{ij} =
-\vec{F}_{ji}`, i.e.  apply Newton's third law.  The force is *always*
stored when the atom is a "local" atom. Index *i* atoms are always "local"
(i.e. *i* < nlocal); index *j* atoms may be "ghost" atoms (*j* >= nlocal).

Depending on the settings used with the :doc:`newton command <newton>`,
those pairs are only listed once globally (newton_pair == 1), then
forces must be stored even with ghost atoms and after all forces are
computed a "reverse communication" is performed to add those ghost atom
forces to their corresponding local atoms.  If the setting is disabled,
then the extra communication is skipped, since for pairs straddling
sub-domain boundaries, the forces are computed twice and only stored
with the local atoms in the domain that *owns* it.

.. code-block:: c++

           f[i][0] += delx * fpair;
           f[i][1] += dely * fpair;
           f[i][2] += delz * fpair;
           if (newton_pair || j < nlocal) {
             f[j][0] -= delx * fpair;
             f[j][1] -= dely * fpair;
             f[j][2] -= delz * fpair;
           }

The ``ev_tally()`` function tallies global or per-atom energy and
virial.  For typical MD simulations, the potential energy is merely a
diagnostic and only needed on output.  Similarly, the pressure may only
be computed for (infrequent) thermodynamic output.  For all timesteps
where this information is not needed either, `eflag` or `evflag` are
zero and the computation and call to the tally function skipped.  Note
that evdwl is initialized to zero at the beginning of the function, so
that it still is valid to access it, even if the energy is not computed
(e.g. when only the virial is needed).

.. code-block:: c++

           if (eflag) evdwl = factor_lj * (aexp - bexp - offset[itype][jtype]);
           if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
         }
       }
     }

If only the global virial is needed and no energy, then calls to
``ev_tally()`` can be avoided altogether, and the global virial can be
computed more efficiently from the dot product of the total per-atom
force vector and the position vector of the corresponding atom,
:math:`\vec{F}\cdot\vec{r}`.  This has to be done *after* all pair-wise
forces are computed and *before* the reverse communication to collect
data from ghost atoms, since the position has to be the position that was
used to compute the force, i.e. *not* the "local" position if that ghost
atom is a periodic copy.

.. code-block:: c++

     if (vflag_fdotr) virial_fdotr_compute();
   }


Computing force and energy for a single pair
""""""""""""""""""""""""""""""""""""""""""""

Certain features in LAMMPS only require computing interactions between
individual pairs of atoms and the (optional) ``single()`` function is
needed to support those features (e.g. for tabulation of force and
energy with :doc:`pair_write <pair_write>`).  This is a repetition of
the force kernel in the ``compute()`` function, but only for a single
pair of atoms, where the (squared) distance is provided as a parameter
(so it may not even be an existing distance between two specific atoms).
The energy is returned as the return value of the function and the force
as the `fforce` reference.  Note, that this is, similar to how *fpair*
is used in the ``compute()`` function, the magnitude of the force along
the vector between the two atoms *divided* by the distance.

The ``single()`` function is optional, but it is expected to be
implemented for any true pair-wise additive potential. Many-body
potentials and special case potentials do not implement it. In a few
special cases (EAM, long-range Coulomb), the ``single()`` function
implements the pairwise additive part of the complete force interaction
and depends on either pre-computed properties (derivative of embedding
term for EAM) or post-computed non-pair-wise force contributions (KSpace
style in case of long-range Coulomb).

The member variable `single_enable` should be set to 0 in the
constructor, if it is not implemented (its default value is 1).

.. code-block:: c++

   /* ---------------------------------------------------------------------- */

   double PairBornGauss::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                                double /*factor_coul*/, double factor_lj, double &fforce)
   {
     double r, dr, aexp, bexp;

     r = sqrt(rsq);
     dr = r - r0[itype][jtype];
     aexp = biga0[itype][jtype] * exp(-alpha[itype][jtype] * r);
     bexp = biga1[itype][jtype] * exp(-beta[itype][jtype] * dr * dr);

     fforce = factor_lj * (alpha[itype][jtype] * aexp - 2.0 * dr * beta[itype][jtype] * bexp) / r;
     return factor_lj * (aexp - bexp - offset[itype][jtype]);
   }


Reading and writing of restart files
""""""""""""""""""""""""""""""""""""

Support for writing and reading binary restart files is provided by the
following four functions.  Writing is only done by MPI processor rank 0.
The output of global (not related to atom types) settings is usually
delegated to the ``write_restart_settings()`` function.  This restart
facility is commonly only used, if there are small number of per-type
parameters.  For potentials that use per-element parameters or tabulated
data and read these from files, those parameters and the name of the
potential file are not written to restart files and the :doc:`pair_coeff
command <pair_coeff>` has to re-issued when restarting.  For pair styles
like "born/gauss" that do support writing to restart files, this is not
required.

Implementing the functions to read and write binary restart files is
optional.  The member variable `restartinfo` should be set to 0 in the
constructor, if they are not implemented (its default value is 1).

.. code-block:: c++

   /* ----------------------------------------------------------------------
      proc 0 writes to restart file
   ------------------------------------------------------------------------- */

   void PairBornGauss::write_restart(FILE *fp)
   {
     write_restart_settings(fp);

     int i, j;
     for (i = 1; i <= atom->ntypes; i++) {
       for (j = i; j <= atom->ntypes; j++) {
         fwrite(&setflag[i][j], sizeof(int), 1, fp);
         if (setflag[i][j]) {
           fwrite(&biga0[i][j], sizeof(double), 1, fp);
           fwrite(&alpha[i][j], sizeof(double), 1, fp);
           fwrite(&biga1[i][j], sizeof(double), 1, fp);
           fwrite(&beta[i][j], sizeof(double), 1, fp);
           fwrite(&r0[i][j], sizeof(double), 1, fp);
           fwrite(&cut[i][j], sizeof(double), 1, fp);
         }
       }
     }
   }

   /* ----------------------------------------------------------------------
      proc 0 writes to restart file
   ------------------------------------------------------------------------- */

   void PairBornGauss::write_restart_settings(FILE *fp)
   {
     fwrite(&cut_global, sizeof(double), 1, fp);
     fwrite(&offset_flag, sizeof(int), 1, fp);
     fwrite(&mix_flag, sizeof(int), 1, fp);
   }

Similarly, on reading, only MPI processor rank 0 has opened the restart
file and will read the data.  The data is then distributed across all
parallel processes using calls to ``MPI_Bcast()``.  Before reading atom
type specific data, the corresponding storage needs to be allocated.
Order and number or storage size of items read must be exactly the same
as when writing, or else the data will be read incorrectly.

Reading uses the :cpp:func:`utils::sfread <LAMMPS_NS::utils::sfread>`
utility function to detect read errors and short reads, so that LAMMPS
can abort if that happens, e.g. when the restart file is corrupted.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */

   void PairBornGauss::read_restart(FILE *fp)
   {
     read_restart_settings(fp);

     allocate();

     int i, j;
     int me = comm->me;
     for (i = 1; i <= atom->ntypes; i++) {
       for (j = i; j <= atom->ntypes; j++) {
         if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
         MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
         if (setflag[i][j]) {
           if (me == 0) {
             utils::sfread(FLERR, &biga0[i][j], sizeof(double), 1, fp, nullptr, error);
             utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
             utils::sfread(FLERR, &biga1[i][j], sizeof(double), 1, fp, nullptr, error);
             utils::sfread(FLERR, &beta[i][j], sizeof(double), 1, fp, nullptr, error);
             utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
             utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
           }
           MPI_Bcast(&biga0[i][j], 1, MPI_DOUBLE, 0, world);
           MPI_Bcast(&alpha[i][j], 1, MPI_DOUBLE, 0, world);
           MPI_Bcast(&biga1[i][j], 1, MPI_DOUBLE, 0, world);
           MPI_Bcast(&beta[i][j], 1, MPI_DOUBLE, 0, world);
           MPI_Bcast(&r0[i][j], 1, MPI_DOUBLE, 0, world);
           MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
         }
       }
     }
   }

   /* ----------------------------------------------------------------------
      proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */

   void PairBornGauss::read_restart_settings(FILE *fp)
   {
     if (comm->me == 0) {
       utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
       utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
       utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
     }
     MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
     MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
     MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
   }

Writing coefficients to data files
""""""""""""""""""""""""""""""""""

The ``write_data()`` and ``write_data_all()`` functions are optional and
write out the current state of the :doc:`pair_coeff
settings<pair_coeff>` as "Pair Coeffs" or "PairIJ Coeffs" sections to a
data file when using the :doc:`write_data command <write_data>`.  The
``write_data()`` only writes out the diagonal elements of the pair
coefficient matrix, as that is required for the format of the "Pair
Coeffs" section.  It is called when the "pair" option of the
:doc:`write_data command <write_data>` is "ii" (the default).  This is
suitable for force fields where *all* off-diagonal terms of the pair
coefficient matrix are generated from mixing.  If explicit settings for
off-diagonal elements were made, LAMMPS will print a warning, as those
would be lost.  To avoid this, the "pair ij" option of :doc:`write_data
<write_data>` can be used which will trigger calling the
``write_data_all()`` function instead, which will write out all settings
of the pair coefficient matrix (regardless of whether they were
originally created from mixing or not).

These data file output functions are only useful for true pair-wise
additive potentials, where the potential parameters can be entered
through *multiple* :doc:`pair_coeff commands <pair_coeff>`.  Pair styles
that require a single "pair_coeff \* \*" command line are not compatible
with reading their parameters from data files.  For pair styles like
*born/gauss* that do support writing to data files, the potential
parameters will be read from the data file, if present and
:doc:`pair_coeff commands <pair_coeff>` may not be needed.

The member variable ``writedata`` should be set to 1 in the constructor,
if these functions are implemented (the default value is 0).

.. code-block:: c++

   /* ----------------------------------------------------------------------
      proc 0 writes to data file
   ------------------------------------------------------------------------- */

   void PairBornGauss::write_data(FILE *fp)
   {
     for (int i = 1; i <= atom->ntypes; i++)
       fprintf(fp, "%d %g %g %g %g %g\n", i, biga0[i][i], alpha[i][i], biga1[i][i], beta[i][i],
               r0[i][i]);
   }

   /* ----------------------------------------------------------------------
      proc 0 writes all pairs to data file
   ------------------------------------------------------------------------- */

   void PairBornGauss::write_data_all(FILE *fp)
   {
     for (int i = 1; i <= atom->ntypes; i++)
       for (int j = i; j <= atom->ntypes; j++)
         fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, biga0[i][j], alpha[i][j], biga1[i][j],
                 beta[i][j], r0[i][j], cut[i][j]);
   }


Give access to internal data
""""""""""""""""""""""""""""

The purpose of the ``extract()`` function is to facilitate access to
internal data of the pair style by other parts of LAMMPS.  One possible
application is to use :doc:`fix adapt <fix_adapt>` to gradually change
potential parameters during a run.  Here, we implement access to the
pair coefficient matrix parameters.

.. code-block:: c++

   /* ---------------------------------------------------------------------- */

   void *PairBornGauss::extract(const char *str, int &dim)
   {
     dim = 2;
     if (strcmp(str, "biga0") == 0) return (void *) biga0;
     if (strcmp(str, "biga1") == 0) return (void *) biga1;
     if (strcmp(str, "r0") == 0) return (void *) r0;
     return nullptr;
   }

Since the mercury potential, for which we have implemented the
born/gauss pair style, has a temperature dependent parameter "biga1", we
can automatically adapt the potential based on the Taylor-MacLaurin
expansion for "biga1" when performing a simulation with a temperature
ramp.  LAMMPS commands for that application are given below:

.. code-block:: LAMMPS

   variable tlo  index 300.0
   variable thi  index 600.0
   variable temp equal ramp(v_tlo,v_thi)
   variable biga1 equal (-2.58717e-8*v_temp+8.40841e-5)*v_temp+1.97475e-2

   fix             1 all nvt temp ${tlo} ${thi} 0.1
   fix             2 all adapt 1 pair born/gauss biga1 * * v_biga1

Case 2: a many-body potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since there is a detailed description of the purpose and general layout
of a pair style in the previous case, we will focus on where the
implementation of a typical many-body potential *differs* from a
pair-wise additive potential.  We will use the implementation of the
Tersoff potential as :doc:`pair_style tersoff <pair_tersoff>` as an
example.  The complete implementation can be found in the files
``src/MANYBODY/pair_tersoff.cpp`` and ``src/MANYBODY/pair_tersoff.h`` of
the LAMMPS source code.

Constructor
"""""""""""

In the constructor, several :doc:`pair style flags <Modify_pair>` must
be set differently for many-body potentials:

- the potential is not pair-wise additive, so the ``single()`` function
  cannot be used. This is indicated by setting the `single_enable`
  member variable to 0 (default value is 1)
- many-body potentials are usually not written to :doc:`binary
  restart files <write_restart>`.  This is indicated by setting the member
  variable `restartinfo` to 0 (default is 1)
- many-body potentials typically read *all* parameters from a file which
  stores parameters indexed with a string (e.g. the element).  For this,
  only a single :doc:`pair_coeff \* \* <pair_coeff>` command is allowed.
  This requirement is set and checked for, when the member variable
  `one_coeff` is set to 1 (default value is 0)
- many-body potentials can produce incorrect results if pairs of atoms
  are excluded from the neighbor list, e.g. explicitly by
  :doc:`neigh_modify exclude <neigh_modify>` or implicitly through
  defining bonds, angles, etc. and having a :doc:`special_bonds setting
  <special_bonds>` that is not "special_bonds lj/coul 1.0 1.0 1.0".
  LAMMPS will check for this and print a suitable warning, when the
  member variable `manybody_flag` is set to 1 (default value is 0).

.. code-block:: c++

   PairTersoff::PairTersoff(LAMMPS *lmp) : Pair(lmp)
   {
     single_enable = 0;
     restartinfo = 0;
     one_coeff = 1;
     manybody_flag = 1;

Neighbor list request
"""""""""""""""""""""

For computing the three-body interactions of the Tersoff potential a
"full" neighbor list (both atoms of a pair are listed in each other's
neighbor list) is required.  By default a "half" neighbor list is
requested (each pair is listed only once).  The request is made in
the ``init_style()`` function.  A more in-depth discussion of neighbor
lists in LAMMPS and how to request them is in :ref:`this section of the
documentation <request-neighbor-list>`

Also, additional conditions must be met for some global settings which
are checked in the ``init_style()`` function.

.. code-block:: c++

   /* ----------------------------------------------------------------------
      init specific to this pair style
   ------------------------------------------------------------------------- */

   void PairTersoff::init_style()
   {
     if (atom->tag_enable == 0)
       error->all(FLERR,"Pair style Tersoff requires atom IDs");
     if (force->newton_pair == 0)
       error->all(FLERR,"Pair style Tersoff requires newton pair on");

     // need a full neighbor list

     neighbor->add_request(this,NeighConst::REQ_FULL);
   }

Computing forces from the neighbor list
"""""""""""""""""""""""""""""""""""""""

Computing forces for a many-body potential is usually more complex than
for a pair-wise additive potential and there are multiple components.
For Tersoff, there is a pair-wise additive two-body term (two nested
loops over indices *i* and *j*) and a three-body term (three nested
loops over indices *i*, *j*, and *k*).  Since the neighbor list has
all neighbors up to the maximum cutoff (for the two-body term), but
the three-body interactions have a significantly shorter cutoff,
a "short neighbor list" is also constructed at the same time while computing
the two-body term and looping over the neighbor list for the first time.

.. code-block:: c++

   if (rsq < cutshortsq) {
     neighshort[numshort++] = j;
     if (numshort >= maxshort) {
       maxshort += maxshort/2;
       memory->grow(neighshort,maxshort,"pair:neighshort");
     }
   }

For the two-body term, only a half neighbor list would be needed, even
though we have requested a full list (for the three-body loops).
Rather than computing all interactions twice, we skip over half of
the entries.  This is done in a slightly complex way to make certain
the same choice is made across all subdomains and so that there is
no load imbalance introduced.

.. code-block:: c++

   jtag = tag[j];
   if (itag > jtag) {
     if ((itag+jtag) % 2 == 0) continue;
   } else if (itag < jtag) {
     if ((itag+jtag) % 2 == 1) continue;
   } else {
     if (x[j][2] < x[i][2]) continue;
     if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
     if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
   }

For the three-body term, there is one additional nested loop and it uses
the "short" neighbor list, accumulated previously.

.. code-block:: c++

   // three-body interactions
   // skip immediately if I-J is not within cutoff
   double fjxtmp,fjytmp,fjztmp;

   for (jj = 0; jj < numshort; jj++) {
     j = neighshort[jj];
     jtype = map[type[j]];

     [...]

     for (kk = 0; kk < numshort; kk++) {
       if (jj == kk) continue;
       k = neighshort[kk];
       ktype = map[type[k]];

       [...]
     }
   [...]


Reading potential parameters
""""""""""""""""""""""""""""

For the Tersoff potential, the parameters are listed in a file and
associated with triples of elements.  Because we have set the
``one_coeff`` flag to 1 in the constructor, there may only be a single
:doc:`pair_coeff \* \* <pair_coeff>` line in the input for this pair
style, and as a consequence the ``coeff()`` function will only be called
once.  Thus, the ``coeff()`` function has to do three tasks, each of
which is delegated to a function in the ``PairTersoff`` class:

#. map elements to atom types.  Those follow the potential file name in the
   command line arguments and are processed by the ``map_element2type()`` function.
#. read and parse the potential parameter file in the ``read_file()`` function.
#. Build data structures where the original and derived parameters are
   indexed by all possible triples of atom types and thus can be looked
   up quickly in the loops for the force computation

.. code-block:: c++

   void PairTersoff::coeff(int narg, char **arg)
   {
     if (!allocated) allocate();

     map_element2type(narg-3,arg+3);

     // read potential file and initialize potential parameters

     read_file(arg[2]);
     setup_params();
   }


Case 3: a potential requiring communication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For some models, the interactions between atoms depends on properties of
their environment which have to be computed *before* the the forces can
be computed.  Since LAMMPS is designed to run in parallel using a
:doc:`domain decomposition strategy <Developer_par_part>`, not all
information of the atoms may be directly available and thus
communication steps may be need to collect data from ghost atoms of
neighboring subdomains or send data to ghost atoms for application
during the pairwise computation.

Specifically, two communication patterns are needed: a "reverse
communication" and a "forward communication".  The reverse communication
collects data added to "ghost" atoms from neighboring sub-domains and
sums it to their corresponding "local" atoms.  This communication is
only required and thus executed when the ``Force::newton_pair`` setting
is 1 (i.e. :doc:`newton on <newton>`, the default).  The forward
communication is used to copy computed per-atom data from "local" atoms
to their corresponding "ghost" atoms in neighboring sub-domains.

For this we will look at how the embedding term of the :doc:`embedded
atom potential EAM <pair_eam>` is implemented in LAMMPS.  The complete
implementation of this pair style can be found in the files
``src/MANYBODY/pair_eam.cpp`` and ``src/MANYBODY/pair_eam.h`` of the
LAMMPS source code.

Allocating additional per-atom storage
""""""""""""""""""""""""""""""""""""""

First suitable (local) per-atom arrays (`rho`, `fp`, `numforce`) are
allocated. These have to be large enough to include ghost atoms, are not
used outside the ``compute()`` function and are re-initialized to zero
once per timestep.

.. code-block:: c++

   if (atom->nmax > nmax) {
     memory->destroy(rho);
     memory->destroy(fp);
     memory->destroy(numforce);
     nmax = atom->nmax;
     memory->create(rho,nmax,"pair:rho");
     memory->create(fp,nmax,"pair:fp");
     memory->create(numforce,nmax,"pair:numforce");
   }

Reverse communication
"""""""""""""""""""""

Then a first loop over all pairs (*i* and *j*) is performed, where data
is stored in the `rho` array representing the electron density at the site of
*i* contributed from all neighbors *j*.  Since the EAM pair style uses
a half neighbor list (for efficiency reasons), a reverse communication is
needed to collect the contributions to `rho` from ghost atoms (only if
:doc:`newton on <newton>` is set for pair styles).

.. code-block:: c++

   if (newton_pair) comm->reverse_comm(this);

To support the reverse communication, two functions must be defined:
``pack_reverse_comm()`` that copies relevant data into a buffer for ghost
atoms and ``unpack_reverse_comm()`` that takes the collected data and adds
it to the `rho` array for the corresponding local atoms that match the
ghost atoms.  In order to allocate sufficiently sized buffers, a flag
must be set in the pair style constructor. Since in this case a single
double precision number is communicated per atom, the `comm_reverse`
member variable is set to 1 (default is 0 = no reverse communication).

.. code-block:: c++

   int PairEAM::pack_reverse_comm(int n, int first, double *buf)
   {
     int i,m,last;

     m = 0;
     last = first + n;
     for (i = first; i < last; i++) buf[m++] = rho[i];
     return m;
   }

   void PairEAM::unpack_reverse_comm(int n, int *list, double *buf)
   {
     int i,j,m;

     m = 0;
     for (i = 0; i < n; i++) {
       j = list[i];
       rho[j] += buf[m++];
     }
   }

Forward communication
"""""""""""""""""""""

From the density array `rho`, the derivative of the embedding energy
`fp` is computed. The computation is only done for "local" atoms, but
for the force computation, that property also is needed on ghost atoms.
For that a forward communication is needed.

.. code-block:: c++

   comm->forward_comm(this);

Similar to the reverse communication, this requires implementing a
``pack_forward_comm()`` and an ``unpack_forward_comm()`` function.
Since there is one double precision number per atom that needs to be
communicated, we must set the `comm_forward` member variable to 1
(default is 0 = no forward communication).

.. code-block:: c++

   int PairEAM::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
   {
     int i,j,m;

     m = 0;
     for (i = 0; i < n; i++) {
       j = list[i];
       buf[m++] = fp[j];
     }
     return m;
   }

   void PairEAM::unpack_forward_comm(int n, int first, double *buf)
   {
     int i,m,last;

     m = 0;
     last = first + n;
     for (i = first; i < last; i++) fp[i] = buf[m++];
   }

Case 4: potentials without a compute() function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A small number of pair style classes do not implement a ``compute()``
function, but instead use that of a different pair style.

Embedded atom variants "eam/fs" and "eam/alloy"
"""""""""""""""""""""""""""""""""""""""""""""""

The pair styles :doc:`eam/fs and eam/alloy <pair_eam>` share the same
model and potential function as the :doc:`eam pair style <pair_eam>`.
They differ in the format of the potential files.  Pair style :doc:`eam
<pair_eam>` supports only potential files for single elements.  For
multi-element systems, the mixed terms are computed from mixed
parameters.  The *eam/fs* and *eam/alloy* pair styles, however,
**require** the use of a single potential file for all elements where
the mixed element potential is included in the tabulation.  That enables
more accurate models for alloys, since the mixed terms can be adjusted
for a better representation of material properties compared to terms
created from mixing of per-element terms in the ``PairEAM`` class.

We take a closer at the *eam/alloy* pair style.  The complete
implementation is in the files ``src/MANYBODY/pair_eam_alloy.cpp`` and
``src/MANYBODY/pair_eam_alloy.h``.

The ``PairEAMAlloy`` class is derived from ``PairEAM`` and not ``Pair``
and overrides only a small number of functions:

.. code-block:: c++

   class PairEAMAlloy : virtual public PairEAM {
    public:
     PairEAMAlloy(class LAMMPS *);
     void coeff(int, char **) override;

    protected:
     void read_file(char *) override;
     void file2array() override;
   };

All other functionality is inherited from the base classes.  In the
constructor we set the ``one_coeff`` flag and the ``many_body`` flag to
1 to indicate the different behavior.

.. code-block:: c++

   PairEAMAlloy::PairEAMAlloy(LAMMPS *lmp) : PairEAM(lmp)
   {
     one_coeff = 1;
     manybody_flag = 1;
   }

The ``coeff()`` function (not shown here) implements the different
behavior when processing the :doc:`pair_coeff command <pair_coeff>`.
The ``read_file()`` and ``file2array()`` replace the corresponding
``PairEAM`` class functions to accommodate the different data and
file format.

AIREBO and AIREBO-M potentials
""""""""""""""""""""""""""""""

The AIREBO-M potential differs from the better known AIREBO potential in
that it use a Morse potential instead of a Lennard-Jones potential for
non-bonded interactions.  Since this difference is very minimal compared
to the entire potential, both potentials are implemented in the
``PairAIREBO`` class and which non-bonded potential is used is
determined by the value of the ``morseflag`` flag, which would be set to
either 0 or 1.

.. code-block:: c++

   class PairAIREBOMorse : public PairAIREBO {
    public:
     PairAIREBOMorse(class LAMMPS *);
     void settings(int, char **) override;
   };

The ``morseflag`` variable defaults to 0 and is set to 1 in the
``PairAIREBOMorse::settings()`` function which is called by the
:doc:`pair_style <pair_style>` command.  This function delegates
all command line processing and setting of other parameters to the
``PairAIREBO::settings()`` function of the base class.

.. code-block:: c++

   void PairAIREBOMorse::settings(int narg, char **arg)
   {
     PairAIREBO::settings(narg, arg);

     morseflag = 1;
   }

The complete implementation is in the files
``src/MANYBODY/pair_airebo.cpp``, ``src/MANYBODY/pair_airebo.h``,
``src/MANYBODY/pair_airebo_morse.cpp``,
``src/MANYBODY/pair_airebo_morse.h``.

--------------

.. _Bomont2:

**(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)
