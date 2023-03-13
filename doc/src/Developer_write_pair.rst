Writing new pair styles
^^^^^^^^^^^^^^^^^^^^^^^

Pair styles are at the core of most simulations with LAMMPS, since they
are used to compute the forces (plus energy and virial contributions, if
needed) on atoms for pairs of atoms within a given cutoff.  This is
often the dominant computation in LAMMPS, and sometimes even the only
one.  Pair styles can be grouped in multiple categories:

#. simple pairwise additive interactions of point particles
   (e.g. :doc:`Lennard-Jones <pair_lj>`, :doc:`Morse <pair_morse>`,
   :doc:`Buckingham <pair_buck>`)
#. pairwise additive interactions of point particles with added
   :doc:`Coulomb <pair_coul>` interactions
#. manybody interactions of point particles (e.g. :doc:`EAM <pair_eam>`,
   :doc:`Tersoff <pair_tersoff>`)
#. complex interactions that include additional per-atom properties
   (e.g. Discrete Element Models (DEM), Peridynamics, Ellipsoids)

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
``Pair::coeff()``.  In addition a custom constructor is needed.  All
other methods are optional and have default implementations in the base
class (most of which do nothing), but they may need to be overridden
depending on the requirements of the model.

Case 1: a pairwise additive model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, we will describe the procedure of adding a simple pair
style to LAMMPS: an empirical model that can be used to model liquid
mercury.

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
shall be added to the :ref:`EXTRA-PAIR <PKG-EXTRA-PAIR>` package.  For
the implementation, we will use :doc:`pair style morse <pair_morse>` as
template.

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

This second segment of the header file will be included by the ``Force``
class in ``force.cpp`` to build a map of "factory functions" that will
create an instance of these classes and return a pointer to it.  The map
connects the name of the pair style, "born/gauss", to the name of the
class, ``PairBornGauss``.  Before including the headers, the ``PAIR_CLASS``
define is set and the ``PairStyle(name,class)`` macro defined as needed.

The list of header files to include is automatically updated by the
build system, so the presence of the file in the ``src/EXTRA-PAIR``
folder and the enabling of the EXTRA-PAIR package will trigger that
LAMMPS includes the new pair style when it is (re-)compiled.  The "//
clang-format" format comments are needed so that running
:ref:`clang-format <clang-format>` on the file will not insert blanks
between "born", "/", and "gauss" which would break the ``PairStyle``
macro.

The third segment of the header is the actual class definition of the
``PairBornGauss`` class.  This has the prototypes for all member
functions that will be implemented by this pair style.  This includes
:doc:`a few required and a number of optional functions <Modify_pair>`.
All functions that were labeled in the base class as "virtual" must be
given the "override" property, as it is done in the code shown below.

The "override" property helps to detect unexpected mismatches because
compilation will stop with an error in case the signature of a function
is changed in the base class without also changing it in all derived
classes.  For example, if this change would add an optional argument
with a default value, then all existing source code *calling* the
function would not need changes and still compile, but the function in
the derived class would no longer override the one in the base class due
to the different number of arguments and the behavior of the pair style
is thus changed in an unintended way.  Using "override" keyword prevents
such issues.

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
parameters are defined.  Since those are internal to the class, they are
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
pointer as argument, but **not** the command line arguments of the
:doc:`pair_style command <pair_style>`.  Instead, those arguments are
processed in the ``Pair::settings()`` function (or rather the version in
the derived class).  The constructor is the place where global defaults
are set and specifically flags are set about which optional features of
a pair style are available.

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
``Pair::compute_(inner|middle|outer)`` are not, but those are also the
default settings and already set in the base class.

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
but moving this to a separate function allows to change global settings
like the default cutoff without having to reissue all pair_coeff
commands or re-read the ``Pair Coeffs`` sections from the data file.
In the ``settings()`` function, also the arrays for storing parameters,
to define cutoffs, track with pairs of parameters have been explicitly
set are allocated and, if needed, initialized.  In this case, the memory
allocation and initialization is moved to a function ``allocate()``.

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
data files.  In the case of the ``Pair Coeffs`` section there is only
one atom type per line and thus the first argument is duplicated.  Since
the atom type arguments of the :doc:`pair_coeff command <pair_coeff>`
may be a range (e.g. \*\ 3 for atom types 1, 2, and 3), the
corresponding arguments are passed to the :cpp:func:`utils::bounds()
<LAMMPS_NS::utils::bounds>` function which will then return the low
and high end of the range.

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

The ``init()`` function is called during the :doc:`"init" phase
<Developer_flow>` of a simulation.  This is where potential parameters
are checked for completeness, derived parameters computed (e.g. the
"offset" of the potential energy at the cutoff distance for use with the
:doc:`pair_modify shift yes <pair_modify>` command.  If a pair style
supports generating "mixed" parameters (i.e. where both atoms of a pair
have a different atom type) from a "mixing rule" from the parameters of
the type with itself, this is the place to compute and store those mixed
values.  The *born/gauss* pair style does not, so we only check for
completeness.  Another purpose of the ``init()`` function is to
symmetrize the potential parameter arrays.  The return value is the
cutoff for the given pair of atom types.  This is used by the neighbor
list code to determine the largest cutoff and then build the neighbor
lists accordingly.

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

The ``compute()`` function is the "workhorse" of a pair style.

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

     // loop over neighbors of my atoms

     for (ii = 0; ii < inum; ii++) {
       i = ilist[ii];
       xtmp = x[i][0];
       ytmp = x[i][1];
       ztmp = x[i][2];
       itype = type[i];
       jlist = firstneigh[i];
       jnum = numneigh[i];

       for (jj = 0; jj < jnum; jj++) {
         j = jlist[jj];
         factor_lj = special_lj[sbmask(j)];
         j &= NEIGHMASK;

         delx = xtmp - x[j][0];
         dely = ytmp - x[j][1];
         delz = ztmp - x[j][2];
         rsq = delx * delx + dely * dely + delz * delz;
         jtype = type[j];

         if (rsq < cutsq[itype][jtype]) {
           r = sqrt(rsq);
           dr = r - r0[itype][jtype];
           aexp = biga0[itype][jtype] * exp(-alpha[itype][jtype] * r);
           bexp = biga1[itype][jtype] * exp(-beta[itype][jtype] * dr * dr);
           fpair = alpha[itype][jtype] * aexp;
           fpair -= 2.0 * beta[itype][jtype] * dr * bexp;
           fpair *= factor_lj / r;

           f[i][0] += delx * fpair;
           f[i][1] += dely * fpair;
           f[i][2] += delz * fpair;
           if (newton_pair || j < nlocal) {
             f[j][0] -= delx * fpair;
             f[j][1] -= dely * fpair;
             f[j][2] -= delz * fpair;
           }

           if (eflag) evdwl = factor_lj * (aexp - bexp - offset[itype][jtype]);
           if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
         }
       }
     }

     if (vflag_fdotr) virial_fdotr_compute();
   }


Computing force and energy for a single pair
""""""""""""""""""""""""""""""""""""""""""""

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
      proc 0 writes to restart file
   ------------------------------------------------------------------------- */

   void PairBornGauss::write_restart_settings(FILE *fp)
   {
     fwrite(&cut_global, sizeof(double), 1, fp);
     fwrite(&offset_flag, sizeof(int), 1, fp);
     fwrite(&mix_flag, sizeof(int), 1, fp);
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

.. code-block::

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

.. code-block::

   /* ---------------------------------------------------------------------- */

   void *PairBornGauss::extract(const char *str, int &dim)
   {
     dim = 2;
     if (strcmp(str, "biga0") == 0) return (void *) biga0;
     if (strcmp(str, "biga1") == 0) return (void *) biga1;
     if (strcmp(str, "r0") == 0) return (void *) r0;
     return nullptr;
   }

--------------

.. _Bomont2:

**(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)
