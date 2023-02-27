Writing new pair styles
^^^^^^^^^^^^^^^^^^^^^^^

Pair styles are at the core of most simulations with LAMMPS since they
are used to compute the forces (plus energy and virial contributions, if
needed) on atoms for pairs of atoms within a given cutoff.  This is
often the dominant computation in LAMMPS and sometimes even the only
one.  Pair styles can be grouped in multiple categories:

#. simple pairwise additive interactions of point particles
   (e.g. :doc:`Lennard-Jones <pair_lj>`, :doc:`Morse <pair_morse>`,
   :doc:`Buckingham <pair_buck>`)
#. pairwise additive interactions of point particles with added
   :doc:`Coulomb <pair_coul>` interactions
#. manybody interactions of point particles (e.g. :doc:`EAM <pair_eam>`,
   :doc:`Tersoff <pair_tersoff>`)
#. complex interactions that include additional per-atom properties
   (e.g. Discrete Element Models (DEM), Peridynamics, Ellipsoids,
   Smoothed Particle Hydrodynamics (SPH))

In the text below we will discuss aspects of implementing pair styles in
LAMMPS by looking at representative case studies.  The design of LAMMPS
allows to focus on the essentials, which is to compute the forces (and
energies or virial contributions), enter and manage the global settings
as well as the potential parameters, and the pair style specific parts
of reading and writing restart and data files.  Most of the complex
tasks like management of the atom positions, domain decomposition and
boundaries, or neighbor list creation are handled transparently by other
parts of the LAMMPS code.

As shown on the page for :doc:`writing or extending pair styles
<Modify_pair>` for the implementation of a pair style a new class must
be written that is either directly or indirectly derived from the
``Pair`` class.  In that derived class there are three *required*
methods in addition to the constructor that must be implemented since
they are "pure" in the base class: ``Pair::compute()``,
``Pair::settings()``, ``Pair::coeff()``.  All other methods are optional
and have default implementations in the base class (most of which do
nothing), but they may need to be overridden depending on the
requirements of the model.

Case 1: a pairwise additive model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we will describe the complete procedure of adding
a simple pair style to LAMMPS: an empirical model that can be used
to model liquid mercury.

Model and general considerations
""""""""""""""""""""""""""""""""

The functional form of the model according to :ref:`(Bomont) <Bomont2>`
consists of a repulsive Born-Mayer exponential term and a temperature
dependent, attractive Gaussian term.

.. math::

   E = A_0 \exp \left( -\alpha r \right) - A_1 \exp\left[ -\beta \left(r - r_0 \right)^2 \right]

For the application to mercury the following parameters are listed:

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
parameters for each pair of atom types.  In addition we need to get the
temperature and a default cutoff value as global settings.

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
extension suggests a C source.

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

Every pair style must be registered in LAMMPS by writing the following
lines of code in the second part of the header before the include guards
for the class definition:

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
class, ``PairBornGauss``.  The list of header files to include is
automatically updated by the build system, so the presence of the file
in the ``src/EXTRA-PAIR`` folder and the enabling of the EXTRA-PAIR
package will trigger that LAMMPS includes the new pair style when it is
(re-)compiled.  The "// clang-format" format comments are needed so that
running :ref:`clang-format <clang-format>` on the file will not insert
blanks between "born", "/", and "gauss" which would break the
``PairStyle`` macro.

The third segment of the header is the actual class definition of the
``PairBornGauss`` class.  This has the prototypes for all member
functions that will be implemented by this pair style.  This includes a
number of optional functions.  All functions that were labeled in the
base class as "virtual" must be given the "override" property as it is
done in the code shown below.  This helps to detect unexpected
mismatches as compile errors in case the signature of a function is
changed in the base class.  For example, if this change would add an
optional argument with a default value, then all existing source code
calling the function would not need changes and still compile, but the
function in the derived class would no longer override the one in the
base class due to the different number of arguments and the behavior of
the pair style is thus changed in an unintended way. Using "override"
prevents such issues.

Also variables and arrays for storing global settings and potential
parameters are defined.  Since those are internal to the class, they are
placed after a "protected:" label.

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

    protected:
     double cut_global, temperature;
     double **cut;
     double **biga0, **alpha, **biga1, **beta, **r0;
     double **a0, **a1, **a2;
     double **offset;

     virtual void allocate();
   };
   }    // namespace LAMMPS_NS
   #endif

Some details of the class definition will be discussed later.

Implementation file
"""""""""""""""""""

We move on to the implementation in the ``pair_born_gauss.cpp`` file.
This file also starts with the LAMMPS copyright and license header.
Below is typically the space where comments may be added with additional
information about this specific file, the author(s) and affiliation(s)
and email address(es) so the author can be easily contacted in case
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
then LAMMPS classes (sorted alphabetically) and system headers and
others, if needed.  Note the standardized C++ notation for headers of
C-library functions. The final statement of this segment imports the
``LAMMPS_NS::`` namespace globally for this file.  This way, all LAMMPS
specific functions and classes do not have to be prefixed with
``LAMMPS_NS::``.

Constructor and destructor (required)
"""""""""""""""""""""""""""""""""""""

The first two functions in the implementation source file are typically
the constructor and the destructor.

.. code-block:: c++

   /* ---------------------------------------------------------------------- */

   PairBornGauss::PairBornGauss(LAMMPS *lmp) : Pair(lmp), temp(nullptr)
   {
     single_enable = 1;
     respa_enable = 0;
     writedata = 1;
   }

   /* ---------------------------------------------------------------------- */

   PairBornGauss::~PairBornGauss()
   {
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

Settings and coefficients (required)
""""""""""""""""""""""""""""""""""""

The next method we need to implement is ``setmask()``:

.. code-block:: c++

   int FixPrintVel::setmask()
   {
     int mask = 0;
     mask |= FixConst::END_OF_STEP;
     return mask;
   }

Computing forces from the neighbor list (required)
""""""""""""""""""""""""""""""""""""""""""""""""""

xxxxxadfasdf

--------------

.. _Bomont2:

**(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)
