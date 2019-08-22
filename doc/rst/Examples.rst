Example scripts
===============

The LAMMPS distribution includes an examples sub-directory with many
sample problems.  Many are 2d models that run quickly and are
straightforward to visualize, requiring at most a couple of minutes to
run on a desktop machine.  Each problem has an input script (in.\*) and
produces a log file (log.\*) when it runs.  Some use a data file
(data.\*) of initial coordinates as additional input.  A few sample log
file run on different machines and different numbers of processors are
included in the directories to compare your answers to.  E.g. a log
file like log.date.crack.foo.P means the "crack" example was run on P
processors of machine "foo" on that date (i.e. with that version of
LAMMPS).

Many of the input files have commented-out lines for creating dump
files and image files.

If you uncomment the :doc:`dump <dump>` command in the input script, a
text dump file will be produced, which can be animated by various
`visualization programs <http://lammps.sandia.gov/viz.html>`_.

If you uncomment the :doc:`dump image <dump>` command in the input
script, and assuming you have built LAMMPS with a JPG library, JPG
snapshot images will be produced when the simulation runs.  They can
be quickly post-processed into a movie using commands described on the
:doc:`dump image <dump_image>` doc page.

Animations of many of the examples can be viewed on the Movies section
of the `LAMMPS web site <lws_>`_.

There are two kinds of sub-directories in the examples dir.  Lowercase
dirs contain one or a few simple, quick-to-run problems.  Uppercase
dirs contain up to several complex scripts that illustrate a
particular kind of simulation method or model.  Some of these run for
longer times, e.g. to measure a particular quantity.

Lists of both kinds of directories are given below.


----------


Lowercase directories
---------------------

+-------------+------------------------------------------------------------------+
| accelerate  | run with various acceleration options (OpenMP, GPU, Phi)         |
+-------------+------------------------------------------------------------------+
| airebo      | polyethylene with AIREBO potential                               |
+-------------+------------------------------------------------------------------+
| atm         | Axilrod-Teller-Muto potential example                            |
+-------------+------------------------------------------------------------------+
| balance     | dynamic load balancing, 2d system                                |
+-------------+------------------------------------------------------------------+
| body        | body particles, 2d system                                        |
+-------------+------------------------------------------------------------------+
| cmap        | CMAP 5-body contributions to CHARMM force field                  |
+-------------+------------------------------------------------------------------+
| colloid     | big colloid particles in a small particle solvent, 2d system     |
+-------------+------------------------------------------------------------------+
| comb        | models using the COMB potential                                  |
+-------------+------------------------------------------------------------------+
| controller  | use of fix controller as a thermostat                            |
+-------------+------------------------------------------------------------------+
| coreshell   | core/shell model using CORESHELL package                         |
+-------------+------------------------------------------------------------------+
| crack       | crack propagation in a 2d solid                                  |
+-------------+------------------------------------------------------------------+
| deposit     | deposit atoms and molecules on a surface                         |
+-------------+------------------------------------------------------------------+
| dipole      | point dipolar particles, 2d system                               |
+-------------+------------------------------------------------------------------+
| dreiding    | methanol via Dreiding FF                                         |
+-------------+------------------------------------------------------------------+
| eim         | NaCl using the EIM potential                                     |
+-------------+------------------------------------------------------------------+
| ellipse     | ellipsoidal particles in spherical solvent, 2d system            |
+-------------+------------------------------------------------------------------+
| flow        | Couette and Poiseuille flow in a 2d channel                      |
+-------------+------------------------------------------------------------------+
| friction    | frictional contact of spherical asperities between 2d surfaces   |
+-------------+------------------------------------------------------------------+
| gcmc        | Grand Canonical Monte Carlo (GCMC) via the fix gcmc command      |
+-------------+------------------------------------------------------------------+
| granregion  | use of fix wall/region/gran as boundary on granular particles    |
+-------------+------------------------------------------------------------------+
| hugoniostat | Hugoniostat shock dynamics                                       |
+-------------+------------------------------------------------------------------+
| hyper       | global and local hyperdynamics of diffusion on Pt surface        |
+-------------+------------------------------------------------------------------+
| indent      | spherical indenter into a 2d solid                               |
+-------------+------------------------------------------------------------------+
| kim         | use of potentials from the `OpenKIM Repository <openkim_>`_      |
+-------------+------------------------------------------------------------------+
| latte       | examples for using fix latte for DFTB via the LATTE library      |
+-------------+------------------------------------------------------------------+
| meam        | MEAM test for SiC and shear (same as shear examples)             |
+-------------+------------------------------------------------------------------+
| melt        | rapid melt of 3d LJ system                                       |
+-------------+------------------------------------------------------------------+
| message     | demos for LAMMPS client/server coupling with the MESSAGE package |
+-------------+------------------------------------------------------------------+
| micelle     | self-assembly of small lipid-like molecules into 2d bilayers     |
+-------------+------------------------------------------------------------------+
| min         | energy minimization of 2d LJ melt                                |
+-------------+------------------------------------------------------------------+
| mscg        | parameterize a multi-scale coarse-graining (MSCG) model          |
+-------------+------------------------------------------------------------------+
| msst        | MSST shock dynamics                                              |
+-------------+------------------------------------------------------------------+
| nb3b        | use of non-bonded 3-body harmonic pair style                     |
+-------------+------------------------------------------------------------------+
| neb         | nudged elastic band (NEB) calculation for barrier finding        |
+-------------+------------------------------------------------------------------+
| nemd        | non-equilibrium MD of 2d sheared system                          |
+-------------+------------------------------------------------------------------+
| obstacle    | flow around two voids in a 2d channel                            |
+-------------+------------------------------------------------------------------+
| peptide     | dynamics of a small solvated peptide chain (5-mer)               |
+-------------+------------------------------------------------------------------+
| peri        | Peridynamic model of cylinder impacted by indenter               |
+-------------+------------------------------------------------------------------+
| pour        | pouring of granular particles into a 3d box, then chute flow     |
+-------------+------------------------------------------------------------------+
| prd         | parallel replica dynamics of vacancy diffusion in bulk Si        |
+-------------+------------------------------------------------------------------+
| python      | using embedded Python in a LAMMPS input script                   |
+-------------+------------------------------------------------------------------+
| qeq         | use of the QEQ package for charge equilibration                  |
+-------------+------------------------------------------------------------------+
| rdf-adf     | computing radial and angle distribution functions for water      |
+-------------+------------------------------------------------------------------+
| reax        | RDX and TATB models using the ReaxFF                             |
+-------------+------------------------------------------------------------------+
| rigid       | rigid bodies modeled as independent or coupled                   |
+-------------+------------------------------------------------------------------+
| shear       | sideways shear applied to 2d solid, with and without a void      |
+-------------+------------------------------------------------------------------+
| snap        | NVE dynamics for BCC tantalum crystal using SNAP potential       |
+-------------+------------------------------------------------------------------+
| srd         | stochastic rotation dynamics (SRD) particles as solvent          |
+-------------+------------------------------------------------------------------+
| streitz     | use of Streitz/Mintmire potential with charge equilibration      |
+-------------+------------------------------------------------------------------+
| tad         | temperature-accelerated dynamics of vacancy diffusion in bulk Si |
+-------------+------------------------------------------------------------------+
| threebody   | regression test input for a variety of manybody potentials       |
+-------------+------------------------------------------------------------------+
| vashishta   | use of the Vashishta potential                                   |
+-------------+------------------------------------------------------------------+
| voronoi     | Voronoi tesselation via compute voronoi/atom command             |
+-------------+------------------------------------------------------------------+

Here is how you can run and visualize one of the sample problems:


.. parsed-literal::

   cd indent
   cp ../../src/lmp_linux .           # copy LAMMPS executable to this dir
   lmp_linux -in in.indent            # run the problem

Running the simulation produces the files *dump.indent* and
*log.lammps*\ .  You can visualize the dump file of snapshots with a
variety of 3rd-party tools highlighted on the
`Visualization <http://lammps.sandia.gov/viz.html>`_ page of the LAMMPS
web site.

If you uncomment the :doc:`dump image <dump_image>` line(s) in the input
script a series of JPG images will be produced by the run (assuming
you built LAMMPS with JPG support; see the
:doc:`Build\_settings <Build_settings>` doc page for details).  These can
be viewed individually or turned into a movie or animated by tools
like ImageMagick or QuickTime or various Windows-based tools.  See the
:doc:`dump image <dump_image>` doc page for more details.  E.g. this
Imagemagick command would create a GIF file suitable for viewing in a
browser.


.. parsed-literal::

   % convert -loop 1 \*.jpg foo.gif


----------


Uppercase directories
---------------------

+------------+--------------------------------------------------------------------------------------------------+
| ASPHERE    | various aspherical particle models, using ellipsoids, rigid bodies, line/triangle particles, etc |
+------------+--------------------------------------------------------------------------------------------------+
| COUPLE     | examples of how to use LAMMPS as a library                                                       |
+------------+--------------------------------------------------------------------------------------------------+
| DIFFUSE    | compute diffusion coefficients via several methods                                               |
+------------+--------------------------------------------------------------------------------------------------+
| ELASTIC    | compute elastic constants at zero temperature                                                    |
+------------+--------------------------------------------------------------------------------------------------+
| ELASTIC\_T | compute elastic constants at finite temperature                                                  |
+------------+--------------------------------------------------------------------------------------------------+
| HEAT       | compute thermal conductivity for LJ and water via fix ehex                                       |
+------------+--------------------------------------------------------------------------------------------------+
| KAPPA      | compute thermal conductivity via several methods                                                 |
+------------+--------------------------------------------------------------------------------------------------+
| MC         | using LAMMPS in a Monte Carlo mode to relax the energy of a system                               |
+------------+--------------------------------------------------------------------------------------------------+
| SPIN       | examples for features of the SPIN package                                                        |
+------------+--------------------------------------------------------------------------------------------------+
| USER       | examples for USER packages and USER-contributed commands                                         |
+------------+--------------------------------------------------------------------------------------------------+
| VISCOSITY  | compute viscosity via several methods                                                            |
+------------+--------------------------------------------------------------------------------------------------+

Nearly all of these directories have README files which give more
details on how to understand and use their contents.

The USER directory has a large number of sub-directories which
correspond by name to a USER package.  They contain scripts that
illustrate how to use the command(s) provided in that package.  Many
of the sub-directories have their own README files which give further
instructions.  See the :doc:`Packages\_details <Packages_details>` doc
page for more info on specific USER packages.

.. _openkim: https://openkim.org




.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
