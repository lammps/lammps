.. index:: suffix

suffix command
==============

Syntax
""""""

.. parsed-literal::

   suffix style args

* style = *off* or *on* or *gpu* or *intel* or *kk* or *omp* or *opt* or *hybrid*
* args = for hybrid style, default suffix to be used and alternative suffix

Examples
""""""""

.. code-block:: LAMMPS

   suffix off
   suffix on
   suffix gpu
   suffix intel
   suffix hybrid intel omp
   suffix kk

Description
"""""""""""

This command allows you to use variants of various styles if they
exist.  In that respect it operates the same as the :doc:`-suffix command-line switch <Run_options>`.  It also has options to turn
off or back on any suffix setting made via the command line.

The specified style can be *gpu*\ , *intel*\ , *kk*\ , *omp*\ , *opt* or
*hybrid*\ . These refer to optional packages that LAMMPS can be built
with, as described on the :doc:`Build package <Build_package>` doc page.
The "gpu" style corresponds to the GPU package, the "intel" style to
the USER-INTEL package, the "kk" style to the KOKKOS package, the
"omp" style to the USER-OMP package, and the "opt" style to the OPT
package.

These are the variants these packages provide:

* GPU = a handful of pair styles and the PPPM kspace_style, optimized to
  run on one or more GPUs or multicore CPU/GPU nodes
* USER-INTEL = a collection of pair styles and neighbor routines
  optimized to run in single, mixed, or double precision on CPUs and
  Intel(R) Xeon Phi(TM) co-processors.
* KOKKOS = a collection of atom, pair, and fix styles optimized to run
  using the Kokkos library on various kinds of hardware, including GPUs
  via CUDA and many-core chips via OpenMP or threading.
* USER-OMP = a collection of pair, bond, angle, dihedral, improper,
  kspace, compute, and fix styles with support for OpenMP
  multi-threading
* OPT = a handful of pair styles, cache-optimized for faster CPU
  performance
* HYBRID = a combination of two packages can be specified (see below)

As an example, all of the packages provide a :doc:`pair_style lj/cut <pair_lj>` variant, with style names lj/cut/opt, lj/cut/omp,
lj/cut/gpu, lj/cut/intel, or lj/cut/kk.  A variant styles
can be specified explicitly in your input script, e.g. pair_style
lj/cut/gpu. If the suffix command is used with the appropriate style,
you do not need to modify your input script.  The specified suffix
(opt,omp,gpu,intel,kk) is automatically appended whenever your
input script command creates a new :doc:`atom <atom_style>`,
:doc:`pair <pair_style>`, :doc:`bond <bond_style>`,
:doc:`angle <angle_style>`, :doc:`dihedral <dihedral_style>`,
:doc:`improper <improper_style>`, :doc:`kspace <kspace_style>`,
:doc:`fix <fix>`, :doc:`compute <compute>`, or :doc:`run <run_style>` style.
If the variant version does not exist, the standard version is
created.

For "hybrid", two packages are specified. The first is used whenever
available. If a style with the first suffix is not available, the style
with the suffix for the second package will be used if available. For
example, "hybrid intel omp" will use styles from the USER-INTEL package
as a first choice and styles from the USER-OMP package as a second choice
if no USER-INTEL variant is available.

If the specified style is *off*\ , then any previously specified suffix
is temporarily disabled, whether it was specified by a command-line
switch or a previous suffix command.  If the specified style is *on*\ ,
a disabled suffix is turned back on.  The use of these 2 commands lets
your input script use a standard LAMMPS style (i.e. a non-accelerated
variant), which can be useful for testing or benchmarking purposes.
Of course this is also possible by not using any suffix commands, and
explicitly appending or not appending the suffix to the relevant
commands in your input script.

.. note::

   The default :doc:`run_style <run_style>` verlet is invoked prior to
   reading the input script and is therefore not affected by a suffix command
   in the input script. The KOKKOS package requires "run_style verlet/kk",
   so when using the KOKKOS package it is necessary to either use the command
   line "-sf kk" command or add an explicit "run_style verlet" command to the
   input script.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`-suffix command-line switch <Run_options>`

**Default:** none
