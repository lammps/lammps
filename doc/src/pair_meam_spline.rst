.. index:: pair_style meam/spline

pair_style meam/spline command
==============================

pair_style meam/spline/omp command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style meam/spline

Examples
""""""""

.. code:: LAMMPS

   pair_style meam/spline
   pair_coeff * * Ti.meam.spline Ti
   pair_coeff * * Ti.meam.spline Ti Ti Ti

Description
"""""""""""

The *meam/spline* style computes pairwise interactions for metals
using a variant of modified embedded-atom method (MEAM) potentials
:ref:`(Lenosky) <Lenosky1>`.  For a single species ("old-style") MEAM,
the total energy E is given by

.. math::

   E & =\sum_{i<j}\phi(r_{ij})+\sum_{i}U(n_{i}) \\
   n_{i} & =\sum_{j}\rho(r_{ij})+\sum_{\substack{j<k,\\j,k\neq i}}f(r_{ij})f(r_{ik})g[\cos(\theta_{jik})]

where :math:`\rho_i` is the density at atom I, :math:`\theta_{jik}` is
the angle between atoms J, I, and K centered on atom I. The five
functions :math:`\phi, U, \rho, f,` and *g* are represented by cubic splines.

The *meam/spline* style also supports a new style multicomponent
modified embedded-atom method (MEAM) potential :ref:`(Zhang) <Zhang4>`, where
the total energy E is given by

.. math::

   E &= \sum_{i<j}\phi_{ij}(r_{ij})+\sum_{i}U_i(n_{i}) \\
   n_{i} & = \sum_{j\ne i}\rho_j(r_{ij})+\sum_{\substack{j<k,\\j,k\neq i}}f_{j}(r_{ij})f_{k}(r_{ik})g_{jk}[\cos(\theta_{jik})]

where the five functions :math:`\phi, U, \rho, f,` and *g* depend on the
chemistry of the atoms in the interaction.  In particular, if there are
N different chemistries, there are N different *U*\ , :math:`\rho`, and
*f* functions, while there are N(N+1)/2 different :math:`\phi` and *g*
functions.  The new style multicomponent MEAM potential files are
indicated by the second line in the file starts with "meam/spline"
followed by the number of elements and the name of each element.

The cutoffs and the coefficients for these spline functions are listed
in a parameter file which is specified by the
:doc:`pair_coeff <pair_coeff>` command.  Parameter files for different
elements are included in the "potentials" directory of the LAMMPS
distribution and have a ".meam.spline" file suffix.  All of these
files are parameterized in terms of LAMMPS :doc:`metal units <units>`.

Note that unlike for other potentials, cutoffs for spline-based MEAM
potentials are not set in the pair_style or pair_coeff command; they
are specified in the potential files themselves.

Unlike the EAM pair style, which retrieves the atomic mass from the
potential file, the spline-based MEAM potentials do not include mass
information; thus you need to use the :doc:`mass <mass>` command to
specify it.

Only a single pair_coeff command is used with the *meam/spline* style
which specifies a potential file with parameters for all needed
elements.  These are mapped to LAMMPS atom types by specifying N
additional arguments after the filename in the pair_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of spline-based MEAM elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine the Ti.meam.spline file has values for Ti (old style).  If
your LAMMPS simulation has 3 atoms types and they are all to be
treated with this potentials, you would use the following pair_coeff
command:

.. code-block:: LAMMPS

   pair_coeff * * Ti.meam.spline Ti Ti Ti

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The three Ti arguments map LAMMPS atom types 1,2,3 to the Ti element
in the potential file.  If a mapping value is specified as NULL, the
mapping is not performed.  This can be used when a *meam/spline*
potential is used as part of the *hybrid* pair style.  The NULL values
are placeholders for atom types that will be used with other
potentials. The old-style potential maps any non-NULL species named
on the command line to that single type.

An example with a two component spline (new style) is TiO.meam.spline, where
the command

.. code-block:: LAMMPS

   pair_coeff * * TiO.meam.spline Ti O

will map the 1st atom type to Ti and the second atom type to O. Note
in this case that the species names need to match exactly with the
names of the elements in the TiO.meam.spline file; otherwise an
error will be raised. This behavior is different than the old style
MEAM files.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

The *meam/spline* pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in an external
potential parameter file.  Thus, you need to re-specify the pair_style
and pair_coeff commands in an input script that reads a restart file.

The *meam/spline* pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

This pair style is only enabled if LAMMPS was built with the USER-MISC
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style meam/c <pair_meamc>`

**Default:** none

----------

.. _Lenosky1:

**(Lenosky)** Lenosky, Sadigh, Alonso, Bulatov, de la Rubia, Kim, Voter,
Kress, Modelling Simulation Materials Science Engineering, 8, 825
(2000).

.. _Zhang4:

**(Zhang)** Zhang and Trinkle, Computational Materials Science, 124, 204-210 (2016).
