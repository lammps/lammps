.. index:: pair_style meam/spline
.. index:: pair_style meam/spline/omp

pair_style meam/spline command
==============================

Accelerator Variants: *meam/spline/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style meam/spline

Examples
""""""""

.. code-block:: LAMMPS

   pair_style meam/spline
   pair_coeff * * Ti.meam.spline Ti
   pair_coeff * * Ti.meam.spline Ti Ti
   pair_coeff * * TiO.meam.spline O Ti Ti

   pair_style hybrid meam/spline lj/cut 10.0
   pair_coeff * * meam/spline TiO.meam.spline Ti O NULL
   pair_coeff * 3 lj/cut 0.01 3.0

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
modified embedded-atom method (MEAM) potential :ref:`(Zhang5) <Zhang4>`, where
the total energy E is given by

.. math::

   E &= \sum_{i<j}\phi_{ij}(r_{ij})+\sum_{i}U_i(n_{i}) \\
   n_{i} & = \sum_{j\ne i}\rho_j(r_{ij})+\sum_{\substack{j<k,\\j,k\neq i}}f_{j}(r_{ij})f_{k}(r_{ik})g_{jk}[\cos(\theta_{jik})]

where the five functions :math:`\phi, U, \rho, f,` and *g* depend on the
chemistry of the atoms in the interaction.  In particular, if there are
N different chemistries, there are N different *U*, :math:`\rho`, and
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

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the potential file.

As an example, the TiO.meam.spline file (new style) has parameters for
Ti and O.  If your LAMMPS simulation has 3 atom types and you want the
first to be O and the second and third to be Ti, you would use the
following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * TiO.meam.spline O Ti Ti

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The O argument maps LAMMPS atom type 1 to the O element in the potential
file.  The two Ti arguments map LAMMPS atom types 2 and 3 to the Ti
element in the potential file.  The element names need to match the
names of the elements in the potential file exactly; otherwise an error
is raised.  Elements of the potential file that are not mapped to any
atom type are not used.

If a mapping value is specified as NULL, the mapping is not performed.
This can be used when a *meam/spline* potential is used as part of the
*hybrid* pair style.  The NULL values are placeholders for atom types
that will be used with other potentials.

Old style potential files, e.g. Ti.meam.spline, have parameters for a
single element and do not list its name.  Any element name can be used
for them in the pair_coeff command, but all atom types that are not
mapped to NULL must use the same name:

.. code-block:: LAMMPS

   pair_coeff * * Ti.meam.spline Ti Ti

.. versionchanged:: 30Sep2026-backport

Previous versions of LAMMPS required each element of the potential file
to be mapped to exactly one atom type, did not support NULL, and
assigned the elements to the atom types in the order of the potential
file, regardless of the order of the element names in the pair_coeff
command.  Also, they used the wrong pair functions for some pairs of
elements with potential files for more than two elements.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

The *meam/spline* pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in an external
potential parameter file.  Thus, you need to re-specify the pair_style
and pair_coeff commands in an input script that reads a restart file.

The *meam/spline* pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

This pair style is only enabled if LAMMPS was built with the MANYBODY
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style meam <pair_meam>`

Default
"""""""

none

----------

.. _Lenosky1:

**(Lenosky)** Lenosky, Sadigh, Alonso, Bulatov, de la Rubia, Kim, Voter,
Kress, Modelling Simulation Materials Science Engineering, 8, 825
(2000).

.. _Zhang4:

**(Zhang5)** Zhang and Trinkle, Computational Materials Science, 124, 204-210 (2016).
