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
   pair_coeff * * Ti.meam.spline Ti O

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

As an example, imagine the Ti.meam.spline file has values for Ti (old style).
In that case your LAMMPS simulation may only have one atom type which has
to be mapped to the Ti element as follows:

.. code-block:: LAMMPS

   pair_coeff * * Ti.meam.spline Ti

The first 2 arguments must be \* \* and there may be only one element
following or NULL. Systems where there would be multiple atom types
assigned to the same element are **not** supported by this pair style
due to limitations in its implementation.  If a mapping value is
specified as NULL, the mapping is not performed.  This can be used when
a *meam/spline* potential is used as part of the *hybrid* pair style.
The NULL values are placeholders for atom types that will be used with
other potentials.

An example with a two component spline (new style) is TiO.meam.spline, where
the command

.. code-block:: LAMMPS

   pair_coeff * * TiO.meam.spline Ti O

will map the first atom type to Ti and the second atom type to O. Note
in this case that the species names need to match exactly with the
names of the elements in the TiO.meam.spline file; otherwise an
error will be raised. This behavior is different than the old style
MEAM files.

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

This pair style does not support mapping multiple atom types to the same element.

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

**(Zhang)** Zhang and Trinkle, Computational Materials Science, 124, 204-210 (2016).
