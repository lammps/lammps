.. index:: pair_style born/gauss

pair_style born/gauss command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style born/gauss cutoff

* born/gauss = name of the pair style
* cutoff = global cutoff (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style born/gauss 10.0
   pair_coeff 1 1 1 1 8.2464e13 12.48 0.042644277 0.44 3.56

Description
"""""""""""

.. versionadded:: 28Mar2023

Pair style *born/gauss* computes pairwise interactions from a combination of a Born-Mayer
repulsive term and a Gaussian attractive term according to :ref:`(Bomont) <Bomont>`:

.. math::

   E = A_0 \exp \left( -\alpha r \right) - A_1 \exp\left[ -\beta \left(r - r_0 \right)^2 \right]
       \qquad r < r_c

:math:`r_c` is the cutoff.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`A_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`A_1` (energy units)
* :math:`\beta` (1/(distance units)^2)
* :math:`r_0` (distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global cutoff is used.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Thus, coefficients for all I,J
pairs must be specified explicitly.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options are not relevant for
this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is only enabled if LAMMPS was built with the EXTRA-PAIR
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style born <pair_born>`

Default
"""""""

none

--------------

.. _Bomont:

**(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)
