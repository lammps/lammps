.. index:: pair_style harmonic/cut
.. index:: pair_style harmonic/cut/omp

pair_style harmonic/cut command
===============================

Accelerator Variants: *harmonic/cut/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *harmonic/cut*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style harmonic/cut
   pair_coeff * * 0.2 2.0
   pair_coeff 1 1 0.5 2.5

Description
"""""""""""

Style *harmonic/cut* computes pairwise repulsive-only harmonic interactions with the formula

.. math::

   E = k (r_c - r)^2  \qquad r < r_c

where :math:`r_c` is the cutoff.  Note that the usual 1/2 factor is
included in :math:`k`.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`k` (energy/distance^2 units)
* :math:`r_c` (distance units)

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the :math:`k` and :math:`r_c`
coefficients can be mixed. The default mix value is *geometric*.
See the "pair_modify" command for details.

Since the potential is zero at and beyond the cutoff parameter by
construction, there is no need to support the :doc:`pair_modify
<pair_modify>` shift or tail options for the energy and pressure of the
pair interaction.

These pair styles write their information to :doc:`binary restart files <restart>`,
so pair_style and pair_coeff commands do not need to be specified in an input script
that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *harmonic/cut* pair style is only enabled if LAMMPS was
built with the EXTRA-PAIR package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
