.. index:: pair_style sph/rhosum

pair_style sph/rhosum command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style sph/rhosum Nstep

* Nstep = timestep interval

Examples
""""""""

.. code-block:: LAMMPS

   pair_style sph/rhosum 10
   pair_coeff * * 2.4

Description
"""""""""""

The sph/rhosum style computes the local particle mass density rho for
SPH particles by kernel function interpolation, every Nstep timesteps.

See `this PDF guide <PDF/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above.

* h (distance units)

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.  Thus, you need to re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command.  It does not support the *inner*,
*middle*, *outer* keywords.

Restrictions
""""""""""""

This pair style is part of the SPH package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, pair_sph/taitwater

Default
"""""""

none
