.. index:: pair\_style sph/heatconduction

pair\_style sph/heatconduction command
======================================

Syntax
""""""


.. parsed-literal::

   pair_style sph/heatconduction

Examples
""""""""


.. parsed-literal::

   pair_style sph/heatconduction
   pair_coeff \* \* 1.0 2.4

Description
"""""""""""

The sph/heatconduction style computes heat transport between SPH particles.
The transport model is the diffusion equation for the internal energy.

See `this PDF guide <USER/sph/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above.

* D diffusion coefficient (length\^2/time units)
* h kernel function cutoff (distance units)


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This style does not support the :doc:`pair\_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.  Thus, you need to re-specify the pair\_style and
pair\_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run\_style respa <run_style>` command.  It does not support the *inner*\ ,
*middle*\ , *outer* keywords.

Restrictions
""""""""""""


This pair style is part of the USER-SPH package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, pair\_sph/rhosum

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
