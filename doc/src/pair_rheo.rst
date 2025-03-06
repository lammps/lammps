.. index:: pair_style rheo

pair_style rheo command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style rheo cutoff keyword values

* cutoff = global cutoff for kernel (distance units)
* zero or more keyword/value pairs may be appended to args
* keyword = *rho/damp* or *artificial/visc* or *harmonic/means*

.. parsed-literal::

     *rho/damp* args = density damping prefactor :math:`\xi`
     *artificial/visc* args = artificial viscosity prefactor :math:`\zeta`
     *harmonic/means* args = none

Examples
""""""""

.. code-block:: LAMMPS

   pair_style rheo 3.0 rho/damp 1.0 artificial/visc 2.0
   pair_coeff * *

Description
"""""""""""

.. versionadded:: 29Aug2024

Pair style *rheo* computes pressure and viscous forces between particles
in the :doc:`rheo package <Howto_rheo>`. If thermal evolution is turned
on in :doc:`fix rheo <fix_rheo>`, then the pair style also calculates
heat exchanged between particles.

The *artificial/viscosity* keyword is used to specify the magnitude
:math:`\zeta` of an optional artificial viscosity contribution to forces.
This factor can help stabilize simulations by smoothing out small length
scale variations in velocity fields. Artificial viscous forces typically
are only exchanged by fluid particles. However, if interfaces are not
reconstructed in fix rheo, fluid particles will also exchange artificial
viscous forces with solid particles to improve stability.

The *rho/damp* keyword is used to specify the magnitude :math:`\xi` of
an optional pairwise damping term between the density of particles. This
factor can help stabilize simulations by smoothing out small length
scale variations in density fields. However, in systems that develop
a density gradient in equilibrium (e.g. in a hydrostatic column underlying
gravity), this option may be inappropriate.

If particles have different viscosities or conductivities, the
*harmonic/means* keyword changes how they are averaged before calculating
pairwise forces or heat exchanges. By default, an arithmetic averaged is
used, however, a harmonic mean may improve stability in systems with multiple
fluid phases with large disparities in viscosities.

No coefficients are defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples
above.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.
Thus, you need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the *inner*,
*middle*, *outer* keywords.

Restrictions
""""""""""""

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo <fix_rheo>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

Density damping and artificial viscous forces are not calculated.
Arithmetic means are used for mixing particle properties.
