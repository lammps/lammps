.. index:: fix rheo

fix rheo command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID rheo cut kstyle keyword values...

* ID, group-ID are documented in :doc:`fix <fix>` command
* rheo = style name of this fix command
* cut = *quintic* or *RK0* or *RK1* or *RK2*
* zero or more keyword/value pairs may be appended to args
* keyword = *shift* or *thermal* or *surface/detection* or *interface/reconstruction* or *rho/sum* or *density* or *sound/squared*

  .. parsed-literal::

       *shift* values = none, turns on velocity shifting
       *thermal* values = none, turns on thermal evolution
       *surface/detection* values = *sdstyle* *limit* *limit/splash*
         *sdstyle* = *coordination* or *divergence*
         *limit* = threshold for surface particles (unitless)
         *limit/splash* = threshold for splash particles (unitless)
       *interface/reconstruct* values = none, reconstructs interfaces with solid particles
       *rho/sum* values = none, uses the kernel to compute the density of particles
       *density* values = *rho0* (density)
       *sound/squared* values = *csq* (velocity\^2)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all rheo 1.0 quintic thermal density 0.1 sound/squared 10.0
   fix 1 all rheo 1.0 RK1 shift surface/detection coordination 40

Description
"""""""""""

Fix description...

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix must be used with an atom style that includes density
such as atom_style rheo or rheo/thermal. This fix must be used in
conjuction with :doc:`fix rheo/pressure <fix_rheo_pressure>`. and
:doc:`fix rheo/viscosity <fix_rheo_viscosity>`, If the *thermal*
setting is used, there must also be an instance of
:doc:`fix rheo/thermal <fix_rheo_thermal>`. The fix group must be
set to all.

This fix is part of the RHEO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix rheo/viscosity <fix_rheo_viscosity>`,
:doc:`fix rheo/pressure <fix_rheo_pressure>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`,
:doc:`pair rheo <pair_rheo>`,
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`

Default
"""""""

*rho0* and *csq* are set to 1.0.
