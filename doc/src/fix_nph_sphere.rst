.. index:: fix nph/sphere

fix nph/sphere command
======================

fix nph/sphere/omp command
==========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nph/sphere args keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nph/sphere = style name of this fix command
* keyword = *disc*
  
  .. parsed-literal::
  
       *disc* value = none = treat particles as 2d discs, not spheres

* additional barostat related keyword/value pairs from the :doc:`fix nph <fix_nh>` command can be appended

Examples
""""""""


.. parsed-literal::

   fix 1 all nph/sphere iso 0.0 0.0 1000.0
   fix 2 all nph/sphere x 5.0 5.0 1000.0
   fix 2 all nph/sphere x 5.0 5.0 1000.0 disc
   fix 2 all nph/sphere x 5.0 5.0 1000.0 drag 0.2
   fix 2 water nph/sphere aniso 0.0 0.0 1000.0 dilate partial

Description
"""""""""""

Perform constant NPH integration to update position, velocity, and
angular velocity each timestep for finite-size spherical particles in
the group using a Nose/Hoover pressure barostat.  P is pressure; H is
enthalpy.  This creates a system trajectory consistent with the
isenthalpic ensemble.

This fix differs from the :doc:`fix nph <fix_nh>` command, which assumes
point particles and only updates their position and velocity.

If the *disc* keyword is used, then each particle is treated as a 2d
disc (circle) instead of as a sphere.  This is only possible for 2d
simulations, as defined by the :doc:`dimension <dimension>` keyword.
The only difference between discs and spheres in this context is their
moment of inertia, as used in the time integration.

Additional parameters affecting the barostat are specified by keywords
and values documented with the :doc:`fix nph <fix_nh>` command.  See,
for example, discussion of the *aniso*\ , and *dilate* keywords.

The particles in the fix group are the only ones whose velocities and
positions are updated by the velocity/position update portion of the
NPH integration.

Regardless of what particles are in the fix group, a global pressure is
computed for all particles.  Similarly, when the size of the simulation
box is changed, all particles are re-scaled to new positions, unless the
keyword *dilate* is specified with a value of *partial*\ , in which case
only the particles in the fix group are re-scaled.  The latter can be
useful for leaving the coordinates of particles in a solid substrate
unchanged and controlling the pressure of a surrounding fluid.


----------


This fix computes a temperature and pressure each timestep.  To do
this, the fix creates its own computes of style "temp/sphere" and
"pressure", as if these commands had been issued:


.. parsed-literal::

   compute fix-ID_temp all temp/sphere
   compute fix-ID_press all pressure fix-ID_temp

See the :doc:`compute temp/sphere <compute_temp_sphere>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + underscore + "temp" or fix\_ID
+ underscore + "press", and the group for the new computes is "all"
since pressure is computed for the entire system.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo\_temp*
and *thermo\_press*.  This means you can change the attributes of this
fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command or print this temperature
or pressure during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* or
*thermo\_press* will have no effect on this fix.


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

**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the state of the Nose/Hoover barostat to :doc:`binary restart files <restart>`.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by this fix.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its thermostatting or barostatting procedure.  If you do this, note
that the kinetic energy derived from the compute temperature should be
consistent with the virial term computed using all atoms for the
pressure.  LAMMPS will warn you if you choose to compute temperature
on a subset of atoms.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Nose/Hoover barostatting to
the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix nph <fix_nh>` command.

This fix can ramp its target pressure over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix requires that atoms store torque and angular velocity (omega)
and a radius as defined by the :doc:`atom_style sphere <atom_style>`
command.

All particles in the group must be finite-size spheres.  They cannot
be point particles.

Use of the *disc* keyword is only allowed for 2d simulations, as
defined by the :doc:`dimension <dimension>` keyword.

Related commands
""""""""""""""""

:doc:`fix nph <fix_nh>`, :doc:`fix nve\_sphere <fix_nve_sphere>`, :doc:`fix nvt\_sphere <fix_nvt_sphere>`, :doc:`fix npt\_sphere <fix_npt_sphere>`,
:doc:`fix_modify <fix_modify>`

**Default:** none


