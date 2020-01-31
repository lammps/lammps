.. index:: fix nvt/sllod

fix nvt/sllod command
=====================

fix nvt/sllod/intel command
===========================

fix nvt/sllod/omp command
=========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID nvt/sllod keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvt/sllod = style name of this fix command
* additional thermostat related keyword/value pairs from the :doc:`fix nvt <fix_nh>` command can be appended

Examples
""""""""


.. parsed-literal::

   fix 1 all nvt/sllod temp 300.0 300.0 100.0
   fix 1 all nvt/sllod temp 300.0 300.0 100.0 drag 0.2

Description
"""""""""""

Perform constant NVT integration to update positions and velocities
each timestep for atoms in the group using a Nose/Hoover temperature
thermostat.  V is volume; T is temperature.  This creates a system
trajectory consistent with the canonical ensemble.

This thermostat is used for a simulation box that is changing size
and/or shape, for example in a non-equilibrium MD (NEMD) simulation.
The size/shape change is induced by use of the :doc:`fix deform <fix_deform>` command, so each point in the simulation box
can be thought of as having a "streaming" velocity.  This
position-dependent streaming velocity is subtracted from each atom's
actual velocity to yield a thermal velocity which is used for
temperature computation and thermostatting.  For example, if the box
is being sheared in x, relative to y, then points at the bottom of the
box (low y) have a small x velocity, while points at the top of the
box (hi y) have a large x velocity.  These velocities do not
contribute to the thermal "temperature" of the atom.

.. note::

   :doc:`Fix deform <fix_deform>` has an option for remapping either
   atom coordinates or velocities to the changing simulation box.  To use
   fix nvt/sllod, fix deform should NOT remap atom positions, because fix
   nvt/sllod adjusts the atom positions and velocities to create a
   velocity profile that matches the changing box size/shape.  Fix deform
   SHOULD remap atom velocities when atoms cross periodic boundaries
   since that is consistent with maintaining the velocity profile created
   by fix nvt/sllod.  LAMMPS will give an error if this setting is not
   consistent.

The SLLOD equations of motion, originally proposed by Hoover and Ladd
(see :ref:`(Evans and Morriss) <Evans3>`), were proven to be equivalent to
Newton's equations of motion for shear flow by :ref:`(Evans and Morriss) <Evans3>`. They were later shown to generate the desired
velocity gradient and the correct production of work by stresses for
all forms of homogeneous flow by :ref:`(Daivis and Todd) <Daivis>`.  As
implemented in LAMMPS, they are coupled to a Nose/Hoover chain
thermostat in a velocity Verlet formulation, closely following the
implementation used for the :doc:`fix nvt <fix_nh>` command.

.. note::

   A recent (2017) book by :ref:`(Daivis and Todd) <Daivis-sllod>`
   discusses use of the SLLOD method and non-equilibrium MD (NEMD)
   thermostatting generally, for both simple and complex fluids,
   e.g. molecular systems.  The latter can be tricky to do correctly.

Additional parameters affecting the thermostat are specified by
keywords and values documented with the :doc:`fix nvt <fix_nh>`
command.  See, for example, discussion of the *temp* and *drag*
keywords.

This fix computes a temperature each timestep.  To do this, the fix
creates its own compute of style "temp/deform", as if this command had
been issued:


.. parsed-literal::

   compute fix-ID_temp group-ID temp/deform

See the :doc:`compute temp/deform <compute_temp_deform>` command for
details.  Note that the ID of the new compute is the fix-ID +
underscore + "temp", and the group for the new compute is the same as
the fix group.

Note that this is NOT the compute used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo\_temp*.
This means you can change the attributes of this fix's temperature
(e.g. its degrees-of-freedom) via the
:doc:`compute_modify <compute_modify>` command or print this temperature
during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* will have no
effect on this fix.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that calculate a temperature
after removing a "bias" from the atom velocities.  E.g. removing the
center-of-mass velocity from a group of atoms or only calculating
temperature on the x-component of velocity or only calculating
temperature for atoms in a geometric region.  This is not done by
default, but only if the :doc:`fix_modify <fix_modify>` command is used
to assign a temperature compute to this fix that includes such a bias
term.  See the doc pages for individual :doc:`compute commands <compute>` to determine which ones include a bias.  In
this case, the thermostat works in the following manner: the current
temperature is calculated taking the bias into account, bias is
removed from each atom, thermostatting is performed on the remaining
thermal degrees of freedom, and the bias is added back in.


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

This fix writes the state of the Nose/Hoover thermostat to :doc:`binary restart files <restart>`.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a :doc:`compute <compute>` you have
defined to this fix which will be used in its thermostatting
procedure.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Nose/Hoover thermostatting to
the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix nvt <fix_nh>` command.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix works best without Nose-Hoover chain thermostats, i.e. using
tchain = 1.  Setting tchain to larger values can result in poor
equilibration.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix langevin <fix_langevin>`,
:doc:`fix_modify <fix_modify>`, :doc:`compute temp/deform <compute_temp_deform>`

Default
"""""""

Same as :doc:`fix nvt <fix_nh>`, except tchain = 1.


----------


.. _Evans3:



**(Evans and Morriss)** Evans and Morriss, Phys Rev A, 30, 1528 (1984).

.. _Daivis:



**(Daivis and Todd)** Daivis and Todd, J Chem Phys, 124, 194103 (2006).

.. _Daivis-sllod:



**(Daivis and Todd)** Daivis and Todd, Nonequilibrium Molecular Dynamics (book),
Cambridge University Press, https://doi.org/10.1017/9781139017848, (2017).


