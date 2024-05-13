Barostats
=========

Barostatting means controlling the pressure in an MD simulation.
:doc:`Thermostatting <Howto_thermostat>` means controlling the
temperature of the particles.  Since the pressure includes a kinetic
component due to particle velocities, both these operations require
calculation of the temperature.  Typically a target temperature (T)
and/or pressure (P) is specified by the user, and the thermostat or
barostat attempts to equilibrate the system to the requested T and/or
P.

Barostatting in LAMMPS is performed by :doc:`fixes <fix>`.  Two
barostatting methods are currently available: Nose-Hoover (npt and
nph) and Berendsen:

* :doc:`fix npt <fix_nh>`
* :doc:`fix npt/sphere <fix_npt_sphere>`
* :doc:`fix npt/asphere <fix_npt_asphere>`
* :doc:`fix nph <fix_nh>`
* :doc:`fix press/berendsen <fix_press_berendsen>`

The :doc:`fix npt <fix_nh>` commands include a Nose-Hoover thermostat
and barostat.  :doc:`Fix nph <fix_nh>` is just a Nose/Hoover barostat;
it does no thermostatting.  Both :doc:`fix nph <fix_nh>` and :doc:`fix press/berendsen <fix_press_berendsen>` can be used in conjunction
with any of the thermostatting fixes.

As with the :doc:`thermostats <Howto_thermostat>`, :doc:`fix npt <fix_nh>`
and :doc:`fix nph <fix_nh>` only use translational motion of the
particles in computing T and P and performing thermo/barostatting.
:doc:`Fix npt/sphere <fix_npt_sphere>` and :doc:`fix npt/asphere <fix_npt_asphere>` thermo/barostat using not only
translation velocities but also rotational velocities for spherical
and aspherical particles.

All of the barostatting fixes use the :doc:`compute pressure <compute_pressure>` compute to calculate a current
pressure.  By default, this compute is created with a simple :doc:`compute temp <compute_temp>` (see the last argument of the :doc:`compute pressure <compute_pressure>` command), which is used to calculated
the kinetic component of the pressure.  The barostatting fixes can
also use temperature computes that remove bias for the purpose of
computing the kinetic component which contributes to the current
pressure.  See the doc pages for the individual fixes and for the
:doc:`fix_modify <fix_modify>` command for instructions on how to assign
a temperature or pressure compute to a barostatting fix.

.. note::

   As with the thermostats, the Nose/Hoover methods (:doc:`fix npt <fix_nh>` and :doc:`fix nph <fix_nh>`) perform time integration.
   :doc:`Fix press/berendsen <fix_press_berendsen>` does NOT, so it should
   be used with one of the constant NVE fixes or with one of the NVT
   fixes.

Thermodynamic output, which can be setup via the
:doc:`thermo_style <thermo_style>` command, often includes pressure
values.  As explained on the page for the
:doc:`thermo_style <thermo_style>` command, the default pressure is
setup by the thermo command itself.  It is NOT the pressure associated
with any barostatting fix you have defined or with any compute you
have defined that calculates a pressure.  The doc pages for the
barostatting fixes explain the ID of the pressure compute they create.
Thus if you want to view these pressures, you need to specify them
explicitly via the :doc:`thermo_style custom <thermo_style>` command.
Or you can use the :doc:`thermo_modify <thermo_modify>` command to
re-define what pressure compute is used for default thermodynamic
output.
