Calculate thermal conductivity
==============================

The thermal conductivity kappa of a material can be measured in at
least 4 ways using various options in LAMMPS.  See the examples/KAPPA
directory for scripts that implement the 4 methods discussed here for
a simple Lennard-Jones fluid model.  Also, see the :doc:`Howto viscosity <Howto_viscosity>` page for an analogous discussion
for viscosity.

The thermal conductivity tensor kappa is a measure of the propensity
of a material to transmit heat energy in a diffusive manner as given
by Fourier's law

J = -kappa grad(T)

where J is the heat flux in units of energy per area per time and
grad(T) is the spatial gradient of temperature.  The thermal
conductivity thus has units of energy per distance per time per degree
K and is often approximated as an isotropic quantity, i.e. as a
scalar.

The first method is to setup two thermostatted regions at opposite
ends of a simulation box, or one in the middle and one at the end of a
periodic box.  By holding the two regions at different temperatures
with a :doc:`thermostatting fix <Howto_thermostat>`, the energy added to
the hot region should equal the energy subtracted from the cold region
and be proportional to the heat flux moving between the regions.  See
the papers by :ref:`Ikeshoji and Hafskjold <howto-Ikeshoji>` and
:ref:`Wirnsberger et al <howto-Wirnsberger>` for details of this idea.  Note
that thermostatting fixes such as :doc:`fix nvt <fix_nh>`, :doc:`fix langevin <fix_langevin>`, and :doc:`fix temp/rescale <fix_temp_rescale>` store the cumulative energy they
add/subtract.

Alternatively, as a second method, the :doc:`fix heat <fix_heat>` or
:doc:`fix ehex <fix_ehex>` commands can be used in place of thermostats
on each of two regions to add/subtract specified amounts of energy to
both regions.  In both cases, the resulting temperatures of the two
regions can be monitored with the "compute temp/region" command and
the temperature profile of the intermediate region can be monitored
with the :doc:`fix ave/chunk <fix_ave_chunk>` and :doc:`compute ke/atom <compute_ke_atom>` commands.

The third method is to perform a reverse non-equilibrium MD simulation
using the :doc:`fix thermal/conductivity <fix_thermal_conductivity>`
command which implements the rNEMD algorithm of Muller-Plathe.
Kinetic energy is swapped between atoms in two different layers of the
simulation box.  This induces a temperature gradient between the two
layers which can be monitored with the :doc:`fix ave/chunk <fix_ave_chunk>` and :doc:`compute ke/atom <compute_ke_atom>` commands.  The fix tallies the
cumulative energy transfer that it performs.  See the :doc:`fix thermal/conductivity <fix_thermal_conductivity>` command for
details.

The fourth method is based on the Green-Kubo (GK) formula which
relates the ensemble average of the auto-correlation of the heat flux
to kappa.  The heat flux can be calculated from the fluctuations of
per-atom potential and kinetic energies and per-atom stress tensor in
a steady-state equilibrated simulation.  This is in contrast to the
two preceding non-equilibrium methods, where energy flows continuously
between hot and cold regions of the simulation box.

The :doc:`compute heat/flux <compute_heat_flux>` command can calculate
the needed heat flux and describes how to implement the Green_Kubo
formalism using additional LAMMPS commands, such as the :doc:`fix ave/correlate <fix_ave_correlate>` command to calculate the needed
auto-correlation.  See the page for the :doc:`compute heat/flux <compute_heat_flux>` command for an example input script
that calculates the thermal conductivity of solid Ar via the GK
formalism.

----------

.. _howto-Ikeshoji:

**(Ikeshoji)** Ikeshoji and Hafskjold, Molecular Physics, 81, 251-261
(1994).

.. _howto-Wirnsberger:

**(Wirnsberger)** Wirnsberger, Frenkel, and Dellago, J Chem Phys, 143, 124104
(2015).
