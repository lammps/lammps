.. index:: compute temp/region

compute temp/region command
===========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID temp/region region-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/region = style name of this compute command
* region-ID = ID of region to use for choosing atoms

Examples
""""""""


.. parsed-literal::

   compute mine flow temp/region boundary

Description
"""""""""""

Define a computation that calculates the temperature of a group of
atoms in a geometric region.  This can be useful for thermostatting
one portion of the simulation box.  E.g. a McDLT simulation where one
side is cooled, and the other side is heated.  A compute of this style
can be used by any command that computes a temperature,
e.g. :doc:`thermo_modify <thermo_modify>`, :doc:`fix temp/rescale <fix_temp_rescale>`, etc.

Note that a *region*\ -style temperature can be used to thermostat with
:doc:`fix temp/rescale <fix_temp_rescale>` or :doc:`fix langevin <fix_langevin>`, but should probably not be used with
Nose/Hoover style fixes (:doc:`fix nvt <fix_nh>`, :doc:`fix npt <fix_nh>`, or :doc:`fix nph <fix_nh>`), if the
degrees-of-freedom included in the computed T varies with time.

The temperature is calculated by the formula KE = dim/2 N k T, where
KE = total kinetic energy of the group of atoms (sum of 1/2 m v\^2),
dim = 2 or 3 = dimensionality of the simulation, N = number of atoms
in both the group and region, k = Boltzmann constant, and T =
temperature.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that v\^2 is replaced by vx\*vy for the xy
component, etc.  The 6 components of the vector are ordered xx, yy,
zz, xy, xz, yz.

The number of atoms contributing to the temperature is calculated each
time the temperature is evaluated since it is assumed atoms can
enter/leave the region.  Thus there is no need to use the *dynamic*
option of the :doc:`compute_modify <compute_modify>` command for this
compute style.

The removal of atoms outside the region by this fix is essentially
computing the temperature after a "bias" has been removed, which in
this case is the velocity of any atoms outside the region.  If this
compute is used with a fix command that performs thermostatting then
this bias will be subtracted from each atom, thermostatting of the
remaining thermal velocity will be performed, and the bias will be
added back in.  Thermostatting fixes that work in this way include
:doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix temp/berendsen <fix_temp_berendsen>`, and :doc:`fix langevin <fix_langevin>`.  This means that when this compute
is used to calculate the temperature for any of the thermostatting
fixes via the :doc:`fix modify temp <fix_modify>` command, the thermostat
will operate only on atoms that are currently in the geometric
region.

Unlike other compute styles that calculate temperature, this compute
does not subtract out degrees-of-freedom due to fixes that constrain
motion, such as :doc:`fix shake <fix_shake>` and :doc:`fix rigid <fix_rigid>`.  This is because those degrees of freedom
(e.g. a constrained bond) could apply to sets of atoms that straddle
the region boundary, and hence the concept is somewhat ill-defined.
If needed the number of subtracted degrees-of-freedom can be set
explicitly using the *extra* option of the
:doc:`compute_modify <compute_modify>` command.

See the :doc:`Howto thermostat <Howto_thermostat>` doc page for a
discussion of different ways to compute temperature and perform
thermostatting.

**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1-6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`, :doc:`compute pressure <compute_pressure>`

**Default:** none
