.. index:: compute temp/eff

compute temp/eff command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/eff

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/eff = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all temp/eff
   compute myTemp mobile temp/eff

Description
"""""""""""

Define a computation that calculates the temperature of a group of
nuclei and electrons in the :doc:`electron force field <pair_eff>`
model.  A compute of this style can be used by commands that compute a
temperature (e.g., :doc:`thermo_modify <thermo_modify>`,
:doc:`fix npt/eff <fix_nh_eff>`).

The temperature is calculated by the formula

.. math::

   \text{KE} = \frac{\text{dim}}{2} N k_B T,

where KE is the total kinetic energy of the group of atoms (sum of
:math:`\frac12 m v^2` for nuclei and sum of
:math:`\frac12 (m v^2 + \frac34 m s^2`) for electrons, where :math:`s`
includes the radial electron velocity contributions), dim = 2 or 3 is the
dimensionality of the simulation, :math:`N` is the number of atoms (only total
number of nuclei in the eFF (see the :doc:`pair_eff <pair_style>`
command) in the group, :math:`k_B` is the Boltzmann constant, and :math:`T` is
the absolute temperature.  This expression is summed over all nuclear and
electronic degrees of freedom, essentially by setting the kinetic contribution
to the heat capacity to :math:`\frac32 k` (where only nuclei contribute). This
subtlety is valid for temperatures well below the Fermi temperature, which for
densities two to five times the density of liquid hydrogen ranges from
86,000 to 170,000 K.

.. note::

   For eFF models, in order to override the default temperature
   reported by LAMMPS in the thermodynamic quantities reported via the
   :doc:`thermo <thermo>` command, the user should apply a
   :doc:`thermo_modify <thermo_modify>` command, as shown in the following
   example:

.. code-block:: LAMMPS

   compute         effTemp all temp/eff
   thermo_style    custom step etotal pe ke temp press
   thermo_modify   temp effTemp

A six-component kinetic energy tensor is also calculated by this compute
for use in the computation of a pressure tensor.  The formula for the
components of the tensor is the same as the above formula, except that
:math:`v^2` is replaced by :math:`v_x v_y` for the :math:`xy` component, etc.
For the eFF, again, the radial electronic velocities are also considered.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

This compute subtracts out degrees-of-freedom due to fixes that
constrain molecular motion, such as :doc:`fix shake <fix_shake>` and
:doc:`fix rigid <fix_rigid>`.  This means the temperature of groups of
atoms that include these constraints will be computed correctly.  If
needed, the subtracted degrees-of-freedom can be altered using the
*extra* option of the :doc:`compute_modify <compute_modify>` command.

See the :doc:`Howto thermostat <Howto_thermostat>` page for a
discussion of different ways to compute temperature and perform
thermostatting.

Output info
"""""""""""

The scalar value calculated by this compute is "intensive", meaning it
is independent of the number of atoms in the simulation.  The vector
values are "extensive", meaning they scale with the number of atoms in
the simulation.

Restrictions
""""""""""""

This compute is part of the EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute temp/partial <compute_temp_partial>`,
:doc:`compute temp/region <compute_temp_region>`,
:doc:`compute pressure <compute_pressure>`

Default
"""""""

none
