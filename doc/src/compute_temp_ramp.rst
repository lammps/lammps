.. index:: compute temp/ramp

compute temp/ramp command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/ramp vdim vlo vhi dim clo chi keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/ramp = style name of this compute command
* vdim = *vx* or *vy* or *vz*
* vlo,vhi = subtract velocities between vlo and vhi (velocity units)
* dim = *x* or *y* or *z*
* clo,chi = lower and upper bound of domain to subtract from (distance units)
* zero or more keyword/value pairs may be appended
* keyword = *units*

.. parsed-literal::

     *units* value = *lattice* or *box*

Examples
""""""""

.. code-block:: LAMMPS

   compute 2nd middle temp/ramp vx 0 8 y 2 12 units lattice

Description
"""""""""""

Define a computation that calculates the temperature of a group of atoms,
after subtracting out an ramped velocity profile before computing the kinetic
energy.  A compute of this style can be used by any command that computes a
temperature (e.g. :doc:`thermo_modify <thermo_modify>`,
:doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix npt <fix_nh>`).

The meaning of the arguments for this command which define the
velocity ramp are the same as for the :doc:`velocity ramp <velocity>`
command which was presumably used to impose the velocity.

After the ramp velocity has been subtracted from the specified
dimension for each atom, the temperature is calculated by the formula

.. math::

   \text{KE} = \frac{\text{dim}}{2} N k_B T,

where KE is the total kinetic energy of the group of atoms (sum of
:math:`\frac12 m v^2`), dim = 2 or 3 is the dimensionality of the simulation,
:math:`N` is the number of atoms in the group, :math:`k_B` is the Boltzmann
constant, and :math:`T` is the absolute temperature.

The *units* keyword determines the meaning of the distance units used
for coordinates (*clo*, *chi*) and velocities (*vlo*, *vhi*).  A *box* value
selects standard distance units as defined by the :doc:`units <units>`
command (e.g., :math:`\mathrm{\mathring{A}}` for units = real or metal).  A
*lattice* value means the distance units are in lattice spacings (i.e.,
velocity in lattice spacings per unit time).  The :doc:`lattice <lattice>`
command must have been previously used to define the lattice spacing.

A kinetic energy tensor, stored as a six-element vector, is also calculated by
this compute for use in the computation of a pressure tensor.  The formula for
the components of the tensor is the same as the above formula, except that
:math:`v^2` is replaced by :math:`v_x v_y` for the :math:`xy` component, and
so on.  The six components of the vector are ordered :math:`xx`, :math:`yy`,
:math:`zz`, :math:`xy`, :math:`xz`, :math:`yz`.

The number of atoms contributing to the temperature is assumed to be constant
for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

The removal of the ramped velocity component by this fix is
essentially computing the temperature after a "bias" has been removed
from the velocity of the atoms.  If this compute is used with a fix
command that performs thermostatting then this bias will be subtracted
from each atom, thermostatting of the remaining thermal velocity will
be performed, and the bias will be added back in.  Thermostatting
fixes that work in this way include :doc:`fix nvt <fix_nh>`,
:doc:`fix temp/rescale <fix_temp_rescale>`,
:doc:`fix temp/berendsen <fix_temp_berendsen>`, and
:doc:`fix langevin <fix_langevin>`.

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

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1--6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
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

:doc:`compute temp <compute_temp>`, :doc:`compute temp/profie <compute_temp_profile>`, :doc:`compute temp/deform <compute_temp_deform>`, :doc:`compute pressure <compute_pressure>`

Default
"""""""

The option default is units = lattice.
