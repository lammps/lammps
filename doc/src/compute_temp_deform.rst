.. index:: compute temp/deform
.. index:: compute temp/deform/kk

compute temp/deform command
===========================

Accelerator Variants: *temp/deform/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/deform

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/deform = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute myTemp all temp/deform

Description
"""""""""""

Define a computation that calculates the temperature of a group of
atoms, after subtracting out a streaming velocity induced by the
simulation box changing size and/or shape, for example in a
non-equilibrium MD (NEMD) simulation.  The size/shape change is
induced by use of the :doc:`fix deform <fix_deform>` command.  A compute
of this style is created by the :doc:`fix nvt/sllod <fix_nvt_sllod>`
command to compute the thermal temperature of atoms for thermostatting
purposes.  A compute of this style can also be used by any command
that computes a temperature (e.g., :doc:`thermo_modify <thermo_modify>`,
:doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix npt <fix_nh>`).

The deformation fix changes the box size and/or shape over time, so
each atom in the simulation box can be thought of as having a
"streaming" velocity.  For example, if the box is being sheared in *x*,
relative to *y*, then atoms at the bottom of the box (low *y*) have a
small *x* velocity, while atoms at the top of the box (high *y*) have a
large *x* velocity.  This position-dependent streaming velocity is
subtracted from each atom's actual velocity to yield a thermal
velocity, which is then used to compute the temperature.

.. note::

   :doc:`Fix deform <fix_deform>` has an option for remapping either
   atom coordinates or velocities to the changing simulation box.  When
   using this compute in conjunction with a deforming box, fix deform
   should NOT remap atom positions, but rather should let atoms respond
   to the changing box by adjusting their own velocities (or let
   :doc:`fix deform <fix_deform>` remap the atom velocities; see its remap
   option).  If fix deform does remap atom positions, then they appear to
   move with the box but their velocity is not changed, and thus they do
   NOT have the streaming velocity assumed by this compute.  LAMMPS will
   warn you if fix deform is defined and its remap setting is not
   consistent with this compute.

After the streaming velocity has been subtracted from each atom, the
temperature is calculated by the formula

.. math::

   \text{KE} = \frac{\text{dim}}{2} N k_B T,

where KE is the total kinetic energy of the group of atoms (sum of
:math:`\frac12 m v^2`, dim = 2 or 3 is the dimensionality of the
simulation, :math:`N` is the number of atoms in the group, :math:`k_B`
is the Boltzmann constant, and :math:`T` is the temperature.  Note that
:math:`v` in the kinetic energy formula is the atom's velocity.

A kinetic energy tensor, stored as a six-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that :math:`v^2` is replaced by :math:`v_x v_y` for
the :math:`xy` component, and so on. The six components of the vector are
ordered :math:`xx`, :math:`yy`, :math:`zz`, :math:`xy`, :math:`xz`,
:math:`yz`.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

The removal of the box deformation velocity component by this fix is
essentially computing the temperature after a "bias" has been removed
from the velocity of the atoms.  If this compute is used with a fix
command that performs thermostatting then this bias will be subtracted
from each atom, thermostatting of the remaining thermal velocity will
be performed, and the bias will be added back in.  Thermostatting
fixes that work in this way include :doc:`fix nvt <fix_nh>`,
:doc:`fix temp/rescale <fix_temp_rescale>`,
:doc:`fix temp/berendsen <fix_temp_berendsen>`, and
:doc:`fix langevin <fix_langevin>`.

.. note::

   The temperature calculated by this compute is only accurate if
   the atoms are indeed moving with a stream velocity profile that
   matches the box deformation.  If not, then the compute will subtract
   off an incorrect stream velocity, yielding a bogus thermal
   temperature.  You should **not** assume that your atoms are streaming at
   the same rate the box is deforming.  Rather, you should monitor their
   velocity profiles (e.g., via the :doc:`fix ave/chunk <fix_ave_chunk>`
   command).  You can also compare the results of this compute to
   :doc:`compute temp/profile <compute_temp_profile>`, which actually
   calculates the stream profile before subtracting it.  If the two computes do
   not give roughly the same temperature, then your atoms are not streaming
   consistent with the box deformation.  See the :doc:`fix deform <fix_deform>`
   command for more details on ways to get atoms to stream consistently with
   the box deformation.

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

The scalar value will be in temperature :doc:`units <units>`.
The vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp/ramp <compute_temp_ramp>`,
:doc:`compute temp/profile <compute_temp_profile>`,
:doc:`fix deform <fix_deform>`,
:doc:`fix nvt/sllod <fix_nvt_sllod>`

Default
"""""""

none
