.. index:: compute temp
.. index:: compute temp/kk

compute temp command
====================

Accelerator Variants: *temp/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all temp
   compute myTemp mobile temp

Description
"""""""""""

Define a computation that calculates the temperature of a group of
atoms.  A compute of this style can be used by any command that
computes a temperature, e.g. :doc:`thermo_modify <thermo_modify>`,
:doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix npt <fix_nh>`,
etc.

The temperature is calculated by the formula

.. math::

   \text{KE} = \frac{\text{dim}}{2} N k_B T,

where KE = total kinetic energy of the group of atoms (sum of
:math:`\frac12 m v^2`), dim = 2 or 3 is the dimensionality of the
simulation, :math:`N` is the number of atoms in the group, :math:`k_B`
is the Boltzmann constant, and :math:`T` is the absolute temperature.

A kinetic energy tensor, stored as a six-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that :math:`v^2` is replaced by :math:`v_x
v_y` for the :math:`xy` component, and so on.  The six components of
the vector are ordered :math:`xx`, :math:`yy`, :math:`zz`, :math:`xy`,
:math:`xz`, :math:`yz`.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

This compute subtracts out degrees-of-freedom due to fixes that
constrain molecular motion, such as :doc:`fix shake <fix_shake>` and
:doc:`fix rigid <fix_rigid>`.  This means the temperature of groups of
atoms that include these constraints will be computed correctly.  If
needed, the subtracted degrees-of-freedom can be altered using the
*extra* option of the :doc:`compute_modify <compute_modify>` command.

A compute of this style with the ID of "thermo_temp" is created when
LAMMPS starts up, as if this command were in the input script:

.. code-block:: LAMMPS

   compute thermo_temp all temp

See the "thermo_style" command for more details.

See the :doc:`Howto thermostat <Howto_thermostat>` page for a
discussion of different ways to compute temperature and perform
thermostatting.

----------

.. include:: accel_styles.rst

----------

Output info
"""""""""""

This compute calculates a global scalar (the temperature) and a global
vector of length six (KE tensor), which can be accessed by indices
1--6.  These values can be used by any command that uses global scalar
or vector values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp/partial <compute_temp_partial>`,
:doc:`compute temp/region <compute_temp_region>`,
:doc:`compute pressure <compute_pressure>`

Default
"""""""

none
