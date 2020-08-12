.. index:: compute cac/nodal/temp

compute cac/nodal/temp command
==============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID temp

* ID, group-ID are documented in :doc:`compute <compute>` command
* cac/nodal/temp = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all cac/nodal/temp
   compute myTemp mobile cac/nodal/temp

Description
"""""""""""

Define a computation that calculates the kinetic temperature of a group of
atoms/elements by using their velocities and nodal velocities.
The temperature is calculated by the formula KE = dim/2 N k T, where
KE = total kinetic energy of the group of atoms/elements (sum of 1/2 m v\^2),
dim = 2 or 3 = dimensionality of the simulation, N = number of atoms
in the group, k = Boltzmann constant, and T = temperature.

The number of atoms/elements contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute_modify <compute_modify>` command if this is not the case.

This compute subtracts out degrees-of-freedom due to fixes that
constrain molecular motion, such as :doc:`fix shake <fix_shake>` and
:doc:`fix rigid <fix_rigid>`.  This means the temperature of groups 
that include these constraints will be computed correctly.  If
needed, the subtracted degrees-of-freedom can be altered using the
*extra* option of the :doc:`compute_modify <compute_modify>` command.

----------

**Output info:**

This compute calculates a global scalar (the temperature).
This value can be used by any command that uses a global scalar 
as input. See the :doc:`Howto output <Howto_output>` doc page for
an overview of LAMMPS output options.

The scalar value calculated by this compute is "intensive".
The scalar value will be in temperature :doc:`units <units>`.

Restrictions
""""""""""""

This compute requires a CAC atom style

**Default:** none
