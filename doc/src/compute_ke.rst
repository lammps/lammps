.. index:: compute ke

compute ke command
==================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID ke

* ID, group-ID are documented in :doc:`compute <compute>` command
* ke = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all ke

Description
"""""""""""

Define a computation that calculates the translational kinetic energy
of a group of particles.

The kinetic energy of each particle is computed as :math:`\frac{1}{2} m
v^2`, where *m* and *v* are the mass and velocity of the particle.

There is a subtle difference between the quantity calculated by this
compute and the kinetic energy calculated by the *ke* or *etotal*
keyword used in thermodynamic output, as specified by the
:doc:`thermo_style <thermo_style>` command.  For this compute, kinetic
energy is "translational" kinetic energy, calculated by the simple
formula above.  For thermodynamic output, the *ke* keyword infers
kinetic energy from the temperature of the system with
:math:`\frac{1}{2} k_B T` of energy for each degree of freedom.  For the
default temperature computation via the :doc:`compute temp
<compute_temp>` command, these are the same.  But different computes
that calculate temperature can subtract out different non-thermal
components of velocity and/or include different degrees of freedom
(translational, rotational, etc).

**Output info:**

This compute calculates a global scalar (the summed KE).  This value
can be used by any command that uses a global scalar value from a
compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute erotate/sphere <compute_erotate_sphere>`

**Default:** none
