.. index:: compute momentum

compute momentum command
========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID momentum

* ID, group-ID are documented in :doc:`compute <compute>` command
* momentum = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all momentum

Description
"""""""""""

Define a computation that calculates the translational momentum
of a group of particles.

The momentum of each particles is computed as m v, where m and v are
the mass and velocity of the particle.

**Output info:**

This compute calculates a global vector (the summed momentum) of
length 3. This value can be used by any command that uses a global
vector value from a compute as input. See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The vector value calculated by this compute is "extensive". The vector
value will be in mass\*velocity :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

**Default:** none
