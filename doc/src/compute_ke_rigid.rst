.. index:: compute ke/rigid

compute ke/rigid command
========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID ke/rigid fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* ke = style name of this compute command
* fix-ID = ID of rigid body fix

Examples
""""""""


.. parsed-literal::

   compute 1 all ke/rigid myRigid

Description
"""""""""""

Define a computation that calculates the translational kinetic energy
of a collection of rigid bodies, as defined by one of the :doc:`fix rigid <fix_rigid>` command variants.

The kinetic energy of each rigid body is computed as 1/2 M Vcm\^2,
where M is the total mass of the rigid body, and Vcm is its
center-of-mass velocity.

The *fix-ID* should be the ID of one of the :doc:`fix rigid <fix_rigid>`
commands which defines the rigid bodies.  The group specified in the
compute command is ignored.  The kinetic energy of all the rigid
bodies defined by the fix rigid command in included in the
calculation.

**Output info:**

This compute calculates a global scalar (the summed KE of all the
rigid bodies).  This value can be used by any command that uses a
global scalar value from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the RIGID package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute erotate/rigid <compute_erotate_rigid>`

**Default:** none
