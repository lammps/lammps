.. index:: compute erotate/rigid

compute erotate/rigid command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID erotate/rigid fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* erotate/rigid = style name of this compute command
* fix-ID = ID of rigid body fix

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all erotate/rigid myRigid

Description
"""""""""""

Define a computation that calculates the rotational kinetic energy of
a collection of rigid bodies, as defined by one of the
:doc:`fix rigid <fix_rigid>` command variants.

The rotational energy of each rigid body is computed as
:math:`\frac12 I \omega_\text{body}^2`,
where :math:`I` is the inertia tensor for the rigid body and
:math:`\omega_\text{body}` is its angular velocity vector.
Both :math:`I` and :math:`\omega_\text{body}` are in the frame of
reference of the rigid body (i.e., :math:`I` is diagonal).

The *fix-ID* should be the ID of one of the :doc:`fix rigid <fix_rigid>`
commands which defines the rigid bodies.  The group specified in the
compute command is ignored.  The rotational energy of all the rigid
bodies defined by the fix rigid command in included in the calculation.

Output info
"""""""""""

This compute calculates a global scalar (the summed rotational energy
of all the rigid bodies).  This value can be used by any command that
uses a global scalar value from a compute as input.  See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the RIGID package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute ke/rigid <compute_ke_rigid>`

Default
"""""""

none
