.. index:: compute mspin/distance

compute mspin/distance command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID mspin/distance fix-ID body1 body2

* ID, group-ID are documented in :doc:`compute <compute>` command
* mspin/distance = style name of this compute command
* fix-ID = ID of rigid/nvt/mspin fix
* body1 = rigid body id of the first nanoparticle
* body2 = rigid body id of the second nanoparticle

Examples
""""""""

.. code-block:: LAMMPS

   compute   1   all mspin/distance fcore 1 2

Description
"""""""""""

Define a computation that calculates the distance between the center-of-masses of
two rigid nanoparticle cores as defined by the
:doc:`fix rigid/nvt/mspin <fix_rigid_mspin>` command.

The *fix-ID* should be the ID of the :doc:`fix rigid/nvt/mspin <fix_rigid_mspin>`
command which defines the rigid bodies. The group specified in the
compute command is ignored.

Output info
"""""""""""

This compute calculates a global scalar (the interparticle distance).
This value can be used by any command that uses a global scalar value from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output options.

The scalar value will be in distance :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the :ref:`MSPIN <PKG-MSPIN>`. It is only enabled if
LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute mspin/energy <compute_mspin_energy>`

Default
"""""""

none
