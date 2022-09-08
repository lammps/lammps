.. index:: compute mspin/energy

compute mspin/energy command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID mspin/energy fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* mspin/energy = style name of this compute command
* fix-ID = ID of rigid/nvt/mspin fix

Examples
""""""""

.. code-block:: LAMMPS

   compute   1   all mspin/energy npcore

Description
"""""""""""

Define a computation that calculates the Zeeman magnetic interaction energy and
magnetic dipolar interaction energy of a collection of rigid magnetic nanoparticles,
as defined by the :doc:`fix rigid/nvt/mspin <fix_rigid_mspin>` command.

The Zeeman energy of each nanoparticle is computed as

.. math::

    - \vec{\mu} \cdot \vec{B}

where :math:`\vec{\mu}` is the magnetic dipole moment of the rigid body, and
:math:`\vec{B}` is the external magnetic field applied to the system.

The dipolar interaction energy of each nanoparticle is computed as

.. math::
    
    - \alpha \frac{\mu_0}{4\pi} \frac{3\left(\vec{\mu_1}\cdot\vec{r}\right)\left(\vec{\mu_2}\cdot\vec{r}\right) - \vec{\mu_1}\cdot\vec{\mu_2}}{r^3}

where magnetic dipole moments are denoted by :math:`\vec{\mu}` and :math:`\vec{r}`
is a center-of-mass distance vector between the particles.

The *fix-ID* should be the ID of one of the :doc:`fix rigid/nvt/mspin <fix_rigid_mspin>`
commands which defines the rigid bodies. The group specified in the
compute command is ignored.  The energies of all the rigid
bodies defined by the *fix rigid/nvt/mspin* command are included in the
calculation.

Output info
"""""""""""

This compute calculates a global vector (the summed Zeeman energy and
the summed interparticle magnetic dipolar interaction energy of all rigid
nanoparticle cores). The values can be used by any command that uses a
global vector value from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output options.

The vector values calculated by this compute is "extensive." The
values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the :ref:`MSPIN <PKG-MSPIN>` package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute mspin/distance <compute_mspin_distance>`

Default
"""""""

none
