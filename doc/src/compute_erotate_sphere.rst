.. index:: compute erotate/sphere
.. index:: compute erotate/sphere/kk

compute erotate/sphere command
==============================

Accelerator Variants: *erotate/sphere/kk*


Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID erotate/sphere

* ID, group-ID are documented in :doc:`compute <compute>` command
* erotate/sphere = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all erotate/sphere

Description
"""""""""""

Define a computation that calculates the rotational kinetic energy of
a group of spherical particles.

The rotational energy is computed as :math:`\frac12 I \omega^2`,
where :math:`I` is the moment of inertia for a sphere and :math:`\omega`
is the particle's angular velocity.

.. note::

   For :doc:`2d models <dimension>`, particles are treated as
   spheres, not disks, meaning their moment of inertia will be the same
   as in 3d.

----------

.. include:: accel_styles.rst

----------

Output info
"""""""""""

This compute calculates a global scalar (the KE).  This value can be
used by any command that uses a global scalar value from a compute as
input.  See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".  The
scalar value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute requires that atoms store a radius and angular velocity
(omega) as defined by the :doc:`atom_style sphere <atom_style>` command.

All particles in the group must be finite-size spheres or point
particles.  They cannot be aspherical.  Point particles will not
contribute to the rotational energy.

Related commands
""""""""""""""""

:doc:`compute erotate/asphere <compute_erotate_asphere>`

Default
"""""""

none
