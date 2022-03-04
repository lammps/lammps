.. index:: compute smd/tlsph/shape

compute smd/tlsph/shape command
===============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/tlsph/shape

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/shape = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/tlsph/shape

Description
"""""""""""

Define a computation that outputs the current shape of the volume
associated with a particle as a rotated ellipsoid.  It is only
meaningful for particles which interact according to the
Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to use Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector of vectors, which can be
accessed by any command that uses per-particle values from a compute
as input. See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

The per-particle vector has 7 entries. The first three entries
correspond to the lengths of the ellipsoid's axes and have units of
length.  These axis values are computed as the contact radius times the
xx, yy, or zz components of the Green-Lagrange strain tensor
associated with the particle.  The next 4 values are quaternions
(order: q, x, y, z) which describe the spatial rotation of the
particle relative to its initial state.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

This quantity will be computed only for particles which interact with
the Total-Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`smd/contact/radius <compute_smd_contact_radius>`

Default
"""""""

none
