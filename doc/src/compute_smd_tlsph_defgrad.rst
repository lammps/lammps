.. index:: compute smd/tlsph/defgrad

compute smd/tlsph/defgrad command
=================================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID smd/tlsph/defgrad

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/defgrad = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all smd/tlsph/defgrad

Description
"""""""""""

Define a computation that calculates the deformation gradient.  It is
only meaningful for particles which interact according to the
Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to use Smooth
Mach Dynamics in LAMMPS.

**Output info:**

This compute outputs a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input. See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The per-particle vector values will be given dimensionless. See
:doc:`units <units>`.  The per-particle vector has 10 entries. The first
nine entries correspond to the xx, xy, xz, yx, yy, yz, zx, zy, zz
components of the asymmetric deformation gradient tensor. The tenth
entry is the determinant of the deformation gradient.

Restrictions
""""""""""""


This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info. TThis compute can
only be used for particles which interact via the total Lagrangian SPH
pair style.

Related commands
""""""""""""""""

:doc:`smd/hourglass/error <compute_smd_hourglass_error>`

**Default:** none
