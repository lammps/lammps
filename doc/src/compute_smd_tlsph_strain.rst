.. index:: compute smd/tlsph/strain

compute smd/tlsph/strain command
================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/tlsph/strain

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/strain = style name of this compute command

Examples
""""""""

.. parsed-literal::

   compute 1 all smd/tlsph/strain

Description
"""""""""""

Define a computation that calculates the Green-Lagrange strain tensor
for particles interacting via the Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

**Output info:**

This compute calculates a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The per-particle tensor values will be given dimensionless. See
:doc:`units <units>`.

The per-particle vector has 6 entries, corresponding to the xx, yy,
zz, xy, xz, yz components of the symmetric strain tensor.

Restrictions
""""""""""""

This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This quantity will be computed only for particles which interact with
the Total-Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`smd/tlsph/strain/rate <compute_smd_tlsph_strain_rate>`,
:doc:`smd/tlsph/stress <compute_smd_tlsph_stress>`

**Default:** none
