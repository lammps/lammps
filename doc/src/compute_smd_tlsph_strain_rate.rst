.. index:: compute smd/tlsph/strain/rate

compute smd/tlsph/strain/rate command
=====================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/tlsph/strain/rate

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/strain/rate = style name of this compute command

Examples
""""""""

.. parsed-literal::

   compute 1 all smd/tlsph/strain/rate

Description
"""""""""""

Define a computation that calculates the rate of the strain tensor for
particles interacting via the Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

**Output info:**

This compute calculates a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input. See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The values will be given in :doc:`units <units>` of one over time.

The per-particle vector has 6 entries, corresponding to the xx, yy,
zz, xy, xz, yz components of the symmetric strain rate tensor.

Restrictions
""""""""""""

This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This quantity will be computed only for particles which interact with
Total-Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`compute smd/tlsph/strain <compute_smd_tlsph_strain>`, :doc:`compute smd/tlsph/stress <compute_smd_tlsph_stress>`

**Default:** none
