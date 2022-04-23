.. index:: compute smd/ulsph/strain/rate

compute smd/ulsph/strain/rate command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID smd/ulsph/strain/rate

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/ulsph/strain/rate = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/ulsph/strain/rate

Description
"""""""""""

Define a computation that outputs the rate of the logarithmic strain
tensor for particles interacting via the updated Lagrangian SPH pair
style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input. See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The values will be given in :doc:`units <units>` of one over time.

The per-particle vector has 6 entries, corresponding to the xx, yy,
zz, xy, xz, yz components of the symmetric strain rate tensor.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This compute can only be used for particles which interact with the
updated Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`compute smd/tlsph/strain/rate <compute_smd_tlsph_strain_rate>`

Default
"""""""

none
