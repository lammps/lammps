.. index:: compute smd/ulsph/stress

compute smd/ulsph/stress command
================================

Syntax
""""""


.. code-block:: LAMMPS

   compute ID group-ID smd/ulsph/stress

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/ulsph/stress = style name of this compute command

Examples
""""""""


.. code-block:: LAMMPS

   compute 1 all smd/ulsph/stress

Description
"""""""""""

Define a computation that outputs the Cauchy stress tensor.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

**Output info:**

This compute calculates a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input. See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The values will be given in :doc:`units <units>` of pressure.

The per-particle vector has 7 entries. The first six entries
correspond to the xx, yy, zz, xy, xz, yz components of the symmetric
Cauchy stress tensor. The seventh entry is the second invariant of the
stress tensor, i.e., the von Mises equivalent stress.

Restrictions
""""""""""""


This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info. This compute can
only be used for particles which interact with the updated Lagrangian
SPH pair style.

Related commands
""""""""""""""""

:doc:`compute smd/ulsph/strain <compute_smd_ulsph_strain>`, :doc:`compute smd/ulsph/strain/rate <compute_smd_ulsph_strain_rate>` :doc:`compute smd/tlsph/stress <compute_smd_tlsph_stress>`

**Default:** none
