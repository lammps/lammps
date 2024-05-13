.. index:: compute smd/tlsph/stress

compute smd/tlsph/stress command
================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID smd/tlsph/stress

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/tlsph/stress = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/tlsph/stress

Description
"""""""""""

Define a computation that outputs the Cauchy stress tensor for
particles interacting via the Total-Lagrangian SPH pair style.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector of vectors (tensors),
which can be accessed by any command that uses per-particle values
from a compute as input. See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The values will be given in :doc:`units <units>` of pressure.

The per-particle vector has 7 entries. The first six entries
correspond to the xx, yy, zz, xy, xz and yz components of the
symmetric Cauchy stress tensor. The seventh entry is the second
invariant of the stress tensor, i.e., the von Mises equivalent stress.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This quantity will be computed only for particles which interact with
the Total-Lagrangian SPH pair style.

Related commands
""""""""""""""""

:doc:`compute smd/tlsph/strain <compute_smd_tlsph_strain>`, :doc:`cmopute smd/tlsph/strain/rate <compute_smd_tlsph_strain_rate>`

Default
"""""""

none
