.. index:: compute cna/atom

compute cna/atom command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID cna/atom cutoff

* ID, group-ID are documented in :doc:`compute <compute>` command
* cna/atom = style name of this compute command
* cutoff = cutoff distance for nearest neighbors (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all cna/atom 3.08

Description
"""""""""""

Define a computation that calculates the CNA (Common Neighbor
Analysis) pattern for each atom in the group.  In solid-state systems
the CNA pattern is a useful measure of the local crystal structure
around an atom.  The CNA methodology is described in :ref:`(Faken) <Faken>`
and :ref:`(Tsuzuki) <Tsuzuki1>`.

Currently, there are five kinds of CNA patterns LAMMPS recognizes:

* fcc = 1
* hcp = 2
* bcc = 3
* icosahedral = 4
* unknown = 5

The value of the CNA pattern will be 0 for atoms not in the specified
compute group.  Note that normally a CNA calculation should only be
performed on mono-component systems.

The CNA calculation can be sensitive to the specified cutoff value.
You should insure the appropriate nearest neighbors of an atom are
found within the cutoff distance for the presumed crystal structure
(e.g., 12 nearest neighbor for perfect FCC and HCP crystals, 14 nearest
neighbors for perfect BCC crystals).  These formulas can be used to
obtain a good cutoff distance:

.. math::

  r_{c}^{\mathrm{fcc}} = & \frac{1}{2} \left(\frac{\sqrt{2}}{2} + 1\right) a
    \approx 0.8536 a \\
  r_{c}^{\mathrm{bcc}} = & \frac{1}{2}(\sqrt{2} + 1) a
    \approx 1.207 a \\
  r_{c}^{\mathrm{hcp}} = & \frac{1}{2}\left(1+\sqrt{\frac{4+2x^{2}}{3}}\right) a

where :math:`a` is the lattice constant for the crystal structure concerned
and in the HCP case, :math:`x = (c/a) / 1.633`, where 1.633 is the ideal
:math:`c/a` for HCP crystals.

Also note that since the CNA calculation in LAMMPS uses the neighbors
of an owned atom to find the nearest neighbors of a ghost atom, the
following relation should also be satisfied:

.. math::

  r_c + r_s > 2*{\rm cutoff}

where :math:`r_c` is the cutoff distance of the potential, :math:`r_s`
is the skin
distance as specified by the :doc:`neighbor <neighbor>` command, and
cutoff is the argument used with the compute cna/atom command.  LAMMPS
will issue a warning if this is not the case.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (e.g. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently or to have multiple compute/dump commands, each with a
*cna/atom* style.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values will be a number from 0 to 5, as explained
above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute centro/atom <compute_centro_atom>`

Default
"""""""

none

----------

.. _Faken:

**(Faken)** Faken, Jonsson, Comput Mater Sci, 2, 279 (1994).

.. _Tsuzuki1:

**(Tsuzuki)** Tsuzuki, Branicio, Rino, Comput Phys Comm, 177, 518 (2007).
