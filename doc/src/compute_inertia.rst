.. index:: compute inertia

compute inertia command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID inertia

* ID, group-ID are documented in :doc:`compute <compute>` command
* inertia = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 molecule inertia

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the symmetric moment of inertia
tensor of the group of atoms around its center of mass, including all
effects due to atoms passing through periodic boundaries.

The inertia tensor is computed as

.. math::

   I = \sum_i \left[ m_i \left( (r_i - r_{\text{cm}})^2 \mathbf{1}
       - (r_i - r_{\text{cm}}) \otimes (r_i - r_{\text{cm}}) \right)
       + I_i \right]

where :math:`r_{\text{cm}}` is the center-of-mass position of the group,
the sum is over all atoms in the group, and :math:`I_i` is the moment of
inertia of atom :math:`i` about its own center.  For point particles
:math:`I_i` is zero; for :doc:`finite-size particles <Howto_spherical>`
(finite-size spheres, ellipsoids, superellipsoids, line segments,
triangles, and body particles) :math:`I_i` is the moment of inertia of
the particle's shape, which is added using the parallel-axis theorem.
For finite spheres, ellipsoids, line segments, and triangles this
matches the treatment used by the :doc:`fix rigid <fix_rigid>` command;
for superellipsoids and body particles the per-particle moment of
inertia stored with each particle is used (these are not supported by
:doc:`fix rigid <fix_rigid>`).

.. note::

   The coordinates of an atom contribute to the inertia tensor in
   "unwrapped" form, by using the image flags associated with each atom.
   See the :doc:`dump custom <dump>` command for a discussion of
   "unwrapped" coordinates.  See the Atoms section of the
   :doc:`read_data <read_data>` command for a discussion of image flags and
   how they are set for each atom.  You can reset the image flags (e.g., to
   0) before invoking this compute by using the :doc:`set image <set>`
   command.

Output info
"""""""""""

This compute calculates a global vector of length 6, which can be accessed
by indices 1--6.  The six components of the symmetric inertia tensor are
ordered :math:`I_{xx}, I_{yy}, I_{zz}, I_{xy}, I_{yz}, I_{xz}`.  These values
can be used by any command that uses global vector values from a compute as
input.  See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The vector values calculated by this compute are "intensive".  The vector
values will be in mass\*distance\ :math:`^2` :doc:`units <units>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`compute inertia/chunk <compute_inertia_chunk>`,
:doc:`compute gyration <compute_gyration>`,
:doc:`variable inertia() function <variable>`

Default
"""""""

none
