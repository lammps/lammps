.. index:: compute gyration

compute gyration command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID gyration

* ID, group-ID are documented in :doc:`compute <compute>` command
* gyration = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 molecule gyration

Description
"""""""""""

Define a computation that calculates the radius of gyration :math:`R_g` of the
group of atoms, including all effects due to atoms passing through
periodic boundaries.

:math:`R_g` is a measure of the size of the group of atoms, and is computed as
the square root of the :math:`R_g^2` value in this formula

.. math::

   R_g^2 = \frac{1}{M} \sum_i m_i (r_i - r_{\text{cm}})^2

where :math:`M` is the total mass of the group, :math:`r_{\text{cm}}` is the
center-of-mass position of the group, and the sum is over all atoms in
the group.

A :math:`R_g^2` tensor, stored as a 6-element vector, is also calculated
by this compute.  The formula for the components of the tensor is the
same as the above formula, except that :math:`(r_i - r_{\text{cm}})^2` is
replaced by
:math:`(r_{i,x} - r_{\text{cm},x}) \cdot (r_{i,y} - r_{\text{cm},y})` for the
:math:`xy` component, and so on.  The six components of the vector are ordered
:math:`xx`, :math:`yy`, :math:`zz`, :math:`xy`, :math:`xz`, :math:`yz`.
Note that unlike the scalar :math:`R_g`, each of the six values of the tensor
is effectively a "squared" value, since the cross-terms may be negative
and taking a square root would be invalid.

.. note::

   The coordinates of an atom contribute to :math:`R_g` in "unwrapped" form,
   by using the image flags associated with each atom.  See the
   :doc:`dump custom <dump>` command for a discussion of "unwrapped" coordinates.
   See the Atoms section of the :doc:`read_data <read_data>` command for a
   discussion of image flags and how they are set for each atom.  You can
   reset the image flags (e.g., to 0) before invoking this compute by
   using the :doc:`set image <set>` command.

Output info
"""""""""""

This compute calculates a global scalar (:math:`R_g`) and a global vector of
length 6 (:math:`R_g^2` tensor), which can be accessed by indices 1--6.  These
values can be used by any command that uses a global scalar value or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The scalar and vector values calculated by this compute are
"intensive".  The scalar and vector values will be in distance and
distance\ :math:`^2` :doc:`units <units>`, respectively.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute gyration/chunk <compute_gyration_chunk>`,
:doc:`compute gyration/shape <compute_gyration_shape>`

Default
"""""""

none
