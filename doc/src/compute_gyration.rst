.. index:: compute gyration

compute gyration command
========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID gyration

* ID, group-ID are documented in :doc:`compute <compute>` command
* gyration = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 molecule gyration

Description
"""""""""""

Define a computation that calculates the radius of gyration Rg of the
group of atoms, including all effects due to atoms passing through
periodic boundaries.

Rg is a measure of the size of the group of atoms, and is computed as
the square root of the Rg\^2 value in this formula

.. image:: Eqs/compute_gyration.jpg
   :align: center

where M is the total mass of the group, Rcm is the center-of-mass
position of the group, and the sum is over all atoms in the group.

A Rg\^2 tensor, stored as a 6-element vector, is also calculated by
this compute.  The formula for the components of the tensor is the
same as the above formula, except that (Ri - Rcm)\^2 is replaced by
(Rix - Rcmx) \* (Riy - Rcmy) for the xy component, etc.  The 6
components of the vector are ordered xx, yy, zz, xy, xz, yz.  Note
that unlike the scalar Rg, each of the 6 values of the tensor is
effectively a "squared" value, since the cross-terms may be negative
and taking a sqrt() would be invalid.

.. note::

   The coordinates of an atom contribute to Rg in "unwrapped" form,
   by using the image flags associated with each atom.  See the :doc:`dump custom <dump>` command for a discussion of "unwrapped" coordinates.
   See the Atoms section of the :doc:`read_data <read_data>` command for a
   discussion of image flags and how they are set for each atom.  You can
   reset the image flags (e.g. to 0) before invoking this compute by
   using the :doc:`set image <set>` command.

**Output info:**

This compute calculates a global scalar (Rg) and a global vector of
length 6 (Rg\^2 tensor), which can be accessed by indices 1-6.  These
values can be used by any command that uses a global scalar value or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The scalar and vector values calculated by this compute are
"intensive".  The scalar and vector values will be in distance and
distance\^2 :doc:`units <units>` respectively.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute gyration/chunk <compute_gyration_chunk>`,
:doc:`compute gyration/shape <compute_gyration_shape>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
