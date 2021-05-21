.. index:: compute com

compute com command
===================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID com

* ID, group-ID are documented in :doc:`compute <compute>` command
* com = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all com

Description
"""""""""""

Define a computation that calculates the center-of-mass of the group
of atoms, including all effects due to atoms passing through periodic
boundaries.

A vector of three quantities is calculated by this compute, which
are the x,y,z coordinates of the center of mass.

.. note::

   The coordinates of an atom contribute to the center-of-mass in
   "unwrapped" form, by using the image flags associated with each atom.
   See the :doc:`dump custom <dump>` command for a discussion of
   "unwrapped" coordinates.  See the Atoms section of the
   :doc:`read_data <read_data>` command for a discussion of image flags and
   how they are set for each atom.  You can reset the image flags
   (e.g. to 0) before invoking this compute by using the :doc:`set image <set>` command.

Output info
"""""""""""

This compute calculates a global vector of length 3, which can be
accessed by indices 1-3 by any command that uses global vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The vector values are "intensive".  The vector values will be in
distance :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute com/chunk <compute_com_chunk>`

Default
"""""""

none
