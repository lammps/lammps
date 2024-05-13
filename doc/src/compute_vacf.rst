.. index:: compute vacf

compute vacf command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID vacf

* ID, group-ID are documented in :doc:`compute <compute>` command
* vacf = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all vacf
   compute 1 upper vacf

Description
"""""""""""

Define a computation that calculates the velocity auto-correlation
function (VACF), averaged over a group of atoms.  Each atom's
contribution to the VACF is its current velocity vector dotted into
its initial velocity vector at the time the compute was specified.

A vector of four quantities is calculated by this compute.  The first three
elements of the vector are :math:`v_x v_{x,0}` (and similar for the
:math:`y` and :math:`z` components), summed and averaged over atoms in the
group, where :math:`v_x` is the current :math:`x`-component of the velocity of
the atom and :math:`v_{x,0}` is the initial :math:`x`-component of the velocity
of the atom.  The fourth element of the vector is the total VACF
(i.e., :math:`(v_x v_{x,0} + v_y v_{y,0} + v_z v_{z,0})`),
summed and averaged over atoms in the group.

The integral of the VACF versus time is proportional to the diffusion
coefficient of the diffusing atoms.  This can be computed in the
following manner, using the :doc:`variable trap() <variable>` function:

.. code-block:: LAMMPS

   compute         2 all vacf
   fix             5 all vector 1 c_2[4]
   variable        diff equal dt*trap(f_5)
   thermo_style    custom step v_diff

.. note::

   If you want the quantities calculated by this compute to be
   continuous when running from a :doc:`restart file <read_restart>`, then
   you should use the same ID for this compute, as in the original run.
   This is so that the fix this compute creates to store per-atom
   quantities will also have the same ID, and thus be initialized
   correctly with time=0 atom velocities from the restart file.

Output info
"""""""""""

This compute calculates a global vector of length 4, which can be
accessed by indices 1--4 by any command that uses global vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The vector values are "intensive".  The vector values will be in
velocity\ :math:`^2` :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`

Default
"""""""

none
