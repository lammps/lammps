.. index:: compute msd

compute msd command
===================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID msd keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* msd = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *com* or *average*

  .. parsed-literal::

       *com* value = *yes* or *no*
       *average* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all msd
   compute 1 upper msd com yes average yes

Description
"""""""""""

Define a computation that calculates the mean-squared displacement
(MSD) of the group of atoms, including all effects due to atoms
passing through periodic boundaries.  For computation of the non-Gaussian
parameter of mean-squared displacement, see the :doc:`compute msd/nongauss <compute_msd_nongauss>` command.

A vector of four quantities is calculated by this compute.  The first 3
elements of the vector are the squared dx,dy,dz displacements, summed
and averaged over atoms in the group.  The 4th element is the total
squared displacement, i.e. (dx\*dx + dy\*dy + dz\*dz), summed and
averaged over atoms in the group.

The slope of the mean-squared displacement (MSD) versus time is
proportional to the diffusion coefficient of the diffusing atoms.

The displacement of an atom is from its reference position. This is
normally the original position at the time
the compute command was issued, unless the *average* keyword is set to *yes*\ .
The value of the displacement will be
0.0 for atoms not in the specified compute group.

If the *com* option is set to *yes* then the effect of any drift in
the center-of-mass of the group of atoms is subtracted out before the
displacement of each atom is calculated.

If the *average* option is set to *yes* then the reference position of
an atom is based on the average position of that atom, corrected for
center-of-mass motion if requested.  The average position is a running
average over all previous calls to the compute, including the current
call. So on the first call it is current position, on the second call
it is the arithmetic average of the current position and the position
on the first call, and so on.  Note that when using this option, the
precise value of the mean square displacement will depend on the
number of times the compute is called. So, for example, changing the
frequency of thermo output may change the computed displacement. Also,
the precise values will be changed if a single simulation is broken up
into two parts, using either multiple run commands or a restart
file. It only makes sense to use this option if the atoms are not
diffusing, so that their average positions relative to the center of
mass of the system are stationary. The most common case is crystalline
solids undergoing thermal motion.

.. note::

   Initial coordinates are stored in "unwrapped" form, by using the
   image flags associated with each atom.  See the :doc:`dump custom <dump>` command for a discussion of "unwrapped" coordinates.
   See the Atoms section of the :doc:`read_data <read_data>` command for a
   discussion of image flags and how they are set for each atom.  You can
   reset the image flags (e.g. to 0) before invoking this compute by
   using the :doc:`set image <set>` command.

.. note::

   If you want the quantities calculated by this compute to be
   continuous when running from a :doc:`restart file <read_restart>`, then
   you should use the same ID for this compute, as in the original run.
   This is so that the fix this compute creates to store per-atom
   quantities will also have the same ID, and thus be initialized
   correctly with atom reference positions from the restart file.  When
   *average* is set to yes, then the atom reference positions are
   restored correctly, but not the number of samples used obtain them. As
   a result, the reference positions from the restart file are combined
   with subsequent positions as if they were from a single sample,
   instead of many, which will change the values of msd somewhat.

**Output info:**

This compute calculates a global vector of length 4, which can be
accessed by indices 1-4 by any command that uses global vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The vector values are "intensive".  The vector values will be in
distance\^2 :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute msd/nongauss <compute_msd_nongauss>`, :doc:`compute displace\_atom <compute_displace_atom>`, :doc:`fix store/state <fix_store_state>`, :doc:`compute msd/chunk <compute_msd_chunk>`

Default
"""""""

The option default are com = no, average = no.
