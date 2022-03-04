.. index:: compute dipole

compute dipole command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID dipole charge-correction

* ID, group-ID are documented in :doc:`compute <compute>` command
* dipole = style name of this compute command
* charge-correction = *mass* or *geometry*, use COM or geometric center for charged chunk correction (optional)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 fluid dipole
   compute dw water dipole geometry

Description
"""""""""""

Define a computation that calculates the dipole vector and total dipole
for a group of atoms.

This compute calculates the x,y,z coordinates of the dipole vector
and the total dipole moment for the atoms in the compute group.
This includes all effects due to atoms passing through periodic boundaries.
For a group with a net charge the resulting dipole is made position independent
by subtracting the position vector of the center of mass or geometric center
times the net charge from the computed dipole vector. Both per-atom charges
and per-atom dipole moments, if present, contribute to the computed dipole.

.. note::

   The coordinates of an atom contribute to the dipole in "unwrapped"
   form, by using the image flags associated with each atom.  See the
   :doc:`dump custom <dump>` command for a discussion of "unwrapped"
   coordinates.  See the Atoms section of the :doc:`read_data
   <read_data>` command for a discussion of image flags and how they are
   set for each atom.  You can reset the image flags (e.g. to 0) before
   invoking this compute by using the :doc:`set image <set>` command.

Output info
"""""""""""

This compute calculations a global scalar containing the magnitude of
the computed dipole moment and a global vector of length 3 with the
dipole vector.  See the :doc:`Howto output <Howto_output>` page for
an overview of LAMMPS output options.

The computed values are "intensive".  The array values will be in
dipole units, i.e. charge units times distance :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute dipole/chunk <compute_dipole_chunk>`

Default
"""""""

Using the center of mass is the default setting for the net charge correction.
