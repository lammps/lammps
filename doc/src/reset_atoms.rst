.. index:: reset_atoms

reset_atoms command
===================

Syntax
""""""

.. code-block:: LAMMPS

   reset_atoms property arguments ...

* property = *id* or *image* or *mol*
* additional arguments depend on the property

  .. code-block:: LAMMPS

     reset_atoms id keyword value ...

  * zero or more keyword/value pairs can be appended
  * keyword = *sort*

    .. parsed-literal::

       *sort* value = *yes* or *no*

  .. code-block:: LAMMPS

     reset_atoms image group-ID

  * group-ID = ID of group of atoms whose image flags will be reset

  .. code-block:: LAMMPS

     reset atoms mol group-ID keyword value ...

  * group-ID = ID of group of atoms whose molecule IDs will be reset
  * zero or more keyword/value pairs can be appended
  * keyword = *compress* or *offset* or *single*

    .. parsed-literal::

       *compress* value = *yes* or *no*
       *offset* value = *Noffset* >= -1
       *single* value = *yes* or *no* to treat single atoms (no bonds) as molecules


Examples
""""""""

.. code-block:: LAMMPS

   reset_atoms id
   reset_atoms id sort yes
   reset_atoms image all
   reset_atoms image mobile
   reset_atoms mol all
   reset_atoms mol all offset 10 single yes
   reset_atoms mol solvent compress yes offset 100
   reset_atoms mol solvent compress no


Description
"""""""""""

.. versionadded:: 22Dec2022

The *reset_atoms* command resets the values of a specified atom
property.  In contrast to the set command, it does this in a
collective manner which resets the values for many atoms in a
self-consistent way.  This is often useful when the simulated system
has undergone significant modifications like adding or removing atoms
or molecules, joining data files, changing bonds, or large-scale
diffusion.

The new values can be thought of as a *reset*, similar to values atoms
would have if a new data file were being read or a new simulation
performed.  Note that the set command also resets atom properties to
new values, but it treats each atom independently.

The *property* setting can be *id* or *image* or *mol*.  For *id*, the
IDs of all the atoms are reset to contiguous values.  For *image*, the
image flags of atoms in the specified *group-ID* are reset so that at
least one atom in each molecule is in the simulation box (image flag =
0).  For *mol*, the molecule IDs of all atoms are reset to contiguous
values.

More details on these operations and their arguments or optional
keyword/value settings are given below.

----------

*Property id*

Reset atom IDs for the entire system, including all the global IDs
stored for bond, angle, dihedral, improper topology data.  This will
create a set of IDs that are numbered contiguously from 1 to N for a N
atoms system.

This can be useful to do after performing a "delete_atoms" command for
a molecular system.  The delete_atoms compress yes option will not
perform this operation due to the existence of bond topology.  It can
also be useful to do after any simulation which has lost atoms,
e.g. due to atoms moving outside a simulation box with fixed
boundaries (see the "boundary command"), or due to evaporation (see
the "fix evaporate" command).

If the *sort* keyword is used with a setting of *yes*, then the
assignment of new atom IDs will be the same no matter how many
processors LAMMPS is running on.  This is done by first doing a
spatial sort of all the atoms into bins and sorting them within each
bin.  Because the set of bins is independent of the number of
processors, this enables a consistent assignment of new IDs to each
atom.

This can be useful to do after using the "create_atoms" command and/or
"replicate" command.  In general those commands do not guarantee
assignment of the same atom ID to the same physical atom when LAMMPS
is run on different numbers of processors.  Enforcing consistent IDs
can be useful for debugging or comparing output from two different
runs.

Note that the spatial sort requires communication of atom IDs and
coordinates between processors in an all-to-all manner.  This is done
efficiently in LAMMPS, but it is more expensive than how atom IDs are
reset without sorting.

Note that whether sorting or not, the resetting of IDs is not a
compression, where gaps in atom IDs are removed by decrementing atom
IDs that are larger.  Instead the IDs for all atoms are erased, and
new IDs are assigned so that the atoms owned by an individual
processor have consecutive IDs, as the :doc:`create_atoms
<create_atoms>` command explains.

.. note::

   If this command is used before a :doc:`pair style <pair_style>` is
   defined, an error about bond topology atom IDs not being found may
   result.  This is because the cutoff distance for ghost atom
   communication was not sufficient to find atoms in bonds, angles, etc
   that are owned by other processors.  The :doc:`comm_modify cutoff
   <comm_modify>` command can be used to correct this issue.  Or you can
   define a pair style before using this command.  If you do the former,
   you should unset the *comm_modify cutoff* after using *reset
   atoms id* so that subsequent communication is not inefficient.

----------

*Property image*

Reset the image flags of atoms so that at least one atom in each
molecule has an image flag of 0.  Molecular topology is respected so
that if the molecule straddles a periodic simulation box boundary, the
images flags of all atoms in the molecule will be consistent.  This
avoids inconsistent image flags that could result from resetting all
image flags to zero with the :doc:`set <set>` command.

.. note::

   If the system has no bonds, there is no reason to use this command,
   since image flags for different atoms do not need to be
   consistent.  Use the :doc:`set <set>` command with its *image*
   keyword instead.

Only image flags for atoms in the specified *group-ID* are reset; all
others remain unchanged.  No check is made for whether the group
covers complete molecule fragments and thus whether the command will
result in inconsistent image flags.

Molecular fragments are identified by the algorithm used by the
:doc:`compute fragment/atom <compute_cluster_atom>` command.  For each
fragment the average of the largest and the smallest image flag in
each direction across all atoms in the fragment is computed and
subtracted from the current image flag in the same direction.

This can be a useful operation to perform after running longer
equilibration runs of mobile systems where molecules would pass
through the system multiple times and thus produce non-zero image
flags.

.. note::

   Same as explained for the :doc:`compute fragment/atom
   <compute_cluster_atom>` command, molecules are identified using the
   current bond topology.  This will **not** account for bonds broken by
   the :doc:`bond_style quartic <bond_quartic>` command, because this
   bond style does not perform a full update of the bond topology data
   structures within LAMMPS.  In that case, using the :doc:`delete_bonds
   all bond 0 remove <delete_bonds>` will permanently delete such
   broken bonds and should thus be used first.

----------

*Property mol*

Reset molecule IDs for a specified group of atoms based on current
bond connectivity.  This will typically create a new set of molecule
IDs for atoms in the group.  Only molecule IDs for atoms in the
specified *group-ID* are reset; molecule IDs for atoms not in the
group are not changed.

For purposes of this operation, molecules are identified by the current
bond connectivity in the system, which may or may not be consistent with
the current molecule IDs.  A molecule in this context is a set of atoms
connected to each other with explicit bonds.  The specific algorithm
used is the one of :doc:`compute fragment/atom <compute_cluster_atom>`
Once the molecules are identified and a new molecule ID computed for
each, this command will update the current molecule ID for all atoms in
the group with the new molecule ID.  Note that if the group excludes
atoms within molecules, one (physical) molecule may become two or more
(logical) molecules.  For example if the group excludes atoms in the
middle of a linear chain, then each end of the chain is considered an
independent molecule and will be assigned a different molecule ID.

This can be a useful operation to perform after running reactive
molecular dynamics run with :doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`, or :doc:`fix bond/break
<fix_bond_break>`, all of which can change molecule topologies. It can
also be useful after molecules have been deleted with the
:doc:`delete_atoms <delete_atoms>` command or after a simulation which
has lost molecules, e.g. via the :doc:`fix evaporate <fix_evaporate>`
command.

The *compress* keyword determines how new molecule IDs are computed.  If
the setting is *yes* (the default) and there are N molecules in the
group, the new molecule IDs will be a set of N contiguous values.  See
the *offset* keyword for details on selecting the range of these values.
If the setting is *no*, the molecule ID of every atom in the molecule
will be set to the smallest atom ID of any atom in the molecule.

The *single* keyword determines whether single atoms (not bonded to
another atom) are treated as one-atom molecules or not, based on the
*yes* or *no* setting.  If the setting is *no* (the default), their
molecule IDs are set to 0.  This setting can be important if the new
molecule IDs will be used as input to other commands such as
:doc:`compute chunk/atom molecule <compute_chunk_atom>` or :doc:`fix
rigid molecule <fix_rigid>`.

The *offset* keyword is only used if the *compress* setting is *yes*.
Its default value is *Noffset* = -1.  In that case, if the specified
group is *all*, then the new compressed molecule IDs will range from 1
to N.  If the specified group is not *all* and the largest molecule ID
of atoms outside that group is M, then the new compressed molecule IDs will
range from M+1 to M+N, to avoid collision with existing molecule
IDs.  If an *Noffset* >= 0 is specified, then the new compressed
molecule IDs will range from *Noffset*\ +1 to *Noffset*\ +N.  If the group
is not *all* there may be collisions with the molecule IDs of other atoms.

.. note::

   Same as explained for the :doc:`compute fragment/atom
   <compute_cluster_atom>` command, molecules are identified using the
   current bond topology.  This will **not** account for bonds broken by
   the :doc:`bond_style quartic <bond_quartic>` command, because this
   bond style does not perform a full update of the bond topology data
   structures within LAMMPS.  In that case, using the :doc:`delete_bonds
   all bond 0 remove <delete_bonds>` will permanently delete such broken
   bonds and should thus be used first.


Restrictions
""""""""""""

The *image* property can only be used when the atom style supports bonds.

Related commands
""""""""""""""""

:doc:`compute fragment/atom <compute_cluster_atom>`
:doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`,
:doc:`fix bond/break <fix_bond_break>`,
:doc:`fix evaporate <fix_evaporate>`,
:doc:`delete_atoms <delete_atoms>`,
:doc:`delete_bonds <delete_bonds>`

Defaults
""""""""

For property *id*, the default keyword setting is sort = no.

For property *mol*, the default keyword settings are compress = yes,
single = no, and offset = -1.
