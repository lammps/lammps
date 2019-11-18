.. index:: atom_modify

atom_modify command
===================

Syntax
""""""


.. code-block:: LAMMPS

   atom_modify keyword values ...

* one or more keyword/value pairs may be appended
* keyword = *id* or *map* or *first* or *sort*
  
  .. parsed-literal::
  
        *id* value = *yes* or *no*
        *map* value = *yes* or *array* or *hash*
        *first* value = group-ID = group whose atoms will appear first in internal atom lists
        *sort* values = Nfreq binsize
          Nfreq = sort atoms spatially every this many time steps
          binsize = bin size for spatial sorting (distance units)



Examples
""""""""


.. code-block:: LAMMPS

   atom_modify map yes
   atom_modify map hash sort 10000 2.0
   atom_modify first colloid

Description
"""""""""""

Modify certain attributes of atoms defined and stored within LAMMPS,
in addition to what is specified by the :doc:`atom\_style <atom_style>`
command.  The *id* and *map* keywords must be specified before a
simulation box is defined; other keywords can be specified any time.

The *id* keyword determines whether non-zero atom IDs can be assigned
to each atom.  If the value is *yes*\ , which is the default, IDs are
assigned, whether you use the :doc:`create atoms <create_atoms>` or
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands to initialize atoms.  If the value is *no* the IDs for all
atoms are assumed to be 0.

If atom IDs are used, they must all be positive integers.  They should
also be unique, though LAMMPS does not check for this.  Typically they
should also be consecutively numbered (from 1 to Natoms), though this
is not required.  Molecular :doc:`atom styles <atom_style>` are those
that store bond topology information (styles bond, angle, molecular,
full).  These styles require atom IDs since the IDs are used to encode
the topology.  Some other LAMMPS commands also require the use of atom
IDs.  E.g. some many-body pair styles use them to avoid double
computation of the I-J interaction between two atoms.

The only reason not to use atom IDs is if you are running an atomic
simulation so large that IDs cannot be uniquely assigned.  For a
default LAMMPS build this limit is 2\^31 or about 2 billion atoms.
However, even in this case, you can use 64-bit atom IDs, allowing 2\^63
or about 9e18 atoms, if you build LAMMPS with the - DLAMMPS\_BIGBIG
switch.  This is described on the :doc:`Build\_settings <Build_settings>`
doc page.  If atom IDs are not used, they must be specified as 0 for
all atoms, e.g. in a data or restart file.

The *map* keyword determines how atoms with specific IDs are found
when required.  An example are the bond (angle, etc) methods which
need to find the local index of an atom with a specific global ID
which is a bond (angle, etc) partner.  LAMMPS performs this operation
efficiently by creating a "map", which is either an *array* or *hash*
table, as described below.

When the *map* keyword is not specified in your input script, LAMMPS
only creates a map for :doc:`atom\_styles <atom_style>` for molecular
systems which have permanent bonds (angles, etc).  No map is created
for atomic systems, since it is normally not needed.  However some
LAMMPS commands require a map, even for atomic systems, and will
generate an error if one does not exist.  The *map* keyword thus
allows you to force the creation of a map.  The *yes* value will
create either an *array* or *hash* style map, as explained in the next
paragraph.  The *array* and *hash* values create an atom-style or
hash-style map respectively.

For an *array*\ -style map, each processor stores a lookup table of
length N, where N is the largest atom ID in the system.  This is a
fast, simple method for many simulations, but requires too much memory
for large simulations.  For a *hash*\ -style map, a hash table is
created on each processor, which finds an atom ID in constant time
(independent of the global number of atom IDs).  It can be slightly
slower than the *array* map, but its memory cost is proportional to
the number of atoms owned by a processor, i.e. N/P when N is the total
number of atoms in the system and P is the number of processors.

The *first* keyword allows a :doc:`group <group>` to be specified whose
atoms will be maintained as the first atoms in each processor's list
of owned atoms.  This in only useful when the specified group is a
small fraction of all the atoms, and there are other operations LAMMPS
is performing that will be sped-up significantly by being able to loop
over the smaller set of atoms.  Otherwise the reordering required by
this option will be a net slow-down.  The :doc:`neigh\_modify include <neigh_modify>` and :doc:`comm\_modify group <comm_modify>`
commands are two examples of commands that require this setting to
work efficiently.  Several :doc:`fixes <fix>`, most notably time
integration fixes like :doc:`fix nve <fix_nve>`, also take advantage of
this setting if the group they operate on is the group specified by
this command.  Note that specifying "all" as the group-ID effectively
turns off the *first* option.

It is OK to use the *first* keyword with a group that has not yet been
defined, e.g. to use the atom\_modify first command at the beginning of
your input script.  LAMMPS does not use the group until a simulation
is run.

The *sort* keyword turns on a spatial sorting or reordering of atoms
within each processor's sub-domain every *Nfreq* timesteps.  If
*Nfreq* is set to 0, then sorting is turned off.  Sorting can improve
cache performance and thus speed-up a LAMMPS simulation, as discussed
in a paper by :ref:`(Meloni) <Meloni>`.  Its efficacy depends on the problem
size (atoms/processor), how quickly the system becomes disordered, and
various other factors.  As a general rule, sorting is typically more
effective at speeding up simulations of liquids as opposed to solids.
In tests we have done, the speed-up can range from zero to 3-4x.

Reordering is performed every *Nfreq* timesteps during a dynamics run
or iterations during a minimization.  More precisely, reordering
occurs at the first reneighboring that occurs after the target
timestep.  The reordering is performed locally by each processor,
using bins of the specified *binsize*\ .  If *binsize* is set to 0.0,
then a binsize equal to half the :doc:`neighbor <neighbor>` cutoff
distance (force cutoff plus skin distance) is used, which is a
reasonable value.  After the atoms have been binned, they are
reordered so that atoms in the same bin are adjacent to each other in
the processor's 1d list of atoms.

The goal of this procedure is for atoms to put atoms close to each
other in the processor's one-dimensional list of atoms that are also
near to each other spatially.  This can improve cache performance when
pairwise interactions and neighbor lists are computed.  Note that if
bins are too small, there will be few atoms/bin.  Likewise if bins are
too large, there will be many atoms/bin.  In both cases, the goal of
cache locality will be undermined.

.. note::

   Running a simulation with sorting on versus off should not
   change the simulation results in a statistical sense.  However, a
   different ordering will induce round-off differences, which will lead
   to diverging trajectories over time when comparing two simulations.
   Various commands, particularly those which use random numbers
   (e.g. :doc:`velocity create <velocity>`, and :doc:`fix langevin <fix_langevin>`), may generate (statistically identical)
   results which depend on the order in which atoms are processed.  The
   order of atoms in a :doc:`dump <dump>` file will also typically change
   if sorting is enabled.

Restrictions
""""""""""""


The *first* and *sort* options cannot be used together.  Since sorting
is on by default, it will be turned off if the *first* keyword is
used with a group-ID that is not "all".

**Related commands:** none

Default
"""""""

By default, *id* is yes.  By default, atomic systems (no bond topology
info) do not use a map.  For molecular systems (with bond topology
info), a map is used.  The default map style is array if no atom ID is
larger than 1 million, otherwise the default is hash.  By default, a
"first" group is not defined.  By default, sorting is enabled with a
frequency of 1000 and a binsize of 0.0, which means the neighbor
cutoff will be used to set the bin size. If no neighbor cutoff is
defined, sorting will be turned off.


----------


.. _Meloni:



**(Meloni)** Meloni, Rosati and Colombo, J Chem Phys, 126, 121102 (2007).
