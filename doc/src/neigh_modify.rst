.. index:: neigh\_modify

neigh\_modify command
=====================

Syntax
""""""


.. parsed-literal::

   neigh_modify keyword values ...

* one or more keyword/value pairs may be listed
  
  .. parsed-literal::
  
     keyword = *delay* or *every* or *check* or *once* or *cluster* or *include* or *exclude* or *page* or *one* or *binsize*
       *delay* value = N
         N = delay building until this many steps since last build
       *every* value = M
         M = build neighbor list every this many steps
       *check* value = *yes* or *no*
         *yes* = only build if some atom has moved half the skin distance or more
         *no* = always build on 1st step that *every* and *delay* are satisfied
       *once*
         *yes* = only build neighbor list once at start of run and never rebuild
         *no* = rebuild neighbor list according to other settings
       *cluster*
         *yes* = check bond,angle,etc neighbor list for nearby clusters
         *no* = do not check bond,angle,etc neighbor list for nearby clusters
       *include* value = group-ID
         group-ID = only build pair neighbor lists for atoms in this group
       *exclude* values:
         type M N
           M,N = exclude if one atom in pair is type M, other is type N
         group group1-ID group2-ID
           group1-ID,group2-ID = exclude if one atom is in 1st group, other in 2nd
         molecule/intra group-ID
           group-ID = exclude if both atoms are in the same molecule and in group
         molecule/inter group-ID
           group-ID = exclude if both atoms are in different molecules and in group
         none
           delete all exclude settings
       *page* value = N
         N = number of pairs stored in a single neighbor page
       *one* value = N
         N = max number of neighbors of one atom
       *binsize* value = size
         size = bin size for neighbor list construction (distance units)



Examples
""""""""


.. parsed-literal::

   neigh_modify every 2 delay 10 check yes page 100000
   neigh_modify exclude type 2 3
   neigh_modify exclude group frozen frozen check no
   neigh_modify exclude group residue1 chain3
   neigh_modify exclude molecule/intra rigid

Description
"""""""""""

This command sets parameters that affect the building and use of
pairwise neighbor lists.  Depending on what pair interactions and
other commands are defined, a simulation may require one or more
neighbor lists.

The *every*\ , *delay*\ , *check*\ , and *once* options affect how often
lists are built as a simulation runs.  The *delay* setting means never
build new lists until at least N steps after the previous build.  The
*every* setting means build lists every M steps (after the delay has
passed).  If the *check* setting is *no*\ , the lists are built on the
first step that satisfies the *delay* and *every* settings.  If the
*check* setting is *yes*\ , then the *every* and *delay* settings
determine when a build may possibly be performed, but an actual build
only occurs if some atom has moved more than half the skin distance
(specified in the :doc:`neighbor <neighbor>` command) since the last
build.

If the *once* setting is yes, then the neighbor list is only built
once at the beginning of each run, and never rebuilt, except on steps
when a restart file is written, or steps when a fix forces a rebuild
to occur (e.g. fixes that create or delete atoms, such as :doc:`fix deposit <fix_deposit>` or :doc:`fix evaporate <fix_evaporate>`).
This setting should only be made if you are certain atoms will not
move far enough that the neighbor list should be rebuilt, e.g. running
a simulation of a cold crystal.  Note that it is not that expensive to
check if neighbor lists should be rebuilt.

When the rRESPA integrator is used (see the :doc:`run_style <run_style>`
command), the *every* and *delay* parameters refer to the longest
(outermost) timestep.

The *cluster* option does a sanity test every time neighbor lists are
built for bond, angle, dihedral, and improper interactions, to check
that each set of 2, 3, or 4 atoms is a cluster of nearby atoms.  It
does this by computing the distance between pairs of atoms in the
interaction and insuring they are not further apart than half the
periodic box length.  If they are, an error is generated, since the
interaction would be computed between far-away atoms instead of their
nearby periodic images.  The only way this should happen is if the
pairwise cutoff is so short that atoms that are part of the same
interaction are not communicated as ghost atoms.  This is an unusual
model (e.g. no pair interactions at all) and the problem can be fixed
by use of the :doc:`comm_modify cutoff <comm_modify>` command.  Note
that to save time, the default *cluster* setting is *no*\ , so that this
check is not performed.

The *include* option limits the building of pairwise neighbor lists to
atoms in the specified group.  This can be useful for models where a
large portion of the simulation is particles that do not interact with
other particles or with each other via pairwise interactions.  The
group specified with this option must also be specified via the
:doc:`atom_modify first <atom_modify>` command.  Note that specifying
"all" as the group-ID effectively turns off the *include* option.

The *exclude* option turns off pairwise interactions between certain
pairs of atoms, by not including them in the neighbor list.  These are
sample scenarios where this is useful:

* In crack simulations, pairwise interactions can be shut off between 2
  slabs of atoms to effectively create a crack.
* When a large collection of atoms is treated as frozen, interactions
  between those atoms can be turned off to save needless
  computation. E.g. Using the :doc:`fix setforce <fix_setforce>` command
  to freeze a wall or portion of a bio-molecule.
* When one or more rigid bodies are specified, interactions within each
  body can be turned off to save needless computation.  See the :doc:`fix rigid <fix_rigid>` command for more details.


The *exclude type* option turns off the pairwise interaction if one
atom is of type M and the other of type N.  M can equal N.  The
*exclude group* option turns off the interaction if one atom is in the
first group and the other is the second.  Group1-ID can equal
group2-ID.  The *exclude molecule/intra* option turns off the
interaction if both atoms are in the specified group and in the same
molecule, as determined by their molecule ID.  The *exclude
molecule/inter* turns off the interaction between pairs of atoms that
have different molecule IDs and are both in the specified group.

Each of the exclude options can be specified multiple times.  The
*exclude type* option is the most efficient option to use; it requires
only a single check, no matter how many times it has been specified.
The other exclude options are more expensive if specified multiple
times; they require one check for each time they have been specified.

Note that the exclude options only affect pairwise interactions; see
the :doc:`delete_bonds <delete_bonds>` command for information on
turning off bond interactions.

.. note::

   Excluding pairwise interactions will not work correctly when
   also using a long-range solver via the
   :doc:`kspace_style <kspace_style>` command.  LAMMPS will give a warning
   to this effect.  This is because the short-range pairwise interaction
   needs to subtract off a term from the total energy for pairs whose
   short-range interaction is excluded, to compensate for how the
   long-range solver treats the interaction.  This is done correctly for
   pairwise interactions that are excluded (or weighted) via the
   :doc:`special_bonds <special_bonds>` command.  But it is not done for
   interactions that are excluded via these neigh\_modify exclude options.

The *page* and *one* options affect how memory is allocated for the
neighbor lists.  For most simulations the default settings for these
options are fine, but if a very large problem is being run or a very
long cutoff is being used, these parameters can be tuned.  The indices
of neighboring atoms are stored in "pages", which are allocated one
after another as they fill up.  The size of each page is set by the
*page* value.  A new page is allocated when the next atom's neighbors
could potentially overflow the list.  This threshold is set by the
*one* value which tells LAMMPS the maximum number of neighbor's one
atom can have.

.. note::

   LAMMPS can crash without an error message if the number of
   neighbors for a single particle is larger than the *page* setting,
   which means it is much, much larger than the *one* setting.  This is
   because LAMMPS doesn't error check these limits for every pairwise
   interaction (too costly), but only after all the particle's neighbors
   have been found.  This problem usually means something is very wrong
   with the way you've setup your problem (particle spacing, cutoff
   length, neighbor skin distance, etc).  If you really expect that many
   neighbors per particle, then boost the *one* and *page* settings
   accordingly.

The *binsize* option allows you to specify what size of bins will be
used in neighbor list construction to sort and find neighboring atoms.
By default, for :doc:`neighbor style bin <neighbor>`, LAMMPS uses bins
that are 1/2 the size of the maximum pair cutoff.  For :doc:`neighbor style multi <neighbor>`, the bins are 1/2 the size of the minimum pair
cutoff.  Typically these are good values for minimizing the time for
neighbor list construction.  This setting overrides the default.
If you make it too big, there is little overhead due to
looping over bins, but more atoms are checked.  If you make it too
small, the optimal number of atoms is checked, but bin overhead goes
up.  If you set the binsize to 0.0, LAMMPS will use the default
binsize of 1/2 the cutoff.

Restrictions
""""""""""""


If the "delay" setting is non-zero, then it must be a multiple of the
"every" setting.

The molecule/intra and molecule/inter exclude options can only be used
with atom styles that define molecule IDs.

The value of the *page* setting must be at least 10x larger than the
*one* setting.  This insures neighbor pages are not mostly empty
space.

Related commands
""""""""""""""""

:doc:`neighbor <neighbor>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

The option defaults are delay = 10, every = 1, check = yes, once = no,
cluster = no, include = all (same as no include option defined),
exclude = none, page = 100000, one = 2000, and binsize = 0.0.


