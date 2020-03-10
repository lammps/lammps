.. index:: compute chunk/atom

compute chunk/atom command
==========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID chunk/atom style args keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* chunk/atom = style name of this compute command

  .. parsed-literal::

     style = *bin/1d* or *bin/2d* or *bin/3d* or *bin/sphere* or *type* or *molecule* or c_ID, c_ID[I], f_ID, f_ID[I], v_name
       *bin/1d* args = dim origin delta
         dim = *x* or *y* or *z*
         origin = *lower* or *center* or *upper* or coordinate value (distance units)
         delta = thickness of spatial bins in dim (distance units)
       *bin/2d* args = dim origin delta dim origin delta
         dim = *x* or *y* or *z*
         origin = *lower* or *center* or *upper* or coordinate value (distance units)
         delta = thickness of spatial bins in dim (distance units)
       *bin/3d* args = dim origin delta dim origin delta dim origin delta
         dim = *x* or *y* or *z*
         origin = *lower* or *center* or *upper* or coordinate value (distance units)
         delta = thickness of spatial bins in dim (distance units)
       *bin/sphere* args = xorig yorig zorig rmin rmax nsbin
         xorig,yorig,zorig = center point of sphere
         srmin,srmax = bin from sphere radius rmin to rmax
         nsbin = # of spherical shell bins between rmin and rmax
       *bin/cylinder* args = dim origin delta c1 c2 rmin rmax ncbin
         dim = *x* or *y* or *z* = axis of cylinder axis
         origin = *lower* or *center* or *upper* or coordinate value (distance units)
         delta = thickness of spatial bins in dim (distance units)
         c1,c2 = coords of cylinder axis in other 2 dimensions (distance units)
         crmin,crmax = bin from cylinder radius rmin to rmax (distance units)
         ncbin = # of concentric circle bins between rmin and rmax
       *type* args = none
       *molecule* args = none
       c_ID, c_ID[I], f_ID, f_ID[I], v_name args = none
         c_ID = per-atom vector calculated by a compute with ID
         c_ID[I] = Ith column of per-atom array calculated by a compute with ID
         f_ID = per-atom vector calculated by a fix with ID
         f_ID[I] = Ith column of per-atom array calculated by a fix with ID
         v_name = per-atom vector calculated by an atom-style variable with name

* zero or more keyword/values pairs may be appended
* keyword = *region* or *nchunk* or *static* or *compress* or *bound* or *discard* or *pbc* or *units*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to be part of a chunk
       *nchunk* value = *once* or *every*
         once = only compute the number of chunks once
         every = re-compute the number of chunks whenever invoked
       *limit* values = 0 or Nc max or Nc exact
         0 = no limit on the number of chunks
         Nc max = limit number of chunks to be <= Nc
         Nc exact = set number of chunks to exactly Nc
       *ids* value = *once* or *nfreq* or *every*
         once = assign chunk IDs to atoms only once, they persist thereafter
         nfreq = assign chunk IDs to atoms only once every Nfreq steps (if invoked by :doc:`fix ave/chunk <fix_ave_chunk>` which sets Nfreq)
         every = assign chunk IDs to atoms whenever invoked
       *compress* value = *yes* or *no*
         yes = compress chunk IDs to eliminate IDs with no atoms
         no = do not compress chunk IDs even if some IDs have no atoms
       *discard* value = *yes* or *no* or *mixed*
         yes = discard atoms with out-of-range chunk IDs by assigning a chunk ID = 0
         no = keep atoms with out-of-range chunk IDs by assigning a valid chunk ID
         mixed = keep or discard such atoms according to spatial binning rule
       *bound* values = x/y/z lo hi
         x/y/z = *x* or *y* or *z* to bound sptial bins in this dimension
         lo = *lower* or coordinate value (distance units)
         hi = *upper* or coordinate value (distance units)
       *pbc* value = *no* or *yes*
         yes = use periodic distance for bin/sphere and bin/cylinder styles
       *units* value = *box* or *lattice* or *reduced*



Examples
""""""""


.. parsed-literal::

   compute 1 all chunk/atom type
   compute 1 all chunk/atom bin/1d z lower 0.02 units reduced
   compute 1 all chunk/atom bin/2d z lower 1.0 y 0.0 2.5
   compute 1 all chunk/atom molecule region sphere nchunk once ids once compress yes
   compute 1 all chunk/atom bin/sphere 5 5 5 2.0 5.0 5 discard yes
   compute 1 all chunk/atom bin/cylinder z lower 2 10 10 2.0 5.0 3 discard yes
   compute 1 all chunk/atom c_cluster

Description
"""""""""""

Define a computation that calculates an integer chunk ID from 1 to
Nchunk for each atom in the group.  Values of chunk IDs are determined
by the *style* of chunk, which can be based on atom type or molecule
ID or spatial binning or a per-atom property or value calculated by
another :doc:`compute <compute>`, :doc:`fix <fix>`, or :doc:`atom-style variable <variable>`.  Per-atom chunk IDs can be used by other
computes with "chunk" in their style name, such as :doc:`compute com/chunk <compute_com_chunk>` or :doc:`compute msd/chunk <compute_msd_chunk>`.  Or they can be used by the :doc:`fix ave/chunk <fix_ave_chunk>` command to sum and time average a
variety of per-atom properties over the atoms in each chunk.  Or they
can simply be accessed by any command that uses per-atom values from a
compute as input, as discussed on the :doc:`Howto output <Howto_output>`
doc page.

See the :doc:`Howto chunk <Howto_chunk>` doc page for an overview of how
this compute can be used with a variety of other commands to tabulate
properties of a simulation.  The page gives several examples of input
script commands that can be used to calculate interesting properties.

Conceptually it is important to realize that this compute does two
simple things.  First, it sets the value of *Nchunk* = the number of
chunks, which can be a constant value or change over time.  Second, it
assigns each atom to a chunk via a chunk ID.  Chunk IDs range from 1
to *Nchunk* inclusive; some chunks may have no atoms assigned to them.
Atoms that do not belong to any chunk are assigned a value of 0.  Note
that the two operations are not always performed together.  For
example, spatial bins can be setup once (which sets *Nchunk*\ ), and
atoms assigned to those bins many times thereafter (setting their
chunk IDs).

All other commands in LAMMPS that use chunk IDs assume there are
*Nchunk* number of chunks, and that every atom is assigned to one of
those chunks, or not assigned to any chunk.

There are many options for specifying for how and when *Nchunk* is
calculated, and how and when chunk IDs are assigned to atoms.  The
details depend on the chunk *style* and its *args*\ , as well as
optional keyword settings.  They can also depend on whether a :doc:`fix ave/chunk <fix_ave_chunk>` command is using this compute, since
that command requires *Nchunk* to remain static across windows of
timesteps it specifies, while it accumulates per-chunk averages.

The details are described below.


----------


The different chunk styles operate as follows.  For each style, how it
calculates *Nchunk* and assigns chunk IDs to atoms is explained.  Note
that using the optional keywords can change both of those actions, as
described further below where the keywords are discussed.


----------


The *binning* styles perform a spatial binning of atoms, and assign an
atom the chunk ID corresponding to the bin number it is in.  *Nchunk*
is set to the number of bins, which can change if the simulation box
size changes.  This also depends on the setting of the *units*
keyword; e.g. for *reduced* units the number of chunks may not change
even if the box size does.

The *bin/1d*\ , *bin/2d*\ , and *bin/3d* styles define bins as 1d layers
(slabs), 2d pencils, or 3d boxes.  The *dim*\ , *origin*\ , and *delta*
settings are specified 1, 2, or 3 times.  For 2d or 3d bins, there is
no restriction on specifying dim = x before dim = y or z, or dim = y
before dim = z.  Bins in a particular *dim* have a bin size in that
dimension given by *delta*\ .  In each dimension, bins are defined
relative to a specified *origin*\ , which may be the lower/upper edge of
the simulation box (in that dimension), or its center point, or a
specified coordinate value.  Starting at the origin, sufficient bins
are created in both directions to completely span the simulation box
or the bounds specified by the optional *bounds* keyword.

For orthogonal simulation boxes, the bins are layers, pencils, or
boxes aligned with the xyz coordinate axes.  For triclinic
(non-orthogonal) simulation boxes, the bin faces are parallel to the
tilted faces of the simulation box.  See the :doc:`Howto triclinic <Howto_triclinic>` doc page for a discussion of the
geometry of triclinic boxes in LAMMPS.  As described there, a tilted
simulation box has edge vectors a,b,c.  In that nomenclature, bins in
the x dimension have faces with normals in the "b" cross "c"
direction.  Bins in y have faces normal to the "a" cross "c"
direction.  And bins in z have faces normal to the "a" cross "b"
direction.  Note that in order to define the size and position of
these bins in an unambiguous fashion, the *units* option must be set
to *reduced* when using a triclinic simulation box, as noted below.

The meaning of *origin* and *delta* for triclinic boxes is as follows.
Consider a triclinic box with bins that are 1d layers or slabs in the
x dimension.  No matter how the box is tilted, an *origin* of 0.0
means start layers at the lower "b" cross "c" plane of the simulation
box and an *origin* of 1.0 means to start layers at the upper "b"
cross "c" face of the box.  A *delta* value of 0.1 in *reduced* units
means there will be 10 layers from 0.0 to 1.0, regardless of the
current size or shape of the simulation box.

The *bin/sphere* style defines a set of spherical shell bins around
the origin (\ *xorig*\ ,\ *yorig*\ ,\ *zorig*\ ), using *nsbin* bins with radii
equally spaced between *srmin* and *srmax*\ .  This is effectively a 1d
vector of bins.  For example, if *srmin* = 1.0 and *srmax* = 10.0 and
*nsbin* = 9, then the first bin spans 1.0 < r < 2.0, and the last bin
spans 9.0 < r 10.0.  The geometry of the bins is the same whether the
simulation box is orthogonal or triclinic; i.e. the spherical shells
are not tilted or scaled differently in different dimensions to
transform them into ellipsoidal shells.

The *bin/cylinder* style defines bins for a cylinder oriented along
the axis *dim* with the axis coordinates in the other two radial
dimensions at (\ *c1*\ ,\ *c2*\ ).  For dim = x, c1/c2 = y/z; for dim = y,
c1/c2 = x/z; for dim = z, c1/c2 = x/y.  This is effectively a 2d array
of bins.  The first dimension is along the cylinder axis, the second
dimension is radially outward from the cylinder axis.  The bin size
and positions along the cylinder axis are specified by the *origin*
and *delta* values, the same as for the *bin/1d*\ , *bin/2d*\ , and
*bin/3d* styles.  There are *ncbin* concentric circle bins in the
radial direction from the cylinder axis with radii equally spaced
between *crmin* and *crmax*\ .  For example, if *crmin* = 1.0 and
*crmax* = 10.0 and *ncbin* = 9, then the first bin spans 1.0 < r <
2.0, and the last bin spans 9.0 < r 10.0.  The geometry of the bins in
the radial dimensions is the same whether the simulation box is
orthogonal or triclinic; i.e. the concentric circles are not tilted or
scaled differently in the two different dimensions to transform them
into ellipses.

The created bins (and hence the chunk IDs) are numbered consecutively
from 1 to the number of bins = *Nchunk*\ .  For *bin2d* and *bin3d*\ , the
numbering varies most rapidly in the first dimension (which could be
x, y, or z), next rapidly in the 2nd dimension, and most slowly in the
3rd dimension.  For *bin/sphere*\ , the bin with smallest radii is chunk
1 and the bni with largest radii is chunk Nchunk = *ncbin*\ .  For
*bin/cylinder*\ , the numbering varies most rapidly in the dimension
along the cylinder axis and most slowly in the radial direction.

Each time this compute is invoked, each atom is mapped to a bin based
on its current position.  Note that between reneighboring timesteps,
atoms can move outside the current simulation box.  If the box is
periodic (in that dimension) the atom is remapping into the periodic
box for purposes of binning.  If the box in not periodic, the atom may
have moved outside the bounds of all bins.  If an atom is not inside
any bin, the *discard* keyword is used to determine how a chunk ID is
assigned to the atom.


----------


The *type* style uses the atom type as the chunk ID.  *Nchunk* is set
to the number of atom types defined for the simulation, e.g. via the
:doc:`create_box <create_box>` or :doc:`read_data <read_data>` commands.


----------


The *molecule* style uses the molecule ID of each atom as its chunk
ID.  *Nchunk* is set to the largest chunk ID.  Note that this excludes
molecule IDs for atoms which are not in the specified group or
optional region.

There is no requirement that all atoms in a particular molecule are
assigned the same chunk ID (zero or non-zero), though you probably
want that to be the case, if you wish to compute a per-molecule
property.  LAMMPS will issue a warning if that is not the case, but
only the first time that *Nchunk* is calculated.

Note that atoms with a molecule ID = 0, which may be non-molecular
solvent atoms, have an out-of-range chunk ID.  These atoms are
discarded (not assigned to any chunk) or assigned to *Nchunk*\ ,
depending on the value of the *discard* keyword.


----------


The *compute/fix/variable* styles set the chunk ID of each atom based
on a quantity calculated and stored by a compute, fix, or variable.
In each case, it must be a per-atom quantity.  In each case the
referenced floating point values are converted to an integer chunk ID
as follows.  The floating point value is truncated (rounded down) to
an integer value.  If the integer value is <= 0, then a chunk ID of 0
is assigned to the atom.  If the integer value is > 0, it becomes the
chunk ID to the atom.  *Nchunk* is set to the largest chunk ID.  Note
that this excludes atoms which are not in the specified group or
optional region.

If the style begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If no bracketed integer is
appended, the per-atom vector calculated by the compute is used.  If a
bracketed integer is appended, the Ith column of the per-atom array
calculated by the compute is used.  Users can also write code for
their own compute styles and :doc:`add them to LAMMPS <Modify>`.

If the style begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If no bracketed integer is
appended, the per-atom vector calculated by the fix is used.  If a
bracketed integer is appended, the Ith column of the per-atom array
calculated by the fix is used.  Note that some fixes only produce
their values on certain timesteps, which must be compatible with the
timestep on which this compute accesses the fix, else an error
results.  Users can also write code for their own fix styles and :doc:`add them to LAMMPS <Modify>`.

If a value begins with "v\_", a variable name for an *atom* or
*atomfile* style :doc:`variable <variable>` must follow which has been
previously defined in the input script.  Variables of style *atom* can
reference thermodynamic keywords and various per-atom attributes, or
invoke other computes, fixes, or variables when they are evaluated, so
this is a very general means of generating per-atom quantities to
treat as a chunk ID.


----------


Normally, *Nchunk* = the number of chunks, is re-calculated every time
this fix is invoked, though the value may or may not change.  As
explained below, the *nchunk* keyword can be set to *once* which means
*Nchunk* will never change.

If a :doc:`fix ave/chunk <fix_ave_chunk>` command uses this compute, it
can also turn off the re-calculation of *Nchunk* for one or more
windows of timesteps.  The extent of the windows, during which Nchunk
is held constant, are determined by the *Nevery*\ , *Nrepeat*\ , *Nfreq*
values and the *ave* keyword setting that are used by the :doc:`fix ave/chunk <fix_ave_chunk>` command.

Specifically, if *ave* = *one*\ , then for each span of *Nfreq*
timesteps, *Nchunk* is held constant between the first timestep when
averaging is done (within the Nfreq-length window), and the last
timestep when averaging is done (multiple of Nfreq).  If *ave* =
*running* or *window*\ , then *Nchunk* is held constant forever,
starting on the first timestep when the :doc:`fix ave/chunk <fix_ave_chunk>` command invokes this compute.

Note that multiple :doc:`fix ave/chunk <fix_ave_chunk>` commands can use
the same compute chunk/atom compute.  However, the time windows they
induce for holding *Nchunk* constant must be identical, else an error
will be generated.


----------


The various optional keywords operate as follows.  Note that some of
them function differently or are ignored by different chunk styles.
Some of them also have different default values, depending on
the chunk style, as listed below.

The *region* keyword applies to all chunk styles.  If used, an atom
must be in both the specified group and the specified geometric
:doc:`region <region>` to be assigned to a chunk.


----------


The *nchunk* keyword applies to all chunk styles.  It specifies how
often *Nchunk* is recalculated, which in turn can affect the chunk IDs
assigned to individual atoms.

If *nchunk* is set to *once*\ , then *Nchunk* is only calculated once,
the first time this compute is invoked.  If *nchunk* is set to
*every*\ , then *Nchunk* is re-calculated every time the compute is
invoked.  Note that, as described above, the use of this compute
by the :doc:`fix ave/chunk <fix_ave_chunk>` command can override
the *every* setting.

The default values for *nchunk* are listed below and depend on the
chunk style and other system and keyword settings.  They attempt to
represent typical use cases for the various chunk styles.  The
*nchunk* value can always be set explicitly if desired.


----------


The *limit* keyword can be used to limit the calculated value of
*Nchunk* = the number of chunks.  The limit is applied each time
*Nchunk* is calculated, which also limits the chunk IDs assigned to
any atom.  The *limit* keyword is used by all chunk styles except the
*binning* styles, which ignore it.  This is because the number of bins
can be tailored using the *bound* keyword (described below) which
effectively limits the size of *Nchunk*\ .

If *limit* is set to *Nc* = 0, then no limit is imposed on *Nchunk*\ ,
though the *compress* keyword can still be used to reduce *Nchunk*\ , as
described below.

If *Nc* > 0, then the effect of the *limit* keyword depends on whether
the *compress* keyword is also used with a setting of *yes*\ , and
whether the *compress* keyword is specified before the *limit* keyword
or after.

In all cases, *Nchunk* is first calculated in the usual way for each
chunk style, as described above.

First, here is what occurs if *compress yes* is not set.  If *limit*
is set to *Nc max*\ , then *Nchunk* is reset to the smaller of *Nchunk*
and *Nc*\ .  If *limit* is set to *Nc exact*\ , then *Nchunk* is reset to
*Nc*\ , whether the original *Nchunk* was larger or smaller than *Nc*\ .
If *Nchunk* shrank due to the *limit* setting, then atom chunk IDs >
*Nchunk* will be reset to 0 or *Nchunk*\ , depending on the setting of
the *discard* keyword.  If *Nchunk* grew, there will simply be some
chunks with no atoms assigned to them.

If *compress yes* is set, and the *compress* keyword comes before the
*limit* keyword, the compression operation is performed first, as
described below, which resets *Nchunk*\ .  The *limit* keyword is then
applied to the new *Nchunk* value, exactly as described in the
preceding paragraph.  Note that in this case, all atoms will end up
with chunk IDs <= *Nc*\ , but their original values (e.g. molecule ID or
compute/fix/variable) may have been > *Nc*\ , because of the compression
operation.

If *compress yes* is set, and the *compress* keyword comes after the
*limit* keyword, then the *limit* value of *Nc* is applied first to
the uncompressed value of *Nchunk*\ , but only if *Nc* < *Nchunk*
(whether *Nc max* or *Nc exact* is used).  This effectively means all
atoms with chunk IDs > *Nc* have their chunk IDs reset to 0 or *Nc*\ ,
depending on the setting of the *discard* keyword.  The compression
operation is then performed, which may shrink *Nchunk* further.  If
the new *Nchunk* < *Nc* and *limit* = *Nc exact* is specified, then
*Nchunk* is reset to *Nc*\ , which results in extra chunks with no atoms
assigned to them.  Note that in this case, all atoms will end up with
chunk IDs <= *Nc*\ , and their original values (e.g. molecule ID or
compute/fix/variable value) will also have been <= *Nc*\ .


----------


The *ids* keyword applies to all chunk styles.  If the setting is
*once* then the chunk IDs assigned to atoms the first time this
compute is invoked will be permanent, and never be re-computed.

If the setting is *nfreq* and if a :doc:`fix ave/chunk <fix_ave_chunk>`
command is using this compute, then in each of the *Nchunk* = constant
time windows (discussed above), the chunk ID's assigned to atoms on
the first step of the time window will persist until the end of the
time window.

If the setting is *every*\ , which is the default, then chunk IDs are
re-calculated on any timestep this compute is invoked.

.. note::

   If you want the persistent chunk-IDs calculated by this compute
   to be continuous when running from a :doc:`restart file <read_restart>`,
   then you should use the same ID for this compute, as in the original
   run.  This is so that the fix this compute creates to store per-atom
   quantities will also have the same ID, and thus be initialized
   correctly with chunk IDs from the restart file.


----------


The *compress* keyword applies to all chunk styles and affects how
*Nchunk* is calculated, which in turn affects the chunk IDs assigned
to each atom.  It is useful for converting a "sparse" set of chunk IDs
(with many IDs that have no atoms assigned to them), into a "dense"
set of IDs, where every chunk has one or more atoms assigned to it.

Two possible use cases are as follows.  If a large simulation box is
mostly empty space, then the *binning* style may produce many bins
with no atoms.  If *compress* is set to *yes*\ , only bins with atoms
will be contribute to *Nchunk*\ .  Likewise, the *molecule* or
*compute/fix/variable* styles may produce large *Nchunk* values.  For
example, the :doc:`compute cluster/atom <compute_cluster_atom>` command
assigns every atom an atom ID for one of the atoms it is clustered
with.  For a million-atom system with 5 clusters, there would only be
5 unique chunk IDs, but the largest chunk ID might be 1 million,
resulting in *Nchunk* = 1 million.  If *compress* is set to *yes*\ ,
*Nchunk* will be reset to 5.

If *compress* is set to *no*\ , which is the default, no compression is
done.  If it is set to *yes*\ , all chunk IDs with no atoms are removed
from the list of chunk IDs, and the list is sorted.  The remaining
chunk IDs are renumbered from 1 to *Nchunk* where *Nchunk* is the new
length of the list.  The chunk IDs assigned to each atom reflect
the new renumbering from 1 to *Nchunk*\ .

The original chunk IDs (before renumbering) can be accessed by the
:doc:`compute property/chunk <compute_property_chunk>` command and its
*id* keyword, or by the :doc:`fix ave/chunk <fix_ave_chunk>` command
which outputs the original IDs as one of the columns in its global
output array.  For example, using the "compute cluster/atom" command
discussed above, the original 5 unique chunk IDs might be atom IDs
(27,4982,58374,857838,1000000).  After compression, these will be
renumbered to (1,2,3,4,5).  The original values (27,...,1000000) can
be output to a file by the :doc:`fix ave/chunk <fix_ave_chunk>` command,
or by using the :doc:`fix ave/time <fix_ave_time>` command in
conjunction with the :doc:`compute property/chunk <compute_property_chunk>` command.

.. note::

   The compression operation requires global communication across
   all processors to share their chunk ID values.  It can require large
   memory on every processor to store them, even after they are
   compressed, if there are a large number of unique chunk IDs with
   atoms assigned to them.  It uses a STL map to find unique chunk IDs
   and store them in sorted order.  Each time an atom is assigned a
   compressed chunk ID, it must access the STL map.  All of this means
   that compression can be expensive, both in memory and CPU time.  The
   use of the *limit* keyword in conjunction with the *compress* keyword
   can affect these costs, depending on which keyword is used first.  So
   use this option with care.


----------


The *discard* keyword applies to all chunk styles.  It affects what
chunk IDs are assigned to atoms that do not match one of the valid
chunk IDs from 1 to *Nchunk*\ .  Note that it does not apply to atoms
that are not in the specified group or optionally specified region.
Those atoms are always assigned a chunk ID = 0.

If the calculated chunk ID for an atom is not within the range 1 to
*Nchunk* then it is a "discard" atom.  Note that *Nchunk* may have
been shrunk by the *limit* keyword.  Or the *compress* keyword may
have eliminated chunk IDs that were valid before the compression took
place, and are now not in the compressed list.  Also note that for the
*molecule* chunk style, if new molecules are added to the system,
their chunk IDs may exceed a previously calculated *Nchunk*\ .
Likewise, evaluation of a compute/fix/variable on a later timestep may
return chunk IDs that are invalid for the previously calculated
*Nchunk*\ .

All the chunk styles except the *binning* styles, must use *discard*
set to either *yes* or *no*\ .  If *discard* is set to *yes*\ , which is
the default, then every "discard" atom has its chunk ID set to 0.  If
*discard* is set to *no*\ , every "discard" atom has its chunk ID set to
*Nchunk*\ .  I.e. it becomes part of the last chunk.

The *binning* styles use the *discard* keyword to decide whether to
discard atoms outside the spatial domain covered by bins, or to assign
them to the bin they are nearest to.

For the *bin/1d*\ , *bin/2d*\ , *bin/3d* styles the details are as
follows.  If *discard* is set to *yes*\ , an out-of-domain atom will
have its chunk ID set to 0.  If *discard* is set to *no*\ , the atom
will have its chunk ID set to the first or last bin in that dimension.
If *discard* is set to *mixed*\ , which is the default, it will only
have its chunk ID set to the first or last bin if bins extend to the
simulation box boundary in that dimension.  This is the case if the
*bound* keyword settings are *lower* and *upper*\ , which is the
default.  If the *bound* keyword settings are numeric values, then the
atom will have its chunk ID set to 0 if it is outside the bounds of
any bin.  Note that in this case, it is possible that the first or
last bin extends beyond the numeric *bounds* settings, depending on
the specified *origin*\ .  If this is the case, the chunk ID of the atom
is only set to 0 if it is outside the first or last bin, not if it is
simply outside the numeric *bounds* setting.

For the *bin/sphere* style the details are as follows.  If *discard*
is set to *yes*\ , an out-of-domain atom will have its chunk ID set to
0.  If *discard* is set to *no* or *mixed*\ , the atom will have its
chunk ID set to the first or last bin, i.e. the innermost or outermost
spherical shell.  If the distance of the atom from the origin is less
than *rmin*\ , it will be assigned to the first bin.  If the distance of
the atom from the origin is greater than *rmax*\ , it will be assigned
to the last bin.

For the *bin/cylinder* style the details are as follows.  If *discard*
is set to *yes*\ , an out-of-domain atom will have its chunk ID set to
0.  If *discard* is set to *no*\ , the atom will have its chunk ID set
to the first or last bin in both the radial and axis dimensions.  If
*discard* is set to *mixed*\ , which is the default, the radial
dimension is treated the same as for *discard* = no.  But for the axis
dimension, it will only have its chunk ID set to the first or last
bin if bins extend to the simulation box boundary in the axis
dimension.  This is the case if the *bound* keyword settings are
*lower* and *upper*\ , which is the default.  If the *bound* keyword
settings are numeric values, then the atom will have its chunk ID set
to 0 if it is outside the bounds of any bin.  Note that in this case,
it is possible that the first or last bin extends beyond the numeric
*bounds* settings, depending on the specified *origin*\ .  If this is
the case, the chunk ID of the atom is only set to 0 if it is outside
the first or last bin, not if it is simply outside the numeric
*bounds* setting.

If *discard* is set to *no* or *mixed*\ , the atom will have its
chunk ID set to the first or last bin, i.e. the innermost or outermost
spherical shell.  If the distance of the atom from the origin is less
than *rmin*\ , it will be assigned to the first bin.  If the distance of
the atom from the origin is greater than *rmax*\ , it will be assigned
to the last bin.


----------


The *bound* keyword only applies to the *bin/1d*\ , *bin/2d*\ , *bin/3d*
styles and to the axis dimension of the *bin/cylinder* style;
otherwise it is ignored.  It can be used one or more times to limit
the extent of bin coverage in a specified dimension, i.e. to only bin
a portion of the box.  If the *lo* setting is *lower* or the *hi*
setting is *upper*\ , the bin extent in that direction extends to the
box boundary.  If a numeric value is used for *lo* and/or *hi*\ , then
the bin extent in the *lo* or *hi* direction extends only to that
value, which is assumed to be inside (or at least near) the simulation
box boundaries, though LAMMPS does not check for this.  Note that
using the *bound* keyword typically reduces the total number of bins
and thus the number of chunks *Nchunk*\ .

The *pbc* keyword only applies to the *bin/sphere* and *bin/cylinder*
styles.  If set to *yes*\ , the distance an atom is from the sphere
origin or cylinder axis is calculated in a minimum image sense with
respect to periodic dimensions, when determining which bin the atom is
in.  I.e. if x is a periodic dimension and the distance between the
atom and the sphere center in the x dimension is greater than 0.5 \*
simulation box length in x, then a box length is subtracted to give a
distance < 0.5 \* simulation box length.  This allosws the sphere or
cylinder center to be near a box edge, and atoms on the other side of
the periodic box will still be close to the center point/axis.  Note
that with a setting of *yes*\ , the outer sphere or cylinder radius must
also be <= 0.5 \* simulation box length in any periodic dimension
except for the cylinder axis dimension, or an error is generated.

The *units* keyword only applies to the *binning* styles; otherwise it
is ignored.  For the *bin/1d*\ , *bin/2d*\ , *bin/3d* styles, it
determines the meaning of the distance units used for the bin sizes
*delta* and for *origin* and *bounds* values if they are coordinate
values.  For the *bin/sphere* style it determines the meaning of the
distance units used for *xorig*\ ,\ *yorig*\ ,\ *zorig* and the radii *srmin*
and *srmax*\ .  For the *bin/cylinder* style it determines the meaning
of the distance units used for *delta*\ ,\ *c1*\ ,\ *c2* and the radii *crmin*
and *crmax*\ .

For orthogonal simulation boxes, any of the 3 options may
be used.  For non-orthogonal (triclinic) simulation boxes, only the
*reduced* option may be used.

A *box* value selects standard distance units as defined by the
:doc:`units <units>` command, e.g. Angstroms for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacing.  A *reduced* value means normalized
unitless values between 0 and 1, which represent the lower and upper
faces of the simulation box respectively.  Thus an *origin* value of
0.5 means the center of the box in any dimension.  A *delta* value of
0.1 means 10 bins span the box in that dimension.

Note that for the *bin/sphere* style, the radii *srmin* and *srmax* are
scaled by the lattice spacing or reduced value of the *x* dimension.

Note that for the *bin/cylinder* style, the radii *crmin* and *crmax*
are scaled by the lattice spacing or reduced value of the 1st
dimension perpendicular to the cylinder axis.  E.g. y for an x-axis
cylinder, x for a y-axis cylinder, and x for a z-axis cylinder.


----------


**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values are unitless chunk IDs, ranging from 1 to
*Nchunk* (inclusive) for atoms assigned to chunks, and 0 for atoms not
belonging to a chunk.

Restrictions
""""""""""""


Even if the *nchunk* keyword is set to *once*\ , the chunk IDs assigned
to each atom are not stored in a restart files.  This means you cannot
expect those assignments to persist in a restarted simulation.
Instead you must re-specify this command and assign atoms to chunks when
the restarted simulation begins.

Related commands
""""""""""""""""

:doc:`fix ave/chunk <fix_ave_chunk>`,
:doc:`compute global/atom <compute_global_atom>`

Default
"""""""

The option defaults are as follows:

* region = none
* nchunk = every, if compress is yes, overriding other defaults listed here
* nchunk = once, for type style
* nchunk = once, for mol style if region is none
* nchunk = every, for mol style if region is set
* nchunk = once, for binning style if the simulation box size is static or units = reduced
* nchunk = every, for binning style if the simulation box size is dynamic and units is lattice or box
* nchunk = every, for compute/fix/variable style
* limit = 0
* ids = every
* compress = no
* discard = yes, for all styles except binning
* discard = mixed, for binning styles
* bound = lower and upper in all dimensions
* pbc = no
* units = lattice
