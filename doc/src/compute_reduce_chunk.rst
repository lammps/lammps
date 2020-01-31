.. index:: compute reduce/chunk

compute reduce/chunk command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID reduce/chunk chunkID mode input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* reduce/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* mode = *sum* or *min* or *max*
* one or more inputs can be listed
* input = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_ID
  
  .. parsed-literal::
  
       c_ID = per-atom vector calculated by a compute with ID
       c_ID[I] = Ith column of per-atom array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = per-atom vector calculated by a fix with ID
       f_ID[I] = Ith column of per-atom array calculated by a fix with ID, I can include wildcard (see below)
       v_name = per-atom vector calculated by an atom-style variable with name



Examples
""""""""


.. parsed-literal::

   compute 1 all reduce/chunk/atom mychunk min c_cluster

Description
"""""""""""

Define a calculation that reduces one or more per-atom vectors into
per-chunk values.  This can be useful for diagnostic output.  Or when
used in conjunction with the :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>` command it can be
used ot create per-atom values that induce a new set of chunks with a
second :doc:`compute chunk/atom <compute_chunk_atom>` command.  An
example is given below.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

For each atom, this compute accesses its chunk ID from the specified
*chunkID* compute.  The per-atom value from an input contributes
to a per-chunk value corresponding the the chunk ID.

The reduction operation is specified by the *mode* setting and is
performed over all the per-atom values from the atoms in each chunk.
The *sum* option adds the pre-atom values to a per-chunk total.  The
*min* or *max* options find the minimum or maximum value of the
per-atom values for each chunk.

Note that only atoms in the specified group contribute to the
reduction operation.  If the *chunkID* compute returns a 0 for the
chunk ID of an atom (i.e. the atom is not in a chunk defined by the
:doc:`compute chunk/atom <compute_chunk_atom>` command), that atom will
also not contribute to the reduction operation.  An input that is a
compute or fix may define its own group which affects the quantities
it returns.  For example, a compute with return a zero value for atoms
that are not in the group specified for that compute.

Each listed input is operated on independently.  Each input can be the
result of a :doc:`compute <compute>` or :doc:`fix <fix>` or the evaluation
of an atom-style :doc:`variable <variable>`.

Note that for values from a compute or fix, the bracketed index I can
be specified using a wildcard asterisk with the index to effectively
specify multiple values.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the size of the vector (for *mode* = scalar) or the
number of columns in the array (for *mode* = vector), then an asterisk
with no numeric values means all indices from 1 to N.  A leading
asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to N (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 compute reduce/chunk
commands are equivalent, since the :doc:`compute property/chunk <compute_property_chunk>` command creates a per-atom
array with 3 columns:


.. parsed-literal::

   compute prop all property/atom vx vy vz
   compute 10 all reduce/chunk mychunk max c_prop[\*]
   compute 10 all reduce/chunk mychunk max c_prop[1] c_prop[2] c_prop[3]


----------


Here is an example of using this compute, in conjunction with the
compute chunk/spread/atom command to identify self-assembled micelles.
The commands below can be added to the examples/in.micelle script.

Imagine a collection of polymer chains or small molecules with
hydrophobic end groups.  All the hydrophobic (HP) atoms are assigned
to a group called "phobic".

These commands will assign a unique cluster ID to all HP atoms within
a specified distance of each other.  A cluster will contain all HP
atoms in a single molecule, but also the HP atoms in nearby molecules,
e.g. molecules that have clumped to form a micelle due to the
attraction induced by the hydrophobicity.  The output of the
chunk/reduce command will be a cluster ID per chunk (molecule).
Molecules with the same cluster ID are in the same micelle.


.. parsed-literal::

   group phobic type 4     # specific to in.micelle model
   compute cluster phobic cluster/atom 2.0
   compute cmol all chunk/atom molecule
   compute reduce phobic reduce/chunk cmol min c_cluster

This per-chunk info could be output in at least two ways:


.. parsed-literal::

   fix 10 all ave/time 1000 1 1000 c_reduce file tmp.phobic mode vector

   compute spread all chunk/spread/atom cmol c_reduce
   dump 1 all custom 1000 tmp.dump id type mol x y z c_cluster c_spread
   dump_modify 1 sort id

In the first case, each snapshot in the tmp.phobic file will contain
one line per molecule.  Molecules with the same value are in the same
micelle.  In the second case each dump snapshot contains all atoms,
each with a final field with the cluster ID of the micelle that the HP
atoms of that atom's molecule belong to.

The result from compute chunk/spread/atom can be used to define a new
set of chunks, where all the atoms in all the molecules in the same
micelle are assigned to the same chunk, i.e. one chunk per micelle.


.. parsed-literal::

   compute micelle all chunk/atom c_spread compress yes

Further analysis on a per-micelle basis can now be performed using any
of the per-chunk computes listed on the :doc:`Howto chunk <Howto_chunk>`
doc page.  E.g. count the number of atoms in each micelle, calculate
its center or mass, shape (moments of inertia), radius of gyration,
etc.


.. parsed-literal::

   compute prop all property/chunk micelle count
   fix 20 all ave/time 1000 1 1000 c_prop file tmp.micelle mode vector

Each snapshot in the tmp.micelle file will have one line per micelle
with its count of atoms, plus a first line for a chunk with all the
solvent atoms.  By the time 50000 steps have elapsed there are a
handful of large micelles.


----------


**Output info:**

This compute calculates a global vector if a single input value is
specified, otherwise a global array is output.  The number of columns
in the array is the number of inputs provided.  The length of the
vector or the number of vector elements or array rows = the number of
chunks *Nchunk* as calculated by the specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  The vector or array can
be accessed by any command that uses global values from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The per-atom values for the vector or each column of the array will be
in whatever :doc:`units <units>` the corresponding input value is in.
The vector or array values are "intensive".

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute chunk/atom <compute_chunk_atom>`, :doc:`compute reduce <compute_reduce>`, :doc:`compute chunk/spread/atom <compute_chunk_spread_atom>`

**Default:** none


