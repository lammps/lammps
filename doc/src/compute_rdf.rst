.. index:: compute rdf

compute rdf command
===================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID rdf Nbin itype1 jtype1 itype2 jtype2 ... keyword/value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* rdf = style name of this compute command
* Nbin = number of RDF bins
* itypeN = central atom type for Nth RDF histogram (see asterisk form below)
* jtypeN = distribution atom type for Nth RDF histogram (see asterisk form below)
* zero or more keyword/value pairs may be appended
* keyword = *cutoff*

  .. parsed-literal::

       *cutoff* value = Rcut
         Rcut = cutoff distance for RDF computation (distance units)



Examples
""""""""


.. parsed-literal::

   compute 1 all rdf 100
   compute 1 all rdf 100 1 1
   compute 1 all rdf 100 \* 3 cutoff 5.0
   compute 1 fluid rdf 500 1 1 1 2 2 1 2 2
   compute 1 fluid rdf 500 1\*3 2 5 \*10 cutoff 3.5

Description
"""""""""""

Define a computation that calculates the radial distribution function
(RDF), also called g(r), and the coordination number for a group of
particles.  Both are calculated in histogram form by binning pairwise
distances into *Nbin* bins from 0.0 to the maximum force cutoff
defined by the :doc:`pair_style <pair_style>` command or the cutoff
distance *Rcut* specified via the *cutoff* keyword.  The bins are of
uniform size in radial distance.  Thus a single bin encompasses a thin
shell of distances in 3d and a thin ring of distances in 2d.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses a neighbor list, it also means
   those pairs will not be included in the RDF. This does not apply when
   using long-range coulomb interactions (\ *coul/long*\ , *coul/msm*\ ,
   *coul/wolf* or similar.  One way to get around this would be to set
   special\_bond scaling factors to very tiny numbers that are not exactly
   zero (e.g. 1.0e-50). Another workaround is to write a dump file, and
   use the :doc:`rerun <rerun>` command to compute the RDF for snapshots in
   the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

By default the RDF is computed out to the maximum force cutoff defined
by the :doc:`pair_style <pair_style>` command.  If the *cutoff* keyword
is used, then the RDF is computed accurately out to the *Rcut* > 0.0
distance specified.

.. note::

   Normally, you should only use the *cutoff* keyword if no pair
   style is defined, e.g. the :doc:`rerun <rerun>` command is being used to
   post-process a dump file of snapshots.  Or if you really want the RDF
   for distances beyond the pair\_style force cutoff and cannot easily
   post-process a dump file to calculate it.  This is because using the
   *cutoff* keyword incurs extra computation and possibly communication,
   which may slow down your simulation.  If you specify a *Rcut* <= force
   cutoff, you will force an additional neighbor list to be built at
   every timestep this command is invoked (or every reneighboring
   timestep, whichever is less frequent), which is inefficient.  LAMMPS
   will warn you if this is the case.  If you specify a *Rcut* > force
   cutoff, you must insure ghost atom information out to *Rcut* + *skin*
   is communicated, via the :doc:`comm_modify cutoff <comm_modify>`
   command, else the RDF computation cannot be performed, and LAMMPS will
   give an error message.  The *skin* value is what is specified with the
   :doc:`neighbor <neighbor>` command.  In this case, you are forcing a
   large neighbor list to be built just for the RDF computation, and
   extra communication to be performed every timestep.

The *itypeN* and *jtypeN* arguments are optional.  These arguments
must come in pairs.  If no pairs are listed, then a single histogram
is computed for g(r) between all atom types.  If one or more pairs are
listed, then a separate histogram is generated for each
*itype*\ ,\ *jtype* pair.

The *itypeN* and *jtypeN* settings can be specified in one of two
ways.  An explicit numeric value can be used, as in the 4th example
above.  Or a wild-card asterisk can be used to specify a range of atom
types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of atom types, then an asterisk with no numeric values means
all types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

If both *itypeN* and *jtypeN* are single values, as in the 4th example
above, this means that a g(r) is computed where atoms of type *itypeN*
are the central atom, and atoms of type *jtypeN* are the distribution
atom.  If either *itypeN* and *jtypeN* represent a range of values via
the wild-card asterisk, as in the 5th example above, this means that a
g(r) is computed where atoms of any of the range of types represented
by *itypeN* are the central atom, and atoms of any of the range of
types represented by *jtypeN* are the distribution atom.

Pairwise distances are generated by looping over a pairwise neighbor
list, just as they would be in a :doc:`pair_style <pair_style>`
computation.  The distance between two atoms I and J is included in a
specific histogram if the following criteria are met:

* atoms I,J are both in the specified compute group
* the distance between atoms I,J is less than the maximum force cutoff
* the type of the I atom matches itypeN (one or a range of types)
* the type of the J atom matches jtypeN (one or a range of types)

It is OK if a particular pairwise distance is included in more than
one individual histogram, due to the way the *itypeN* and *jtypeN*
arguments are specified.

The g(r) value for a bin is calculated from the histogram count by
scaling it by the idealized number of how many counts there would be
if atoms of type *jtypeN* were uniformly distributed.  Thus it
involves the count of *itypeN* atoms, the count of *jtypeN* atoms, the
volume of the entire simulation box, and the volume of the bin's thin
shell in 3d (or the area of the bin's thin ring in 2d).

A coordination number coord(r) is also calculated, which is the number
of atoms of type *jtypeN* within the current bin or closer, averaged
over atoms of type *itypeN*\ .  This is calculated as the area- or
volume-weighted sum of g(r) values over all bins up to and including
the current bin, multiplied by the global average volume density of
atoms of type jtypeN.

The simplest way to output the results of the compute rdf calculation
to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for
example:


.. parsed-literal::

   compute myRDF all rdf 50
   fix 1 all ave/time 100 1 100 c_myRDF[\*] file tmp.rdf mode vector

**Output info:**

This compute calculates a global array with the number of rows =
*Nbins*\ , and the number of columns = 1 + 2\*Npairs, where Npairs is the
number of I,J pairings specified.  The first column has the bin
coordinate (center of the bin), Each successive set of 2 columns has
the g(r) and coord(r) values for a specific set of *itypeN* versus
*jtypeN* interactions, as described above.  These values can be used
by any command that uses a global values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "intensive".

The first column of array values will be in distance
:doc:`units <units>`.  The g(r) columns of array values are normalized
numbers >= 0.0.  The coordination number columns of array values are
also numbers >= 0.0.

Restrictions
""""""""""""


The RDF is not computed for distances longer than the force cutoff,
since processors (in parallel) don't know about atom coordinates for
atoms further away than that distance.  If you want an RDF for larger
distances, you can use the :doc:`rerun <rerun>` command to post-process
a dump file and set the cutoff for the potential to be longer in the
rerun script.  Note that in the rerun context, the force cutoff is
arbitrary, since you aren't running dynamics and thus are not changing
your model.  The definition of g(r) used by LAMMPS is only appropriate
for characterizing atoms that are uniformly distributed throughout the
simulation cell. In such cases, the coordination number is still
correct and meaningful.  As an example, if a large simulation cell
contains only one atom of type *itypeN* and one of *jtypeN*\ , then g(r)
will register an arbitrarily large spike at whatever distance they
happen to be at, and zero everywhere else.  Coord(r) will show a step
change from zero to one at the location of the spike in g(r).

.. note::

   compute rdf can handle dynamic groups and systems where atoms
   are added or removed, but this causes that certain normalization
   parameters need to be re-computed in every step and include collective
   communication operations. This will reduce performance and limit
   parallel efficiency and scaling. For systems, where only the type
   of atoms changes (e.g. when using :doc:`fix atom/swap <fix_atom_swap>`),
   you need to explicitly request the dynamic normalization updates
   via :doc:`compute_modify dynamic yes <compute_modify>`

Related commands
""""""""""""""""

:doc:`fix ave/time <fix_ave_time>`, :doc:`compute_modify <compute_modify>`,
:doc:`compute adf <compute_adf>`

Default
"""""""

The keyword defaults are cutoff = 0.0 (use the pairwise force cutoff).
