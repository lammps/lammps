.. index:: compute rdf

compute rdf command
===================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID rdf Nbin itype1 jtype1 itype2 jtype2 ... keyword/value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* rdf = style name of this compute command
* Nbin = number of RDF bins
* itypeN = central atom type for Nth RDF histogram (integer, type label, or asterisk form)
* jtypeN = distribution atom type for Nth RDF histogram (integer, type label, or asterisk form)
* zero or more keyword/value pairs may be appended
* keyword = *cutoff*

  .. parsed-literal::

       *cutoff* value = Rcut
         Rcut = cutoff distance for RDF computation (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all rdf 100
   compute 1 all rdf 100 1 1
   compute 1 all rdf 100 * 3 cutoff 5.0
   compute 1 fluid rdf 500 1 1 1 2 2 1 2 2
   compute 1 fluid rdf 500 1*3 2 5 *10 cutoff 3.5

Description
"""""""""""

Define a computation that calculates the radial distribution function
(RDF), also called :math:`g(r)`, and the coordination number for a group of
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
   using long-range coulomb interactions (\ *coul/long*, *coul/msm*,
   *coul/wolf* or similar.  One way to get around this would be to set
   special_bond scaling factors to very tiny numbers that are not exactly
   zero (e.g., :math:`1.0 \times 10^{-50}`).  Another workaround is to write a
   dump file, and use the :doc:`rerun <rerun>` command to compute the RDF for
   snapshots in the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

By default the RDF is computed out to the maximum force cutoff defined
by the :doc:`pair_style <pair_style>` command.  If the *cutoff* keyword
is used, then the RDF is computed accurately out to the *Rcut* :math:`> 0.0`
distance specified.

.. note::

   Normally, you should only use the *cutoff* keyword if no pair
   style is defined (e.g., the :doc:`rerun <rerun>` command is being used to
   post-process a dump file of snapshots) or if you really want the RDF
   for distances beyond the pair_style force cutoff and cannot easily
   post-process a dump file to calculate it.  This is because using the
   *cutoff* keyword incurs extra computation and possibly communication,
   which may slow down your simulation.  If you specify *Rcut* :math:`\le`
   force cutoff, you will force an additional neighbor list to be built at
   every timestep this command is invoked (or every reneighboring
   timestep, whichever is less frequent), which is inefficient.  LAMMPS
   will warn you if this is the case.  If you specify a *Rcut* > force
   cutoff, you must ensure ghost atom information out to *Rcut* + *skin*
   is communicated, via the :doc:`comm_modify cutoff <comm_modify>`
   command, else the RDF computation cannot be performed, and LAMMPS will
   give an error message.  The *skin* value is what is specified with the
   :doc:`neighbor <neighbor>` command.  In this case, you are forcing a
   large neighbor list to be built just for the RDF computation, and
   extra communication to be performed every timestep.

The *itypeN* and *jtypeN* arguments are optional.  These arguments
must come in pairs.  If no pairs are listed, then a single histogram
is computed for :math:`g(r)` between all atom types.  If one or more pairs are
listed, then a separate histogram is generated for each
*itype*,\ *jtype* pair.

The *itypeN* and *jtypeN* settings can be specified in one of three
ways.  One or both of the types in the I,J pair can be a
:doc:`type label <Howto_type_labels>`.  Or an explicit numeric value can be
used, as in the fourth example above. Or a wild-card asterisk can be used
to specify a range of atom types.  This takes the form "\*" or "\*n" or
"m\*" or "m\*n". If :math:`N` is the number of atom types, then an asterisk
with no numeric values means all types from 1 to :math:`N`. A leading
asterisk means all types from 1 to n (inclusive).  A trailing asterisk
means all types from m to :math:`N` (inclusive).  A middle asterisk means
all types from m to n (inclusive).

If both *itypeN* and *jtypeN* are single values, as in the fourth example
above, this means that a :math:`g(r)` is computed where atoms of type *itypeN*
are the central atom, and atoms of type *jtypeN* are the distribution
atom.  If either *itypeN* and *jtypeN* represent a range of values via
the wild-card asterisk, as in the fifth example above, this means that a
:math:`g(r)` is computed where atoms of any of the range of types represented
by *itypeN* are the central atom, and atoms of any of the range of
types represented by *jtypeN* are the distribution atom.

Pairwise distances are generated by looping over a pairwise neighbor
list, just as they would be in a :doc:`pair_style <pair_style>`
computation.  The distance between two atoms :math:`I` and :math:`J` is
included in a specific histogram if the following criteria are met:

* atoms :math:`I` and :math:`J` are both in the specified compute group
* the distance between atoms :math:`I` and :math:`J` is less than the maximum
  force cutoff
* the type of the :math:`I` atom matches *itypeN* (one or a range of types)
* the type of the :math:`J` atom matches *jtypeN* (one or a range of types)

It is OK if a particular pairwise distance is included in more than
one individual histogram, due to the way the *itypeN* and *jtypeN*
arguments are specified.

The :math:`g(r)` value for a bin is calculated from the histogram count by
scaling it by the idealized number of how many counts there would be
if atoms of type *jtypeN* were uniformly distributed.  Thus it
involves the count of *itypeN* atoms, the count of *jtypeN* atoms, the
volume of the entire simulation box, and the volume of the bin's thin
shell in 3d (or the area of the bin's thin ring in 2d).

A coordination number :math:`\mathrm{coord}(r)` is also calculated, which is
the number of atoms of type *jtypeN* within the current bin or closer, averaged
over atoms of type *itypeN*\ .  This is calculated as the area- or
volume-weighted sum of :math:`g(r)` values over all bins up to and including
the current bin, multiplied by the global average volume density of
atoms of type *jtypeN*.

The simplest way to output the results of the compute rdf calculation
to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for
example:

.. code-block:: LAMMPS

   compute myRDF all rdf 50
   fix 1 all ave/time 100 1 100 c_myRDF[*] file tmp.rdf mode vector

Output info
"""""""""""

This compute calculates a global array in which the number of rows is
*Nbins* and the number of columns is :math:`1 + 2N_\text{pairs}`, where
:math:`N_\text{pairs}` is the number of :math:`I,J` pairings specified.
The first column has the bin coordinate (center of the bin), and each
successive set of two columns has the :math:`g(r)` and :math:`\text{coord}(r)`
values for a specific set of *itypeN* versus *jtypeN* interactions,
as described above.  These values can be used
by any command that uses a global values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "intensive".

The first column of array values will be in distance
:doc:`units <units>`.  The :math:`g(r)` columns of array values are normalized
numbers :math:`\ge 0.0`.  The coordination number columns of array values are
also numbers :math:`\ge 0.0`.

Restrictions
""""""""""""

By default, the RDF is not computed for distances longer than the
largest force cutoff, since the neighbor list creation will only contain
pairs up to that distance (plus neighbor list skin).  This distance can
be increased using the *cutoff* keyword but this keyword is only valid
with :doc:`neighbor styles 'bin' and 'nsq' <neighbor>`.

If you want an RDF for larger distances, you can also use the
:doc:`rerun <rerun>` command to post-process a dump file, use :doc:`pair
style zero <pair_zero>` and set the force cutoff to be longer in the
rerun script.  Note that in the rerun context, the force cutoff is
arbitrary and with pair style zero you are not computing any forces, and
you are not running dynamics you are not changing the model that
generated the trajectory.

The definition of :math:`g(r)` used by LAMMPS is only appropriate for
characterizing atoms that are uniformly distributed throughout the
simulation cell. In such cases, the coordination number is still correct
and meaningful.  As an example, if a large simulation cell contains only
one atom of type *itypeN* and one of *jtypeN*, then :math:`g(r)` will
register an arbitrarily large spike at whatever distance they happen to
be at, and zero everywhere else.  The function :math:`\text{coord}(r)`
will show a step change from zero to one at the location of the spike in
:math:`g(r)`.

.. note::

   compute rdf can handle dynamic groups and systems where atoms
   are added or removed, but this causes that certain normalization
   parameters need to be re-computed in every step and include collective
   communication operations. This will reduce performance and limit
   parallel efficiency and scaling. For systems, where only the type
   of atoms changes (e.g., when using :doc:`fix atom/swap <fix_atom_swap>`),
   you need to explicitly request the dynamic normalization updates
   via :doc:`compute_modify dynamic/dof yes <compute_modify>`

Related commands
""""""""""""""""

:doc:`fix ave/time <fix_ave_time>`, :doc:`compute_modify <compute_modify>`,
:doc:`compute adf <compute_adf>`

Default
"""""""

The keyword defaults are cutoff = 0.0 (use the pairwise force cutoff).
