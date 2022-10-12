.. index:: neighbor

neighbor command
================

Syntax
""""""

.. parsed-literal::

   neighbor skin style

* skin = extra distance beyond force cutoff (distance units)
* style = *bin* or *nsq* or *multi* or *multi/old*

Examples
""""""""

.. code-block:: LAMMPS

   neighbor 0.3 bin
   neighbor 2.0 nsq

Description
"""""""""""

This command sets parameters that affect the building of pairwise
neighbor lists.  All atom pairs within a neighbor cutoff distance
equal to the their force cutoff plus the *skin* distance are stored in
the list.  Typically, the larger the skin distance, the less often
neighbor lists need to be built, but more pairs must be checked for
possible force interactions every timestep.  The default value for
*skin* depends on the choice of units for the simulation; see the
default values below.

The *skin* distance is also used to determine how often atoms migrate
to new processors if the *check* option of the
:doc:`neigh_modify <neigh_modify>` command is set to *yes*\ .  Atoms are
migrated (communicated) to new processors on the same timestep that
neighbor lists are re-built.

The *style* value selects what algorithm is used to build the list.
The *bin* style creates the list by binning which is an operation that
scales linearly with N/P, the number of atoms per processor where N =
total number of atoms and P = number of processors.  It is almost
always faster than the *nsq* style which scales as (N/P)\^2.  For
unsolvated small molecules in a non-periodic box, the *nsq* choice can
sometimes be faster.  Either style should give the same answers.

The *multi* style is a modified binning algorithm that is useful for
systems with a wide range of cutoff distances, e.g. due to different
size particles. For granular pair styles, cutoffs are set to the sum of
the maximum atomic radii for each atom type.  For the *bin* style, the
bin size is set to 1/2 of the largest cutoff distance between any pair
of atom types and a single set of bins is defined to search over for all
atom types.  This can be inefficient if one pair of types has a very
long cutoff, but other type pairs have a much shorter cutoff. The
*multi* style uses different sized bins for collections of different
sized particles, where "size" may mean the physical size of the particle
or its cutoff distance for interacting with other particles. Different
sets of bins are then used to construct the neighbor lists as as further
described by Shire, Hanley, and Stratford :ref:`(Shire) <bytype-Shire>`.
This imposes some extra setup overhead, but the searches themselves may
be much faster. By default, each atom type defines a separate collection
of particles. For systems where two or more atom types have the same
size (either physical size or cutoff distance), the definition of
collections can be customized, which can result in less overhead and
faster performance. See the :doc:`neigh_modify <neigh_modify>` command
for how to define custom collections. Whether the collection definition
is customized or not, also see the :doc:`comm_modify mode multi
<comm_modify>` command for communication options that further improve
performance in a manner consistent with neighbor style multi.

An alternate style, *multi/old*, sets the bin size to 1/2 of the shortest
cutoff distance and multiple sets of bins are defined to search over for
different atom types. This algorithm used to be the default *multi*
algorithm in LAMMPS but was found to be significantly slower than the new
approach. For now we are keeping the old option in case there are use cases
where multi/old outperforms the new multi style.

.. note::

   If there are multiple sub-styles in a :doc:`hybrid/overlay pair style
   <pair_hybrid>` that cover the same atom types, but have significantly
   different cutoffs, the *multi* style does not apply.  Instead, the
   :doc:`pair_modify neigh/trim <pair_modify>` setting applies (which is
   *yes* by default).  Please check the neighbor list summary printed at
   the beginning of a calculation to verify that the desired set of
   neighbor list builds is performed.


The :doc:`neigh_modify <neigh_modify>` command has additional options
that control how often neighbor lists are built and which pairs are
stored in the list.

When a run is finished, counts of the number of neighbors stored in
the pairwise list and the number of times neighbor lists were built
are printed to the screen and log file.  See the :doc:`Run output <Run_output>` page for details.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`neigh_modify <neigh_modify>`, :doc:`units <units>`,
:doc:`comm_modify <comm_modify>`

Default
"""""""

| 0.3 bin for units = lj, skin = 0.3 sigma
| 2.0 bin for units = real or metal, skin = 2.0 Angstroms
| 0.001 bin for units = si, skin = 0.001 meters = 1.0 mm
| 0.1 bin for units = cgs, skin = 0.1 cm = 1.0 mm
|

----------

.. _bytype-Shire:

**(Shire)** Shire, Hanley and Stratford, Comp Part Mech, (2020).
