.. index:: neighbor

neighbor command
================

Syntax
""""""


.. parsed-literal::

   neighbor skin style

* skin = extra distance beyond force cutoff (distance units)
* style = *bin* or *nsq* or *multi*

Examples
""""""""


.. parsed-literal::

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
size particles.  For the *bin* style, the bin size is set to 1/2 of
the largest cutoff distance between any pair of atom types and a
single set of bins is defined to search over for all atom types.  This
can be inefficient if one pair of types has a very long cutoff, but
other type pairs have a much shorter cutoff.  For style *multi* the
bin size is set to 1/2 of the shortest cutoff distance and multiple
sets of bins are defined to search over for different atom types.
This imposes some extra setup overhead, but the searches themselves
may be much faster for the short-cutoff cases.  See the :doc:`comm_modify mode multi <comm_modify>` command for a communication option
that may also be beneficial for simulations of this kind.

The :doc:`neigh_modify <neigh_modify>` command has additional options
that control how often neighbor lists are built and which pairs are
stored in the list.

When a run is finished, counts of the number of neighbors stored in
the pairwise list and the number of times neighbor lists were built
are printed to the screen and log file.  See the :doc:`Run output <Run_output>` doc page for details.

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
