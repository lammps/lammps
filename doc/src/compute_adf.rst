.. index:: compute adf

compute adf command
===================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID adf Nbin itype1 jtype1 ktype1 Rjinner1 Rjouter1 Rkinner1 Rkouter1 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* adf = style name of this compute command
* Nbin = number of ADF bins
* itypeN = central atom type for Nth ADF histogram (see asterisk form below)
* jtypeN = J atom type for Nth ADF histogram (see asterisk form below)
* ktypeN = K atom type for Nth ADF histogram (see asterisk form below)
* RjinnerN =  inner radius of J atom shell for Nth ADF histogram (distance units)
* RjouterN =  outer radius of J atom shell for Nth ADF histogram (distance units)
* RkinnerN = inner radius of K atom shell for Nth ADF histogram (distance units)
* RkouterN =  outer radius of K atom shell for Nth ADF histogram (distance units)
* zero or one keyword/value pairs may be appended
* keyword = *ordinate*
  
  .. parsed-literal::
  
       *ordinate* value = *degree* or *radian* or *cosine*
         Choose the ordinate parameter for the histogram



Examples
""""""""


.. parsed-literal::

   compute 1 fluid adf 32 1 1 1 0.0 1.2 0.0 1.2 &
                          1 1 2 0.0 1.2 0.0 1.5 &
                          1 2 2 0.0 1.5 0.0 1.5 &
                          2 1 1 0.0 1.2 0.0 1.2 &
                          2 1 2 0.0 1.5 2.0 3.5 &
                          2 2 2 2.0 3.5 2.0 3.5
   compute 1 fluid adf 32 1\*2 1\*2 1\*2 0.5 3.5
   compute 1 fluid adf 32

Description
"""""""""""

Define a computation that calculates one or more angular distribution functions
(ADF) for a group of particles.  Each ADF is calculated in histogram form
by measuring the angle formed by a central atom and two neighbor atoms and
binning these angles into *Nbin* bins.
Only neighbors for which *Rinner* < *R* < *Router* are counted, where
*Rinner* and *Router* are specified separately for the first and second
neighbor atom in each requested ADF.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses a neighbor list, it also means
   those pairs will not be included in the ADF. This does not apply when
   using long-range coulomb interactions (\ *coul/long*\ , *coul/msm*\ ,
   *coul/wolf* or similar.  One way to get around this would be to set
   special\_bond scaling factors to very tiny numbers that are not exactly
   zero (e.g. 1.0e-50). Another workaround is to write a dump file, and
   use the :doc:`rerun <rerun>` command to compute the ADF for snapshots in
   the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

.. note::

   If you request any outer cutoff *Router* > force cutoff, or if no
   pair style is defined,  e.g. the :doc:`rerun <rerun>` command is being used to
   post-process a dump file of snapshots you must insure ghost atom information
   out to the largest value of *Router* + *skin* is communicated, via the
   :doc:`comm_modify cutoff <comm_modify>` command, else the ADF computation
   cannot be performed, and LAMMPS will give an error message.  The *skin* value
   is what is specified with the :doc:`neighbor <neighbor>` command.

The *itypeN*\ ,\ *jtypeN*\ ,\ *ktypeN* settings can be specified in one of two
ways.  An explicit numeric value can be used, as in the 1st example
above.  Or a wild-card asterisk can be used to specify a range of atom
types as in the 2nd example above.
This takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the
number of atom types, then an asterisk with no numeric values means
all types from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

If *itypeN*\ , *jtypeN*\ , and *ktypeN* are single values, as in the 1st example
above, this means that the ADF is computed where atoms of type *itypeN*
are the central atom, and neighbor atoms of type *jtypeN* and *ktypeN*
are forming the angle.  If any of *itypeN*\ , *jtypeN*\ , or *ktypeN*
represent a range of values via
the wild-card asterisk, as in the 2nd example above, this means that the
ADF is computed where atoms of any of the range of types represented
by *itypeN* are the central atom, and the angle is formed by two neighbors,
one neighbor in the range of types represented by *jtypeN* and another neighbor
in the range of types represented by *ktypeN*\ .

If no *itypeN*\ , *jtypeN*\ , *ktypeN* settings are specified, then
LAMMPS will generate a single ADF for all atoms in the group.
The inner cutoff is set to zero and the outer cutoff is set
to the force cutoff. If no pair\_style is specified, there is no
force cutoff and LAMMPS will give an error message. Note that
in most cases, generating an ADF for all atoms is not a good thing.
Such an ADF is both uninformative and
extremely expensive to compute.  For example, with liquid water
with a 10 A force cutoff, there are 80,000 angles per atom.
In addition, most of the interesting angular structure occurs for
neighbors that are the closest to the central atom, involving
just a few dozen angles.

Angles for each ADF are generated by double-looping over the list of
neighbors of each central atom I,
just as they would be in the force calculation for
a three-body potential such as :doc:`Stillinger-Weber <pair_sw>`.
The angle formed by central atom I and neighbor atoms J and K is included in an
ADF if the following criteria are met:

* atoms I,J,K are all in the specified compute group
* the distance between atoms I,J is between Rjinner and Rjouter
* the distance between atoms I,K is between Rkinner and Rkouter
* the type of the I atom matches itypeN (one or a range of types)
* atoms I,J,K are distinct
* the type of the J atom matches jtypeN (one or a range of types)
* the type of the K atom matches ktypeN (one or a range of types)

Each unique angle satisfying the above criteria is counted only once, regardless
of whether either or both of the neighbor atoms making up the
angle appear in both the J and K lists.
It is OK if a particular angle is included in more than
one individual histogram, due to the way the *itypeN*\ , *jtypeN*\ , *ktypeN*
arguments are specified.

The first ADF value for a bin is calculated from the histogram count by
dividing by the total number of triples satisfying the criteria,
so that the integral of the ADF w.r.t. angle is 1, i.e. the ADF
is a probability density function.

The second ADF value is reported as a cumulative sum of
all bins up to the current bins, averaged
over atoms of type *itypeN*\ . It represents the
number of angles per central atom with angle less
than or equal to the angle of the current bin,
analogous to the coordination
number radial distribution function.

The *ordinate* optional keyword determines
whether the bins are of uniform angular size from zero
to 180 (\ *degree*\ ), zero to Pi (\ *radian*\ ), or the
cosine of the angle uniform in the range [-1,1] (\ *cosine*\ ).
*cosine* has the advantage of eliminating the *acos()* function
call, which speeds up the compute by 2-3x, and it is also preferred
on physical grounds, because the for uniformly distributed particles
in 3D, the angular probability density w.r.t dtheta is
sin(theta)/2, while for d(cos(theta)), it is 1/2,
Regardless of which ordinate is chosen, the first column of ADF
values is normalized w.r.t. the range of that ordinate, so that
the integral is 1.

The simplest way to output the results of the compute adf calculation
to a file is to use the :doc:`fix ave/time <fix_ave_time>` command, for
example:


.. parsed-literal::

   compute myADF all adf 32 2 2 2 0.5 3.5 0.5 3.5
   fix 1 all ave/time 100 1 100 c_myADF[\*] file tmp.adf mode vector

**Output info:**

This compute calculates a global array with the number of rows =
*Nbins*\ , and the number of columns = 1 + 2\*Ntriples, where Ntriples is the
number of I,J,K triples specified.  The first column has the bin
coordinate (angle-related ordinate at midpoint of bin). Each subsequent column has
the two ADF values for a specific set of (\ *itypeN*\ ,\ *jtypeN*\ ,\ *ktypeN*\ )
interactions, as described above.  These values can be used
by any command that uses a global values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The array values calculated by this compute are all "intensive".

The first column of array values is the angle-related ordinate, either
the angle in degrees or radians, or the cosine of the angle.  Each
subsequent pair of columns gives the first and second kinds of ADF
for a specific set of (\ *itypeN*\ ,\ *jtypeN*\ ,\ *ktypeN*\ ). The values
in the first ADF column are normalized numbers >= 0.0,
whose integral w.r.t. the ordinate is 1,
i.e. the first ADF is a normalized probability distribution.
The values in the second ADF column are also numbers >= 0.0.
They are the cumulative density distribution of angles per atom.
By definition, this ADF is monotonically increasing from zero to
a maximum value equal to the average total number of
angles per atom satisfying the ADF criteria.

Restrictions
""""""""""""


The ADF is not computed for neighbors outside the force cutoff,
since processors (in parallel) don't know about atom coordinates for
atoms further away than that distance.  If you want an ADF for larger
distances, you can use the :doc:`rerun <rerun>` command to post-process
a dump file and set the cutoff for the potential to be longer in the
rerun script.  Note that in the rerun context, the force cutoff is
arbitrary, since you aren't running dynamics and thus are not changing
your model.

Related commands
""""""""""""""""

:doc:`compute rdf <compute_rdf>`, :doc:`fix ave/time <fix_ave_time>`, :doc:`compute_modify <compute_modify>`

Default
"""""""

The keyword default is ordinate = degree.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
