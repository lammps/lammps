.. index:: pair_style none

pair_style none command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style none

Examples
""""""""

.. code-block:: LAMMPS

   pair_style none

Description
"""""""""""

Using a pair style of none means pair forces and energies are not
computed.

With this choice, the force cutoff is 0.0, which means that only atoms
within the neighbor skin distance (see the :doc:`neighbor <neighbor>`
command) are communicated between processors.  You must insure the
skin distance is large enough to acquire atoms needed for computing
bonds, angles, etc.

A pair style of *none* will also prevent pairwise neighbor lists from
being built.  However if the :doc:`neighbor <neighbor>` style is *bin*\ ,
data structures for binning are still allocated.  If the neighbor skin
distance is small, then these data structures can consume a large
amount of memory.  So you should either set the neighbor style to
*nsq* or set the skin distance to a larger value.

See the :doc:`pair_style zero <pair_zero>` for a way to trigger the
building of a neighbor lists, but compute no pairwise interactions.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_style zero <pair_zero>`

**Default:** none
