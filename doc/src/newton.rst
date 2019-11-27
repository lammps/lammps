.. index:: newton

newton command
==============

Syntax
""""""


.. parsed-literal::

   newton flag
   newton flag1 flag2

* flag = *on* or *off* for both pairwise and bonded interactions
* flag1 = *on* or *off* for pairwise interactions
* flag2 = *on* or *off* for bonded interactions

Examples
""""""""


.. parsed-literal::

   newton off
   newton on off

Description
"""""""""""

This command turns Newton's 3rd law *on* or *off* for pairwise and
bonded interactions.  For most problems, setting Newton's 3rd law to
*on* means a modest savings in computation at the cost of two times
more communication.  Whether this is faster depends on problem size,
force cutoff lengths, a machine's compute/communication ratio, and how
many processors are being used.

Setting the pairwise newton flag to *off* means that if two
interacting atoms are on different processors, both processors compute
their interaction and the resulting force information is not
communicated.  Similarly, for bonded interactions, newton *off* means
that if a bond, angle, dihedral, or improper interaction contains
atoms on 2 or more processors, the interaction is computed by each
processor.

LAMMPS should produce the same answers for any newton flag settings,
except for round-off issues.

With :doc:`run\_style <run_style>` *respa* and only bonded interactions
(bond, angle, etc) computed in the innermost timestep, it may be
faster to turn newton *off* for bonded interactions, to avoid extra
communication in the innermost loop.

Restrictions
""""""""""""


The newton bond setting cannot be changed after the simulation box is
defined by a :doc:`read\_data <read_data>` or
:doc:`create\_box <create_box>` command.

Related commands
""""""""""""""""

:doc:`run\_style <run_style>` respa

Default
"""""""


.. parsed-literal::

   newton on


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
