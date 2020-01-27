.. index:: compute pair

compute pair command
====================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID pair pstyle [nstyle] [evalue]

* ID, group-ID are documented in :doc:`compute <compute>` command
* pair = style name of this compute command
* pstyle = style name of a pair style that calculates additional values
* nsub = *n*\ -instance of a sub-style, if a pair style is used multiple times in a hybrid style
* *evalue* = *epair* or *evdwl* or *ecoul* or blank (optional)


Examples
""""""""


.. parsed-literal::

   compute 1 all pair gauss
   compute 1 all pair lj/cut/coul/cut ecoul
   compute 1 all pair tersoff 2 epair
   compute 1 all pair reax/c

Description
"""""""""""

Define a computation that extracts additional values calculated by a
pair style, and makes them accessible for output or further processing
by other commands.

.. note::

   The group specified for this command is **ignored**\ .

The specified *pstyle* must be a pair style used in your simulation
either by itself or as a sub-style in a :doc:`pair_style hybrid or hybrid/overlay <pair_hybrid>` command. If the sub-style is
used more than once, an additional number *nsub* has to be specified
in order to choose which instance of the sub-style will be used by
the compute. Not specifying the number in this case will cause the
compute to fail.

The *evalue* setting is optional.  All
pair styles tally a potential energy *epair* which may be broken into
two parts: *evdwl* and *ecoul* such that *epair* = *evdwl* + *ecoul*\ .
If the pair style calculates Coulombic interactions, their energy will
be tallied in *ecoul*\ .  Everything else (whether it is a Lennard-Jones
style van der Waals interaction or not) is tallied in *evdwl*\ .  If
*evalue* is blank or specified as *epair*\ , then *epair* is stored
as a global scalar by this compute.  This is useful when using
:doc:`pair_style hybrid <pair_hybrid>` if you want to know the portion
of the total energy contributed by one sub-style.  If *evalue* is
specified as *evdwl* or *ecoul*\ , then just that portion of the energy
is stored as a global scalar.

.. note::

   The energy returned by the *evdwl* keyword does not include tail
   corrections, even if they are enabled via the
   :doc:`pair_modify <pair_modify>` command.

Some pair styles tally additional quantities, e.g. a breakdown of
potential energy into 14 components is tallied by the :doc:`pair_style reax/c <pair_reaxc>` command.  These values (1 or more)
are stored as a global vector by this compute.  See the doc page for
:doc:`individual pair styles <pair_style>` for info on these values.

**Output info:**

This compute calculates a global scalar which is *epair* or *evdwl* or
*ecoul*\ .  If the pair style supports it, it also calculates a global
vector of length >= 1, as determined by the pair style.  These values
can be used by any command that uses global scalar or vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The scalar and vector values calculated by this compute are
"extensive".

The scalar value will be in energy :doc:`units <units>`.  The vector
values will typically also be in energy :doc:`units <units>`, but see
the doc page for the pair style for details.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute pe <compute_pe>`, :doc:`compute bond <compute_bond>`

Default
"""""""

The keyword defaults are *evalue* = *epair*\ , nsub = 0.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
