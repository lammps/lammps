.. index:: compute phase/atom
.. index:: compute phase/atom/kk

compute phase/atom command
================================

Accelerator Variants: *phase/atom/kk*

Syntax
""""""

.. parsed-literal::

   compute ID group-ID phase/atom keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* phase/atom = style name of this compute command
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *cutoff*
       *cutoff* value = distance cutoff

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all phase/atom

   compute 1 all phase/atom cutoff 5.0
   comm_modify cutoff 5.0

Description
"""""""""""

Define a computation that calculates the local density and temperature
for each atom and neighbors inside a spherical cutoff.

The optional keyword *cutoff* defines the distance cutoff
used when searching for neighbors. The default value is the cutoff
specified by the pair style. If no pair style is defined, then a cutoff
must be defined using this keyword. If the specified cutoff is larger than
that of the pair_style plus neighbor skin (or no pair style is defined),
the *comm_modify cutoff* option must also be set to match that of the
*cutoff* keyword.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included in the order parameter.  This
   difficulty can be circumvented by writing a dump file, and using the
   :doc:`rerun <rerun>` command to compute the order parameter for
   snapshots in the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

----------


.. include:: accel_styles.rst


----------

Output info
"""""""""""

This compute calculates a per-atom array with two columns: density and temperature.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`comm_modify <comm_modify>`

Default
"""""""

The option defaults are *cutoff* = pair style cutoff

