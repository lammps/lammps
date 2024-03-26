.. index:: compute composition/atom
.. index:: compute composition/atom/kk

compute composition/atom command
================================

Accelerator Variants: *composition/atom/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID composition/atom keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* composition/atom = style name of this compute command
* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *cutoff*
       *cutoff* value = distance cutoff

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all composition/atom

   compute 1 all composition/atom cutoff 9.0
   comm_modify cutoff 9.0


Description
"""""""""""

.. versionadded:: 21Nov2023

Define a computation that calculates a local composition vector for each
atom. For a central atom with :math:`M` neighbors within the neighbor cutoff sphere,
composition is defined as the number of atoms of a given type
(including the central atom) divided by (:math:`M+1`).  For a given central atom,
the sum of all compositions equals one.

.. note::

   This compute uses the number of atom types, not chemical species, assigned in
   :doc:`pair_coeff <pair_coeff>` command.  If an interatomic potential has two
   species (i.e., Cu and Ni) assigned to four different atom types in
   :doc:`pair_coeff <pair_coeff>` (i.e., 'Cu Cu Ni Ni'), the compute will
   output four fractional values.  In those cases, the user may desire an extra
   calculation step to consolidate per-type fractions into per-species fractions.
   This calculation can be conducted within LAMMPS using another compute such as
   :doc:`compute reduce <compute_reduce>`, an atom-style :doc:`variable`, or as a
   post-processing step.

----------

The optional keyword *cutoff* defines the distance cutoff used when
searching for neighbors. The default value is the cutoff specified by
the pair style. If no pair style is defined, then a cutoff must be
defined using this keyword. If the specified cutoff is larger than
that of the pair_style plus neighbor skin (or no pair style is
defined), the *comm_modify cutoff* option must also be set to match
that of the *cutoff* keyword.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.
   This is the default setting for the :doc:`special_bonds
   <special_bonds>` command, and means those pairwise interactions do
   not appear in the neighbor list.  Because this compute uses the
   neighbor list, it also means those pairs will not be included in
   the order parameter.  This difficulty can be circumvented by
   writing a dump file, and using the :doc:`rerun <rerun>` command to
   compute the order parameter for snapshots in the dump file.  The
   rerun script can use a :doc:`special_bonds <special_bonds>` command
   that includes all pairs in the neighbor list.

----------

Output info
"""""""""""

This compute calculates a per-atom array with :math:`1 + N` columns, where :math:`N`
is the number of atom types. The first column is a count of the number of atoms
used to calculate composition (including the central atom), and each subsequent
column indicates the fraction of that atom type within the cutoff sphere.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>`
doc page for an overview of LAMMPS output options.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This compute requires :doc:`neighbor styles 'bin' or 'nsq' <neighbor>`.

Related commands
""""""""""""""""

:doc:`comm_modify <comm_modify>`

Default
"""""""

The option defaults are *cutoff* = pair style cutoff.
