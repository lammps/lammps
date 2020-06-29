.. index:: compute coord/atom

compute coord/atom command
==========================

compute coord/atom/kk command
===================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID coord/atom cstyle args ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* coord/atom = style name of this compute command
* cstyle = *cutoff* or *orientorder*

  .. parsed-literal::

       *cutoff* args = cutoff [group group2-ID] typeN
         cutoff = distance within which to count coordination neighbors (distance units)
         group *group2-ID* = select group-ID to restrict which atoms to consider for coordination number (optional)
         typeN = atom type for Nth coordination count (see asterisk form below)
       *orientorder* args = orientorderID threshold
         orientorderID = ID of an orientorder/atom compute
         threshold = minimum value of the product of two "connected" atoms

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all coord/atom cutoff 2.0
   compute 1 all coord/atom cutoff 6.0 1 2
   compute 1 all coord/atom cutoff 6.0 2*4 5*8 *
   compute 1 solute coord/atom cutoff 2.0 group solvent
   compute 1 all coord/atom orientorder 2 0.5

Description
"""""""""""

This compute performs calculations between neighboring atoms to
determine a coordination value.  The specific calculation and the
meaning of the resulting value depend on the *cstyle* keyword used.

The *cutoff* cstyle calculates one or more traditional coordination
numbers for each atom.  A coordination number is defined as the number
of neighbor atoms with specified atom type(s), and optionally within
the specified group, that are within the specified cutoff distance from
the central atom. The compute group selects only the central atoms; all
neighboring atoms, unless selected by type, type range, or group option,
are included in the coordination number tally.

The optional *group* keyword allows to specify from which group atoms
contribute to the coordination number. Default setting is group 'all'.

The *typeN* keywords allow specification of which atom types
contribute to each coordination number.  One coordination number is
computed for each of the *typeN* keywords listed.  If no *typeN*
keywords are listed, a single coordination number is calculated, which
includes atoms of all types (same as the "\*" format, see below).

The *typeN* keywords can be specified in one of two ways.  An explicit
numeric value can be used, as in the second example above.  Or a
wild-card asterisk can be used to specify a range of atom types.  This
takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of
atom types, then an asterisk with no numeric values means all types
from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

The *orientorder* cstyle calculates the number of "connected" neighbor
atoms J around each central atom I.  For this *cstyle*\ , connected is
defined by the orientational order parameter calculated by the
:doc:`compute orientorder/atom <compute_orientorder_atom>` command.
This *cstyle* thus allows one to apply the ten Wolde's criterion to
identify crystal-like atoms in a system, as discussed in :ref:`ten Wolde <tenWolde1>`.

The ID of the previously specified :doc:`compute orientorder/atom <compute_orientorder_atom>` command is specified as
*orientorderID*\ .  The compute must invoke its *components* option to
calculate components of the *Ybar_lm* vector for each atoms, as
described in its documentation.  Note that orientorder/atom compute
defines its own criteria for identifying neighboring atoms.  If the
scalar product (*Ybar_lm(i)*,*Ybar_lm(j)*), calculated by the
orientorder/atom compute is larger than the specified *threshold*\ ,
then I and J are connected, and the coordination value of I is
incremented by one.

For all *cstyle* settings, all coordination values will be 0.0 for
atoms not in the specified compute group.

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
   those pairs will not be included in the coordination count.  One way
   to get around this, is to write a dump file, and use the
   :doc:`rerun <rerun>` command to compute the coordination for snapshots
   in the dump file.  The rerun script can use a
   :doc:`special_bonds <special_bonds>` command that includes all pairs in
   the neighbor list.

----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------

**Output info:**

For *cstyle* cutoff, this compute can calculate a per-atom vector or
array.  If single *type1* keyword is specified (or if none are
specified), this compute calculates a per-atom vector.  If multiple
*typeN* keywords are specified, this compute calculates a per-atom
array, with N columns.

For *cstyle* orientorder, this compute calculates a per-atom vector.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

The per-atom vector or array values will be a number >= 0.0, as
explained above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute cluster/atom <compute_cluster_atom>`
:doc:`compute orientorder/atom <compute_orientorder_atom>`

Default
"""""""

group = all

----------

.. _tenWolde1:

**(tenWolde)** P. R. ten Wolde, M. J. Ruiz-Montero, D. Frenkel,
J. Chem. Phys. 104, 9932 (1996).
