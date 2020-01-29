.. index:: delete\_bonds

delete\_bonds command
=====================

Syntax
""""""


.. parsed-literal::

   delete_bonds group-ID style arg keyword ...

* group-ID = group ID
* style = *multi* or *atom* or *bond* or *angle* or *dihedral* or
          *improper* or *stats*
  
  .. parsed-literal::
  
       *multi* arg = none
       *atom* arg = an atom type or range of types (see below)
       *bond* arg = a bond type or range of types (see below)
       *angle* arg = an angle type or range of types (see below)
       *dihedral* arg = a dihedral type or range of types (see below)
       *improper* arg = an improper type or range of types (see below)
       *stats* arg = none

* zero or more keywords may be appended
* keyword = *any* or *undo* or *remove* or *special*


Examples
""""""""


.. parsed-literal::

   delete_bonds frozen multi remove
   delete_bonds all atom 4 special
   delete_bonds all bond 0\*3 special
   delete_bonds all stats

Description
"""""""""""

Turn off (or on) molecular topology interactions, i.e. bonds, angles,
dihedrals, impropers.  This command is useful for deleting
interactions that have been previously turned off by bond-breaking
potentials.  It is also useful for turning off topology interactions
between frozen or rigid atoms.  Pairwise interactions can be turned
off via the :doc:`neigh_modify exclude <neigh_modify>` command.  The
:doc:`fix shake <fix_shake>` command also effectively turns off certain
bond and angle interactions.

For all styles, by default, an interaction is only turned off (or on)
if all the atoms involved are in the specified group.  See the *any*
keyword to change the behavior.

Several of the styles (\ *atom*\ , *bond*\ , *angle*\ , *dihedral*\ ,
*improper*\ ) take a *type* as an argument.  The specified *type* should
be an integer from 0 to N, where N is the number of relevant types
(atom types, bond types, etc).  A value of 0 is only relevant for
style *bond*\ ; see details below.  In all cases, a wildcard asterisk
can be used in place of or in conjunction with the *type* argument to
specify a range of types.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the number of types, then an asterisk with no numeric
values means all types from 0 to N.  A leading asterisk means all
types from 0 to n (inclusive).  A trailing asterisk means all types
from n to N (inclusive).  A middle asterisk means all types from m to
n (inclusive).  Note that it is fine to include a type of 0 for
non-bond styles; it will simply be ignored.

For style *multi* all bond, angle, dihedral, and improper interactions
of any type, involving atoms in the group, are turned off.

Style *atom* is the same as style *multi* except that in addition, one
or more of the atoms involved in the bond, angle, dihedral, or
improper interaction must also be of the specified atom type.

For style *bond*\ , only bonds are candidates for turn-off, and the bond
must also be of the specified type.  Styles *angle*\ , *dihedral*\ , and
*improper* are treated similarly.

For style *bond*\ , you can set the type to 0 to delete bonds that have
been previously broken by a bond-breaking potential (which sets the
bond type to 0 when a bond is broken); e.g. see the :doc:`bond_style quartic <bond_style>` command.

For style *stats* no interactions are turned off (or on); the status
of all interactions in the specified group is simply reported.  This
is useful for diagnostic purposes if bonds have been turned off by a
bond-breaking potential during a previous run.

The default behavior of the delete\_bonds command is to turn off
interactions by toggling their type to a negative value, but not to
permanently remove the interaction.  E.g. a bond\_type of 2 is set to
-2.  The neighbor list creation routines will not include such an
interaction in their interaction lists.  The default is also to not
alter the list of 1-2, 1-3, 1-4 neighbors computed by the
:doc:`special_bonds <special_bonds>` command and used to weight pairwise
force and energy calculations.  This means that pairwise computations
will proceed as if the bond (or angle, etc) were still turned on.

Several keywords can be appended to the argument list to alter the
default behaviors.

The *any* keyword changes the requirement that all atoms in the bond
(angle, etc) must be in the specified group in order to turn-off the
interaction.  Instead, if any of the atoms in the interaction are in
the specified group, it will be turned off (or on if the *undo*
keyword is used).

The *undo* keyword inverts the delete\_bonds command so that the
specified bonds, angles, etc are turned on if they are currently
turned off.  This means a negative value is toggled to positive.  For
example, for style *angle*\ , if *type* is specified as 2, then all
angles with current type = -2, are reset to type = 2.  Note that the
:doc:`fix shake <fix_shake>` command also sets bond and angle types
negative, so this option should not be used on those interactions.

The *remove* keyword is invoked at the end of the delete\_bonds
operation.  It causes turned-off bonds (angles, etc) to be removed
from each atom's data structure and then adjusts the global bond
(angle, etc) counts accordingly.  Removal is a permanent change;
removed bonds cannot be turned back on via the *undo* keyword.
Removal does not alter the pairwise 1-2, 1-3, 1-4 weighting list.

The *special* keyword is invoked at the end of the delete\_bonds
operation, after (optional) removal.  It re-computes the pairwise 1-2,
1-3, 1-4 weighting list.  The weighting list computation treats
turned-off bonds the same as turned-on.  Thus, turned-off bonds must
be removed if you wish to change the weighting list.

Note that the choice of *remove* and *special* options affects how
1-2, 1-3, 1-4 pairwise interactions will be computed across bonds that
have been modified by the delete\_bonds command.

Restrictions
""""""""""""


This command requires inter-processor communication to acquire ghost
atoms, to coordinate the deleting of bonds, angles, etc between atoms
shared by multiple processors.  This means that your system must be
ready to perform a simulation before using this command (force fields
setup, atom masses set, etc).  Just as would be needed to run
dynamics, the force field you define should define a cutoff
(e.g. through a :doc:`pair_style <pair_style>` command) which is long
enough for a processor to acquire the ghost atoms its needs to compute
bond, angle, etc interactions.

If deleted bonds (angles, etc) are removed but the 1-2, 1-3, 1-4
weighting list is not re-computed, this can cause a later :doc:`fix shake <fix_shake>` command to fail due to an atom's bonds being
inconsistent with the weighting list.  This should only happen if the
group used in the fix command includes both atoms in the bond, in
which case you probably should be recomputing the weighting list.

Related commands
""""""""""""""""

:doc:`neigh_modify <neigh_modify>` exclude,
:doc:`special_bonds <special_bonds>`, :doc:`fix shake <fix_shake>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
