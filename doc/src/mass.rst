.. index:: mass

mass command
============

Syntax
""""""

.. parsed-literal::

   mass I value

* I = atom type (see asterisk form below)
* value = mass

Examples
""""""""

.. code-block:: LAMMPS

   mass 1 1.0
   mass * 62.5
   mass 2* 62.5

Description
"""""""""""

Set the mass for all atoms of one or more atom types.  Per-type mass
values can also be set in the :doc:`read_data <read_data>` data file
using the "Masses" keyword.  See the :doc:`units <units>` command for
what mass units to use.

The I index can be specified in one of two ways.  An explicit numeric
value can be used, as in the first example above.  Or a wild-card
asterisk can be used to set the mass for multiple atom types.  This
takes the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of
atom types, then an asterisk with no numeric values means all types
from 1 to N.  A leading asterisk means all types from 1 to n
(inclusive).  A trailing asterisk means all types from n to N
(inclusive).  A middle asterisk means all types from m to n
(inclusive).

A line in a :doc:`data file <read_data>` that follows the "Masses"
keyword specifies mass using the same format as the arguments of the
mass command in an input script, except that no wild-card asterisk can
be used.  For example, under the "Masses" section of a data file, the
line that corresponds to the first example above would be listed as

.. parsed-literal::

   1 1.0

Note that the mass command can only be used if the :doc:`atom style <atom_style>` requires per-type atom mass to be set.
Currently, all but the *sphere* and *ellipsoid* and *peri* styles do.
They require mass to be set for individual particles, not types.
Per-atom masses are defined in the data file read by the
:doc:`read_data <read_data>` command, or set to default values by the
:doc:`create_atoms <create_atoms>` command.  Per-atom masses can also be
set to new values by the :doc:`set mass <set>` or :doc:`set density <set>`
commands.

Also note that :doc:`pair_style eam <pair_eam>` and :doc:`pair_style bop <pair_bop>` commands define the masses of atom types in their
respective potential files, in which case the mass command is normally
not used.

If you define a :doc:`hybrid atom style <atom_style>` which includes one
(or more) sub-styles which require per-type mass and one (or more)
sub-styles which require per-atom mass, then you must define both.
However, in this case the per-type mass will be ignored; only the
per-atom mass will be used by LAMMPS.

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

All masses must be defined before a simulation is run.  They must also
all be defined before a :doc:`velocity <velocity>` or :doc:`fix shake <fix_shake>` command is used.

The mass assigned to any type or atom must be > 0.0.

Related commands
""""""""""""""""

none


Default
"""""""

none
