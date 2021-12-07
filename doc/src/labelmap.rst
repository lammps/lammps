.. index:: labelmap

labelmap command
==================

Syntax
""""""

.. code-block:: LAMMPS

   labelmap type-kind numeric-type type-label ... keyword arg

* type-kind = *atom* or *bond* or *angle* or *dihedral* or *improper*
* one or more numeric-type/type-label pairs may be specified
* zero or more keyword/arg pairs may be appended
* keyword = *mapID*

  .. parsed-literal::

       *mapID* arg = name
         name = ID for the label map

Examples
""""""""

.. code-block:: LAMMPS

   labelmap atom 3 carbon
   labelmap bond 1 [c1][c2] 2 [c1][hc]
   labelmap bond 1 [c1][c2] 2 [c1][hc] mapID myMap
   labelmap atom $(label(carbon)) C // change type label from 'carbon' to 'C'

Description
"""""""""""

Define alphanumeric type labels to associate with one or more numeric
atom, bond, angle, dihedral or improper types.  A collection of type
labels for all atom types, bond types, etc is stored as a label map.
As explained below, this command also allows definition of multiple
label maps by use of the mapID keyword.

The default label maps can also be defined by the :doc:`read_data <read_data>`
command when it reads these sections in a data file: Atom Type Labels,
Bond Type Labels, etc.  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for a general discussion of how type
labels can be used.

As explained on the Howto page, valid type labels can contain any
alphanumeric character, but cannot start with a number.  They can also
contain standard characters such as square brackets "[" and "]", dash
"-", underscore "_", plus "+" and equals "=" signs.  Note that type
labels cannot contain the symbols '#' or '*'.

A *labelmap* command can only modify the label map for one type-kind
(atom types, bond types, etc).  Any number of numeric-type/type-label
pairs may follow.  If a type label already exists for a given numeric
type, it will be overwritten.  Type labels must be unique; assigning
the same type label to multiple numeric types is not allowed.  Note
that there is no requirement that a type label be defined for every
numeric type.

It may sometimes be useful to define multiple label maps.  This can be
done using the *mapID* keyword.  The specified *name* is the mapID
assigned to the label map.  The ID of a label map can only contain
alphanumeric characters and underscores.  If the *mapID* keyword is
not used, the specified type labels are assigned to the default label
map, which has no mapID.

----------

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`write_data <write_data>`,
:doc:`molecule <molecule>`, :doc:`fix bond/react <fix_bond_react>`

Default
"""""""

If the mapID keyword is not used, specified type labels are assigned
to the default map.
