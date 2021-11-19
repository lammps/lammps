.. index:: labelmap

labelmap command
==================

Syntax
""""""

.. code-block:: LAMMPS

   labelmap interaction numeric-type type-label ... keyword values

* interaction = *atom* or *bond* or *angle* or *dihedral* or *improper*
* one or more numeric-type/type-label pairs may be appended
* keyword = *map/assign*

  .. parsed-literal::

       *map/assign* value = mapID
         mapID = name of auxiliary label map

Examples
""""""""

.. code-block:: LAMMPS

   labelmap atom 3 carbon
   labelmap bond 1 [c1][c2] 2 [c1][hc]

Description
"""""""""""

Define or replace the alphanumeric type label associated with a
numeric atom, bond, angle, dihedral or improper type. Type labels are
strings that can be assigned to each interaction type. Pairing each
numeric type used by LAMMPS with a character-based type can be
helpful in a variety of situations. For example, type labels can make
various inputs more readable and compatible with other inputs (data
files, molecule templates, etc.), particularly if your model's force
field uses alphanumeric atom types. The default type label map can
also be defined by :doc:`data files <read_data>` using the relevant
Type Label sections.

.. note::

   If substituting numeric types with type labels is currently
   supported by a given command, this feature will be mentioned on
   that command's doc page. Or, for input script commands, type labels
   can be processed using :doc:`variable <variable>` syntax using
   labelmap functions.

.. note::

   Type labels cannot start with a number.

The label map of only one type of interaction (atom, bond, etc.) may
be modified by a given *labelmap* command. Any number of
numeric-type/type-label pairs may follow. If a type label already
exists for a given numeric type, it will be overwritten. Types labels
must be unique; assigning the same type label to multiple numeric
types is not permitted. There does not need to be a type label
associated with every numeric type; in this case, those types without
a label must be referenced by their numeric type.

For certain features, it is necessary to be able to extract the
atom types that make up a given bond, angle, dihedral or improper. The
standard way to write multi-atom interaction types is to enclose the
constituent atom types in brackets. For example, a bond type between
atom type 'C' and 'H' is written as:

.. code-block:: LAMMPS

   bond_style harmonic
   bond_coeff [C][H] 80.0 1.2

In some circumstances, it may be helpful to define extra label maps in
addition to the default one. An auxiliary label map can be created and
populated using the *map/assign* keyword. Type labels of auxiliary
label maps can be referenced by prefixing the type label with the
mapID followed by "::". For example, a type label "C" assigned to a
label map named "Map2" can be created and referenced as follows:

.. code-block:: LAMMPS

   labelmap bond 1 [C][H] map/assign Map2
   bond_coeff Map2::[C][H] 80.0 1.2

Auxiliary label maps are not written to restart or data files. Due to
the syntax for auxiliary maps, type labels cannot contain the
substring "::". The default label map has no prefix.

----------

Restrictions
""""""""""""

Type labels cannot contain the symbols '#' or '*'. This command must
come after the simulation box is defined by a :doc:`read_data <read_data>`,
:doc:`read_restart <read_restart>`, or :doc:`create_box <create_box>`
command.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`write_data <write_data>`,
:doc:`molecule <molecule>`, :doc:`fix bond/react <fix_bond_react>`

Default
"""""""

The option default is no type label map.
