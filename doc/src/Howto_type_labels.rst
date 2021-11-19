Type labels
===========

Each atom in LAMMPS has an associated atom type.  Likewise each bond,
angle, dihedral, improper enumerated in a data file read by the
:doc:`read_data <read_data>` command has a bond type, angle type, etc.

By default, type values are integers from 1 to Ntypes, wherever they
appear in LAMMPS input or output files.  However, LAMMPS also has the
option to use alphanumeric strings, called type labels, for these
values in its input and output.

This can be useful in various scenarios.  For example, type labels can
make inputs more readable and compatible with other inputs (e.g. data
files, molecule template files read by the :doc:`molecule <molecule>`
command), particularly if a force field uses alphanumeric atom types.
See a list of other commands below which can use label types in
interesting ways.

A collection of one or more type labels for one category of types
(atom types, bond types, etc) is stored as a "label map" which is
simply a list of numeric types and associated type labels.  Within a
map, each type label must be unique.  It can be assigned to only one
numeric type.  Not all numeric types need have a type label assigned.

There can be mutliple label maps defined for a single type category,
e.g. atom types.  There is a default label map which has no mapID.
Additional label maps each have a mapID, which is a string containing
only alphanumeric characters and underscores.

Valid type labels can contain any alphanumeric character, but cannot
start with a number.  They can also contain any of these characters:
square brackets "[" and "]", dash "-", underscore "_", JAKE what other
chars are allowed?  Note that type labels cannot contain the symbols
'#' or '*'.

As explained below, for certain command in LAMMPS, it is useful to
define type labels so that the atom types that make up a given bond,
angle, dihedral, or improper can be deduced from the type label.  A
standard way to do that is to define the type label for the bond by
enclosing the constituent atom types in square brackets.  E.g. define
a C-H bond with a type label "[C][H]".

There are two ways to define label maps.  One is via the
:doc:`labelmap <labelmap>` command, which has an optional *mapID*
keyword to allow creation of type labels in either a default map or an
additional map with a mapID.  The other is via the :doc:`read_data
<read_data>` command.  A data file can have sections such as Atom Type
Labels, Bond Type Lables, etc, which associate type labels with
numeric types.  Only default label maps can be defined in this manner.

If defined, default label maps can be written out to data files by the
:doc:`write_data <write_data>` command.  They are also written to and
read from restart files, by the :doc:`write_restart <write_restart>`
and :doc:`read_restart <read_restart>` commands.  Label maps with
mapIDs cannot be written to either data or restart files by these
commands.

----------

Use of type labels in LAMMPS input or output
""""""""""""""""""""""""""""""""""""""""""""

Any LAMMPS input script command which takes a numeric type as an
argument, can use the associated type label instead, with the optional
mapID prepended, followed by a double colon "::".  If a type label is
not defined for a particular numeric type, only its numeric type can
be used.

The first example uses a default label map for bond types.  The second
uses a label map with mapID = Map2.

.. code-block:: LAMMPS

   bond_coeff 2 80.0 1.2               # numeric type
   labelmap bond 2 [C][H]
   bond_coeff [C][H] 80.0 1.2          # type label

or

.. code-block:: LAMMPS

   bond_coeff 2 80.0 1.2               # numeric type
   labelmap bond 2 [C][H] mapID Map2
   bond_coeff Map2::[C][H] 80.0 1.2    # type label

Support for type labels is a work-in-progress for LAMMPS as of
Nov 2021.  If an input script command allows substituting for a
numeric type argument with a type label, it will be mentioned on that
command's doc page.  If not yet supported, any input script command
that requires a numeric type can instead use a variable with a
labelmap function that converts a type label to a numeric type, as in
the last example above.  See the :doc:`variable <variable>` command
for details.

For example, here is how the bond_coeff command could be used with a
type label if it did not yet support them, either with an explicit
variable command or an implicit variable used in the bond_coeff
command.

.. code-block:: LAMMPS

   labelmap bond 2 [C][H]
   variable bond2 equal blabel([C][H])
   bond_coeff ${bond2} 80.0 1.2

.. code-block:: LAMMPS

   labelmap bond 2 [C][H]
   bond_coeff $(blabel([C][H])) 80.0 1.2

Support for output of type labels in dump files will be added to
LAMMPS soon (as of Nov 2021).

----------

Commands that can use label types in interesting ways
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Provide a few details on these to spark user's interest?  If not yet
implemented, just say this is planned, as of Nov 2021.

fix bond/react
pair_style kim
other?
