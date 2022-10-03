Type labels
===========

.. versionadded:: 15Sep2022

Each atom in LAMMPS has an associated numeric atom type. Similarly,
each bond, angle, dihedral, and improper is assigned a bond type,
angle type, and so on.  The primary use of these types is to map
potential (force field) parameters to the interactions of the atom,
bond, angle, dihedral, and improper.

By default, type values are entered as integers from 1 to Ntypes
wherever they appear in LAMMPS input or output files.  The total number
Ntypes for each interaction is "locked in" when the simulation box
is created.

A recent addition to LAMMPS is the option to use strings - referred
to as type labels - as an alternative.  Using type labels instead of
numeric types can be advantageous in various scenarios.  For example,
type labels can make inputs more readable and generic (i.e. usable through
the :doc:`include command <include>` for different systems with different
numerical values assigned to types.  This generality also applies to
other inputs like data files read by :doc:`read_data <read_data>` or
molecule template files read by the :doc:`molecule <molecule>`
command.  See below for a list of other commands that can use
type labels in different ways.

LAMMPS will *internally* continue to use numeric types, which means
that many previous restrictions still apply.  For example, the total
number of types is locked in when creating the simulation box, and
potential parameters for each type must be provided even if not used
by any interactions.

A collection of type labels for all type-kinds (atom types, bond types,
etc.) is stored as a "label map" which is simply a list of numeric types
and their associated type labels.  Within a type-kind, each type label
must be unique.  It can be assigned to only one numeric type.  To read
and write type labels to data files for a given type-kind, *all*
associated numeric types need have a type label assigned.  Partial
maps can be saved with the :doc:`labelmap write <labelmap>` command
and read back with the :doc:`include <include>` command.

Valid type labels can contain most ASCII characters, but cannot start
with a number, a '#', or a '*'.  Also, labels must not contain whitespace
characters.  When using the :doc:`labelmap command <labelmap>` in the
LAMMPS input, if certain characters appear in the type label, such as
the single (') or double (") quote or the '#' character, the label
must be put in either double, single, or triple (""") quotes.  Triple
quotes allow for the most generic type label strings, but they require
to have a leading and trailing blank space.  When defining type labels
the blanks will be ignored. Example:

.. code-block:: LAMMPS

   labelmap angle 1 """ C1'-C2"-C3# """

This command will map the string ```C1'-C2"-C3#``` to the angle type 1.

There are two ways to define label maps.  One is via the :doc:`labelmap
<labelmap>` command.  The other is via the :doc:`read_data <read_data>`
command.  A data file can have sections such as *Atom Type Labels*, *Bond
Type Labels*, etc., which assign type labels to numeric types.  The
label map can be written out to data files by the :doc:`write_data
<write_data>` command.  This map is also written to and read from
restart files, by the :doc:`write_restart <write_restart>` and
:doc:`read_restart <read_restart>` commands.

----------

Use of type labels in LAMMPS input or output
""""""""""""""""""""""""""""""""""""""""""""

Many LAMMPS input script commands that take a numeric type as an
argument can use the associated type label instead.  If a type label
is not defined for a particular numeric type, only its numeric type
can be used.

This example assigns labels to the atom types, and then uses the type
labels to redefine the pair coefficients.

.. code-block:: LAMMPS

   pair_coeff 1 2 1.0 1.0              # numeric types
   labelmap atom 1 C 2 H
   pair_coeff C H 1.0 1.0              # type labels

Adding support for type labels to various commands is an ongoing
project.  If an input script command (or a section in a file read by a
command) allows substituting a type label for a numeric type argument,
it will be explicitly mentioned in that command's documentation page.

As a temporary measure, input script commands can take advantage of
variables and how they can be expanded during processing of the input.
The variables can use functions that will translate type label strings
to their respective number as defined in the current label map.  See the
:doc:`variable <variable>` command for details.

For example, here is how the pair_coeff command could be used with
type labels if it did not yet support them, either with an explicit
variable command or an implicit variable used in the pair_coeff
command.

.. code-block:: LAMMPS

   labelmap atom 1 C 2 H
   variable atom1 equal label2type(atom,C)
   variable atom2 equal label2type(atom,H)
   pair_coeff ${atom1} ${atom2} 1.0 1.0

.. code-block:: LAMMPS

   labelmap atom 1 C 2 H
   pair_coeff $(label2type(atom,C)) $(label2type(atom,H)) 80.0 1.2

----------

Commands that can use label types
"""""""""""""""""""""""""""""""""

Any workflow that involves reading multiple data files, molecule
templates or a combination of the two can be streamlined by using type
labels instead of numeric types, because types are automatically synced
between the files.  The creation of simulation-ready reaction templates
for :doc:`fix bond/react <fix_bond_react>` is much simpler when using
type labels, and results in templates that can be used without
modification in multiple simulations or different systems.
