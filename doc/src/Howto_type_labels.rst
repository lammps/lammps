Type labels
===========

Each atom in LAMMPS has an associated atom type.  Likewise each bond,
angle, dihedral, improper enumerated in a data file read by the
:doc:`read_data <read_data>` command has a bond type, angle type, etc.

By default, type values are integers from 1 to Ntypes, wherever they
appear in LAMMPS input or output files.  However, LAMMPS also has the
option to use alphanumeric strings, called type labels, for these
values in its input and output.

Using labels instead of numeric types can be useful in various
scenarios.  For example, type labels can make inputs more readable and
compatible with other inputs (e.g., data files, molecule template
files read by the :doc:`molecule <molecule>` command, etc.),
particularly if a force field uses alphanumeric atom types. See a list
of other commands below which can use label types in interesting ways.

A collection of type labels for all type-kinds (atom types, bond
types, etc.) is stored as a "label map" which is simply a list of
numeric types and associated type labels.  Within a type-kind, each
type label must be unique.  It can be assigned to only one numeric
type. To read and write type labels to data files for a given
type-kind, all associated numeric types need have a type label
assigned.

Valid type labels can contain any alphanumeric character, but cannot
start with a number.  They can also contain standard characters, such
as square brackets "[" and "]", underscore "_", dash "-", plus "+" and
equals "=" signs.  Note that type labels cannot contain the comment
symbol '#'.

There are two ways to define label maps.  One is via the
:doc:`labelmap <labelmap>` command.  The other is via the
:doc:`read_data <read_data>` command.  A data file can have sections
such as Atom Type Labels, Bond Type Labels, etc, which associate type
labels with numeric types.  The label map can be written out to data
files by the :doc:`write_data <write_data>` command.  This map is also
written to and read from restart files, by the :doc:`write_restart
<write_restart>` and :doc:`read_restart <read_restart>` commands.

----------

Use of type labels in LAMMPS input or output
""""""""""""""""""""""""""""""""""""""""""""

Any LAMMPS input script command which takes a numeric type as an
argument, can use the associated type label instead.  If a type label
is not defined for a particular numeric type, only its numeric type
can be used.

This example assigns labels to the atom types, and then uses the type
labels to redefine the pair coefficients.

.. code-block:: LAMMPS

   pair_coeff 1 2 1.0 1.0              # numeric types
   labelmap atom 1 C 2 H
   pair_coeff C H 1.0 1.0              # type labels

Support for type labels is a work-in-progress for LAMMPS as of
Nov 2021.  If an input script command allows substituting for a
numeric type argument with a type label, it will be mentioned on that
command's doc page.  If not yet supported, any input script command
that requires a numeric type can instead use a variable with a
labelmap function that converts a type label to a numeric type, as in
the last example above.  See the :doc:`variable <variable>` command
for details.

For example, here is how the pair_coeff command could be used with
type labels if it did not yet support them, either with an explicit
variable command or an implicit variable used in the pair_coeff
command.

.. code-block:: LAMMPS

   labelmap atom 1 C 2 H
   variable atom1 equal label(C)
   variable atom2 equal label(H)
   pair_coeff ${atom1} ${atom2} 1.0 1.0

.. code-block:: LAMMPS

   labelmap atom 1 C 2 H
   pair_coeff $(label(C)) $(label(H)) 80.0 1.2

Support for output of type labels in dump files will be added to
LAMMPS soon (as of Nov 2021).

----------

Commands that can use label types in interesting ways
"""""""""""""""""""""""""""""""""""""""""""""""""""""

As of Nov 2021, efforts are underway to utilize type labels in various
commands.

Any workflow that involves reading multiple data files, molecule
templates or a combination of the two will be greatly streamlined by
using type labels instead of numeric types, because types are
automatically synced between the files.  For example, the creation of
simulation-ready reaction templates for :doc:`fix bond/react <fix_bond_react>`
is much simpler when using type labels, and results in templates that
can be used without modification in new simulations.  Additional fix
bond/react features enabled by type labels are in progress, such as
using wildcards to further increase the portability of reaction
templates, as well as automatically inferring the types of newly
created bond, angle, etc. interactions.

LAMMPS type labels will be used in a planned extension of OpenKIM to
support bonded force fields (FFs) (such as CHARMM, AMBER, IFF, etc.).
Users will be able to use a bonded FF, packaged as an OpenKIM
Simulator Model (SM), using the `kim init` command.  The SM will
include all required interaction parameters (pair, bond, angle,
dihedral, improper) defined in terms of the standard atom types for
that FF.  Molecular configurations can then be specified within a
LAMMPS script or read in from a data file by using type labels that
match the atom types for that FF.
