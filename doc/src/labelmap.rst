.. index:: labelmap

labelmap command
==================

Syntax
""""""

.. code-block:: LAMMPS

   labelmap option args

* *option* = *atom* or *bond* or *angle* or *dihedral* or *improper* or *clear* or *write* or *check_labels*

  .. parsed-literal::

     *clear* = no args
     *write* arg = filename
     *atom* or *bond* or *angle* or *dihedral* or *improper*
       args = list of one or more numeric-type/type-label pairs
     *check_labels* value = string
       string = string containing any of the following characters: 'b', 'a', 'd', or 'i'

Examples
""""""""

.. code-block:: LAMMPS

   labelmap atom 1 c1 2 hc 3 cp 4 nt
   labelmap atom 3 carbon 4 'c3"' 5 "c1'" 6 "c#"
   labelmap atom $(label2type(atom,carbon)) C  # change type label from 'carbon' to 'C'
   labelmap clear
   labelmap write mymap.include
   labelmap bond 1 carbonyl 2 nitrile 3 """ c1'-c2" """

Description
"""""""""""

.. versionadded:: 15Sep2022

Define alphanumeric type labels to associate with one or more numeric
atom, bond, angle, dihedral or improper types.  A collection of type
labels for all atom types, bond types, etc. is stored as a label map.

The label map can also be defined by the :doc:`read_data <read_data>`
command when it reads these sections in a data file: Atom Type Labels,
Bond Type Labels, etc.  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for a general discussion of how type
labels can be used.  See :ref:`(Gissinger) <Typelabel1>` for a
discussion of the type label implementation in LAMMPS and its uses.

Valid type labels can contain any alphanumeric character, but must not
start with a number, a '#', or a '*' character.  They can contain other
standard ASCII characters such as angular or square brackets '<' and '>'
or '[' and ']', parenthesis '(' and ')', dash '-', underscore '_', plus
'+' and equals '=' signs and more.  They must not contain blanks or any
other whitespace.  Note that type labels must be put in single or double
quotation marks if they contain the '#' character or if they contain a
double (") or single quotation mark (').  If the label contains both a
single and a double quotation mark, then triple quotation (""") must be
used.  When enclosing a type label with quotation marks, the LAMMPS
input parser may require adding leading or trailing blanks around the
type label so it can identify the enclosing quotation marks.  Those
blanks will be removed when defining the label.

A *labelmap* command can only modify the label map for one type-kind
(atom types, bond types, etc).  Any number of numeric-type/type-label
pairs may follow.  If a type label already exists for the same numeric
type, it will be overwritten.  Type labels must be unique; assigning the
same type label to multiple numeric types within the same type-kind is
not allowed.  When reading and writing data files, it is required that
there is a label defined for *every* numeric type within a given
type-kind in order to write out the type label section for that
type-kind.

The *clear* option resets the label map and thus discards all previous
settings.

The *write* option takes a filename as argument and writes the current
label mappings to a file as a sequence of *labelmap* commands, so the
file can be copied into a new LAMMPS input file or read in using the
:doc:`include <include>` command.

.. versionadded:: 30Mar2026

The *check_labels* keyword provides a warning if the type label of a
bond, angle, dihedral, or improper defined in the simulation is not
consistent with the atom types of its constituent atoms.  This
consistency check is performed only once, when the simulation is
initialized by the first :doc:`run <run>` or :doc:`minimize <minimize>`
command after this *labelmap* command.  The *check_labels* value is a
single string that should contain one or more of the characters 'b',
'a', 'd', and 'i', which correspond to bonds, angles, dihedrals, and
impropers, respectively.  For example, the keyword/value pair
'check_labels badi' will check all the type labels of all higher-order
interactions, while 'check_labels adi' will only check type labels for
angles, dihedrals, and impropers.  The *check_labels* keyword requires a
specific :doc:`type label<Howto_type_labels>` format to infer the types
of higher-order interactions.  Bond, angle, dihedral, and improper type
labels must contain their constituent atom types delimited by hyphens,
e.g., 'c2-c2-c2-n' for a dihedral that contains three atoms of type 'c2'
and one atom of 'n'.  If the constituent atoms do not have these atom
types in the proper order, a warning will be generated when using this
*check_types* keyword.  Interactions that have been disabled, e.g., via
:doc:`fix_shake <fix_shake>`, will not be checked.

----------

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

Label maps are considered experimental when using the KOKKOS package.

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`write_data <write_data>`,
:doc:`molecule <molecule>`, :doc:`fix bond/react <fix_bond_react>`

Default
"""""""

none

-----------

.. _Typelabel1:

**(Gissinger)** J. R. Gissinger, I. Nikiforov, Y. Afshar, B. Waters, M. Choi, D. S. Karls, A. Stukowski, W. Im, H. Heinz, A. Kohlmeyer, and E. B. Tadmor, J Phys Chem B, 128, 3282-3297 (2024).
