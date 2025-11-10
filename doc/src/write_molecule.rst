.. index:: write_molecule

write_molecule command
======================

Syntax
""""""

.. code-block:: LAMMPS

   write_molecule mol-ID file

* mol-ID = ID of the molecule template to be written
* file = name of file to write the molecule template to

Examples
""""""""

.. code-block:: LAMMPS

   write_molecule mol1 molecule1.mol
   write_molecule mol1 molecule1.json
   write_molecule twomols template_set%.mol

Description
"""""""""""

.. versionadded:: TBD

Write the data from a :doc:`molecule template<molecule>` to a molecule file.

The molecule file format is determined by the file name: if the file
name ends in ``.json`` the file will be written in `JSON format
<https://www.json.org/>`_, otherwise the file is written in the native
LAMMPS molecule file format.

When the molecule template contains multiple molecules, as defined by a
:doc:`molecule command <molecule>` with multiple molecule files, the
filename *must* contain a '%' character.  That '%' character will be
replaced by the molecule number (starting from 1) and each molecule is
written to a separate file.

----------

Output Format
"""""""""""""

The output format follows the description of the molecule file format in
the :doc:`molecule command documentation <molecule>`.  When parsing
molecule files there are a few requirement about the order of sections
for the native format, but generally the order of list entries can be
chosen freely.  On writing, however, a specific ordering is enforced
following the lists given below.  All sections with per-atom data are
sorted by the atom-ID value starting from 1.  The list of bonds is
ordered by the atom-ID of the first atom in the bond.  The lists of
angles, dihedrals, and impropers are ordered by the atom-ID of the
second atom in the definition.  For entries with the same atom-ID, the
order in which they were defined in the original input data is
maintained.

Native Format
^^^^^^^^^^^^^

The native format starts with a title line similar to the following where
the angular brackets are replaced with the actual data:

.. parsed-literal::

   # MOLECULE <molecule-ID>, unit = <units setting>, set <set-ID> of <nsets>, <original title>

This followed by a header section with the following keywords in the
listed order:

- `atoms` (required keyword)
- `bonds` (if non-zero)
- `angles` (if non-zero)
- `dihedrals` (if non-zero)
- `impropers` (if non-zero)
- `fragments` (if non-zero)
- `mass` (if present in original input)
- `body` (if present in original input)
- `com` (if present in original input)
- `inertia` (if present in original input)

The header section is followed by individual data sections in the order
given below if the corresponding data was provided in the original input
or explicitly added later:

- `Coords`
- `Types`
- `Molecules`
- `Fragments`
- `Charges`
- `Diameters`
- `Dipoles`
- `Masses`
- `Bonds`
- `Angles`
- `Dihedrals`
- `Impropers`
- `Special Bond Counts`
- `Special Bonds`
- `Shake Flags`
- `Shake Atoms`
- `Shake Bond Types`
- `Body Integers`
- `Body Doubles`

JSON Format
^^^^^^^^^^^

The JSON output is written as a text file in strict `JSON format
<https://www.json.org/>`_ with an indentation level of 2 and following
the layout given in the `JSON schema file <https://json-schema.org/>`_
available for download at
https://download.lammps.org/json/molecule-schema.json

Unlike for the native format output, there is no need to provide any
counts for the number of entries for keys containing lists of settings
with a specified format, since those can be directly inferred from the
data structures after parsing a JSON file. The top level keys are
ordered as follows:

- `application`
- `format`
- `revision`
- `schema`
- `title`
- `units`
- `masstotal` (if present in original input)
- `com` (if present in original input)
- `inertia` (if present in original input)
- `coords` (if present in molecule data)
- `types`
- `molecules` (if present in molecule data)
- `fragments` (if present in molecule data)
- `charges` (if present in molecule data)
- `diameters` (if present in molecule data)
- `dipoles` (if present in molecule data)
- `masses` (if present in molecule data)
- `bonds` (if present in molecule data)
- `angles` (if present in molecule data)
- `dihedrals` (if present in molecule data)
- `impropers` (if present in molecule data)
- `special` (if present in original input)
- `shake` (if present in molecule data)
- `body` (if present in molecule data)

----------

Restrictions
""""""""""""

None

Related commands
""""""""""""""""

:doc:`molecule <molecule>`

Defaults
""""""""

None
