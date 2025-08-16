.. index:: molecule

molecule command
================

Syntax
""""""

.. code-block:: LAMMPS

   molecule ID file1 keyword values ... file2 keyword values ... fileN ...

* ID = user-assigned name for the molecule template
* file1,file2,... = names of files containing molecule descriptions
* zero or more keyword/value pairs may be appended after each file
* keyword = *offset* or *toff* or *boff* or *aoff* or *doff* or *ioff* or *scale*

  .. parsed-literal::

       *offset* values = Toff Boff Aoff Doff Ioff
         Toff = offset to add to atom types
         Boff = offset to add to bond types
         Aoff = offset to add to angle types
         Doff = offset to add to dihedral types
         Ioff = offset to add to improper types
       *toff* value = Toff
         Toff = offset to add to atom types
       *boff* value = Boff
         Boff = offset to add to bond types
       *aoff* value = Aoff
         Aoff = offset to add to angle types
       *doff* value = Doff
         Doff = offset to add to dihedral types
       *ioff* value = Ioff
         Ioff = offset to add to improper types
       *scale* value = sfactor
         sfactor = scale factor to apply to the size, mass, and dipole of the molecule

Examples
""""""""

.. code-block:: LAMMPS

   molecule 1 mymol.txt
   molecule water tip3p.json
   molecule 1 co2.txt h2o.txt
   molecule CO2 co2.txt boff 3 aoff 2
   molecule 1 mymol.txt offset 6 9 18 23 14
   molecule objects file.1 scale 1.5 file.1 scale 2.0 file.2 scale 1.3

Description
"""""""""""

Define a molecule template that can be used as part of other LAMMPS
commands, typically to define a collection of particles as a bonded
molecule or a rigid body.  Commands that currently use molecule
templates include:

* :doc:`fix deposit <fix_deposit>`
* :doc:`fix pour <fix_pour>`
* :doc:`fix rigid/small <fix_rigid>`
* :doc:`fix shake <fix_shake>`
* :doc:`fix gcmc <fix_gcmc>`
* :doc:`fix bond/react <fix_bond_react>`
* :doc:`create_atoms <create_atoms>`
* :doc:`atom_style template <atom_style>`

The ID of a molecule template can only contain alphanumeric characters
and underscores, same as other IDs in LAMMPS.

A single template can contain multiple molecules, listed one per file.
Some of the commands listed above currently use only the first
molecule in the template, and will issue a warning if the template
contains multiple molecules.  The :doc:`atom_style template
<atom_style>` command allows multiple-molecule templates to define a
system with more than one templated molecule.

The molecule file can be either in a *native* format or in `JSON format
<https://www.json.org/>`_.  JSON formal filenames **must** have the
extension ".json".  Files with any other name will be assumed to be in
the "native" format.  The details of the two formats are described
below.  When referencing multiple molecule files in a single *molecule*
command, each of those files may be either format.

Each filename can be followed by optional keywords which are applied
only to the molecule in the file as used in this template.  This is to
make it easy to use the same molecule file in different molecule
templates or in different simulations.  You can specify the same file
multiple times with different optional keywords.

The *offset*, *toff*, *boff*, *aoff*, *doff*, *ioff* keywords
add the specified offset values to the atom types, bond types, angle
types, dihedral types, and/or improper types as they are read from the
molecule file.  E.g. if *toff* = 2, and the file uses atom types
1,2,3, then each created molecule will have atom types 3,4,5.  For the
*offset* keyword, all five offset values must be specified, but
individual values will be ignored if the molecule template does not
use that attribute (e.g. no bonds).

.. note::

   Offsets are **ignored** on lines using type labels, as the type
   labels will determine the actual types directly depending on the
   current :doc:`labelmap <labelmap>` settings.

The *scale* keyword scales the size of the molecule.  This can be useful
for modeling polydisperse granular rigid bodies.  The scale factor is
applied to each of these properties in the molecule file, if they are
defined: the individual particle coordinates (Coords or "coords"
section), the individual mass of each particle (Masses or "masses"
section), the individual diameters of each particle (Diameters or
"diameters" section), the per-atom dipoles (Dipoles or "dipoles"
section) the total mass of the molecule (header keyword = mass), the
center-of-mass of the molecule (header keyword = com), and the moments
of inertia of the molecule (header keyword = inertia).

.. note::

   The molecule command can be used to define molecules with bonds,
   angles, dihedrals, impropers, and special bond lists of neighbors
   within a molecular topology, so that you can later add the molecules
   to your simulation, via one or more of the commands listed above.
   Since this topology-related information requires that suitable
   storage is reserved when LAMMPS creates the simulation box (e.g. when
   using the :doc:`create_box <create_box>` command or the
   :doc:`read_data <read_data>` command) suitable space has to be
   reserved at that step so you do not overflow those pre-allocated data
   structures when adding molecules later.  Both the :doc:`create_box
   <create_box>` command and the :doc:`read_data <read_data>` command
   have "extra" options which ensure extra space is allocated for
   storing topology info for molecules that are added later.  This
   feature is *not* available for the :doc:`read_restart command
   <read_restart>`, thus binary restart files need to be converted
   to data files first.

----------

Format of a native molecule file
""""""""""""""""""""""""""""""""

The format of an "native" individual molecule file looks similar but is
*different* from that of a data file read by the :doc:`read_data
<read_data>` commands.  Here is a simple example for a TIP3P water
molecule:

.. code-block::

   # Water molecule. TIP3P geometry
   # header section:
   3 atoms
   2 bonds
   1 angles

   # body section:
   Coords

   1    0.00000  -0.06556   0.00000
   2    0.75695   0.52032   0.00000
   3   -0.75695   0.52032   0.00000

   Types

   1        1   # O
   2        2   # H
   3        2   # H

   Charges

   1       -0.834
   2        0.417
   3        0.417

   Bonds

   1   1      1      2
   2   1      1      3

   Angles

   1   1      2      1      3

A molecule file has a header and a body.  The header appears first.  The
first line of the header and thus of the molecule file is *always*
skipped; it typically contains a description of the file or a comment
from the software that created the file.

Then lines are read one line at a time.  Lines can have a trailing
comment starting with '#' that is ignored.  There *must* be at least one
blank between any valid content and the comment.  If the line is blank
(i.e. contains only white-space after comments are deleted), it is
skipped.  If the line contains a header keyword, the corresponding
value(s) is/are read from the line.  A line that is *not* blank and does
*not* contains a header keyword begins the body of the file.

The body of the file contains zero or more sections.  The first line
of a section has only a keyword.  The next line is skipped.  The
remaining lines of the section contain values.  The number of lines
depends on the section keyword as described below.  Zero or more blank
lines can be used between sections.  Sections can appear in any order,
with a few exceptions as noted below.

These are the recognized header keywords.  Header lines can come in
any order.  The numeric value(s) are read from the beginning of the
line.  The keyword should appear at the end of the line.  All these
settings have default values, as explained below.  A line need only
appear if the value(s) are different than the default, except when
defining a *body* particle, which requires setting the number of
*atoms* to 1, and setting the *inertia* in a specific section (see below).

.. list-table::
      :header-rows: 1
      :widths: 21 12 47 20

      * - Number(s)
        - Keyword
        - Meaning
        - Default Value
      * - N
        - atoms
        - # of atoms N in molecule
        - keyword is *required*
      * - Nb
        - bonds
        - # of bonds Nb in molecule
        - 0
      * - Na
        - angles
        - # of angles Na in molecule
        - 0
      * - Nd
        - dihedrals
        - # of dihedrals Nd in molecule
        - 0
      * - Ni
        - impropers
        - # of impropers Ni in molecule
        - 0
      * - Nf
        - fragments
        - # of fragments Nf in molecule
        - 0
      * - Ninteger Ndouble
        - body
        - # of integer and floating-point values in :doc:`body particle <Howto_body>`
        - 0 0
      * - Mtotal
        - mass
        - total mass of molecule
        - computed
      * - Xc Yc Zc
        - com
        - coordinates of center-of-mass of molecule
        - computed
      * - Ixx Iyy Izz Ixy Ixz Iyz
        - inertia
        - 6 components of inertia tensor of molecule
        - computed

For *mass*, *com*, and *inertia*, the default is for LAMMPS to calculate
this quantity itself if needed, assuming the molecules consist of a set
of point particles or finite-size particles (with a non-zero diameter)
that do **not** overlap.  If finite-size particles in the molecule
**do** overlap, LAMMPS will not account for the overlap effects when
calculating any of these 3 quantities, so you should pre-compute them
yourself and list the values in the file.

The mass and center-of-mass coordinates (Xc,Yc,Zc) are
self-explanatory.  The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz)
should be the values consistent with the current orientation of the
rigid body around its center of mass.  The values are with respect to
the simulation box XYZ axes, not with respect to the principal axes of
the rigid body itself.  LAMMPS performs the latter calculation
internally.

These are the allowed section keywords for the body of the file.

* *Coords, Types, Molecules, Fragments, Charges, Diameters, Dipoles, Masses* = atom-property sections
* *Bonds, Angles, Dihedrals, Impropers* = molecular topology sections
* *Special Bond Counts, Special Bonds* = special neighbor info
* *Shake Flags, Shake Atoms, Shake Bond Types* = SHAKE info
* *Body Integers, Body Doubles* = body-property sections

For the Types, Bonds, Angles, Dihedrals, and Impropers sections, each
atom/bond/angle/etc type can be specified either as a number (numeric
type) or as an alphanumeric type label.  The latter is only allowed if
type labels have been defined, either by the :doc:`labelmap
<labelmap>` command or in data files read by the :doc:`read_data
<read_data>` command which have sections for Atom Type Labels, Bond
Type Labels, Angle Type Labels, etc.  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for the allowed syntax of type labels
and a general discussion of how type labels can be used.
When using type labels, any values specified as *offset* are ignored.

If a Bonds section is specified then the Special Bond Counts and
Special Bonds sections can also be used, if desired, to explicitly
list the 1-2, 1-3, 1-4 neighbors within the molecule topology (see
details below).  This is optional since if these sections are not
included, LAMMPS will auto-generate this information.  Note that
LAMMPS uses this info to properly exclude or weight bonded pairwise
interactions between bonded atoms.  See the :doc:`special_bonds
<special_bonds>` command for more details.  One reason to list the
special bond info explicitly is for the :doc:`thermalized Drude
oscillator model <Howto_drude>` which treats the bonds between nuclear
cores and Drude electrons in a different manner.

.. note::

   Whether a section is required depends on how the molecule template
   is used by other LAMMPS commands.  For example, to add a molecule
   via the :doc:`fix deposit <fix_deposit>` command, the Coords and
   Types sections are required.  To add a rigid body via the :doc:`fix
   pour <fix_pour>` command, the Bonds (Angles, etc) sections are not
   required, since the molecule will be treated as a rigid body.  Some
   sections are optional.  For example, the :doc:`fix pour <fix_pour>`
   command can be used to add "molecules" which are clusters of
   finite-size granular particles.  If the Diameters section is not
   specified, each particle in the molecule will have a default
   diameter of 1.0.  See the doc pages for LAMMPS commands that use
   molecule templates for more details.

Each section is listed below in alphabetic order.  The format of each
section is described including the number of lines it must contain and
rules (if any) for whether it can appear in the data file.  For per-
atom sections, entries should be numbered from 1 to Natoms (where
Natoms is the number of atoms in the template), indicating which atom
(or bond, etc) the entry applies to.  Per-atom sections need to
include a setting for every atom, but the atoms can be listed in any
order.

----------

*Coords* section:

* one line per atom
* line syntax: ID x y z
* x,y,z = coordinate of atom

----------

*Types* section:

* one line per atom
* line syntax: ID type
* type = atom type of atom (1-Natomtype, or type label)

----------

*Molecules* section:

* one line per atom
* line syntax: ID molecule-ID
* molecule-ID = molecule ID of atom

----------

*Fragments* section:

* one line per fragment
* line syntax: ID a b c d ...
* a,b,c,d,... = IDs of atoms in fragment

The ID of a fragment can only contain alphanumeric characters and
underscores.  The atom IDs should be values from 1 to Natoms, where
Natoms = # of atoms in the molecule.

----------

*Charges* section:

* one line per atom
* line syntax: ID q
* q = charge on atom

This section is only allowed for :doc:`atom styles <atom_style>` that
support charge.  If this section is not included, the default charge
on each atom in the molecule is 0.0.

----------

*Diameters* section:

* one line per atom
* line syntax: ID diam
* diam = diameter of atom

This section is only allowed for :doc:`atom styles <atom_style>` that
support finite-size spherical particles, e.g. atom_style sphere.  If
not listed, the default diameter of each atom in the molecule is 1.0.

----------

.. versionadded:: 7Feb2024

*Dipoles* section:

* one line per atom
* line syntax: ID mux muy muz
* mux,muy,muz = x-, y-, and z-component of point dipole vector of atom

This section is only allowed for :doc:`atom styles <atom_style>` that
support particles with point dipoles, e.g. atom_style dipole.  If not
listed, the default dipole component of each atom in the molecule is set
to 0.0.

----------

*Masses* section:

* one line per atom
* line syntax: ID mass
* mass = mass of atom

This section is only allowed for :doc:`atom styles <atom_style>` that
support per-atom mass, as opposed to per-type mass.  See the
:doc:`mass <mass>` command for details.  If this section is not
included, the default mass for each atom is derived from its volume
(see Diameters section) and a default density of 1.0, in
:doc:`units <units>` of mass/volume.

----------

*Bonds* section:

* one line per bond
* line syntax: ID type atom1 atom2
* type = bond type (1-Nbondtype, or type label)
* atom1,atom2 = IDs of atoms in bond

The IDs for the two atoms in each bond should be values
from 1 to Natoms, where Natoms = # of atoms in the molecule.

----------

*Angles* section:

* one line per angle
* line syntax: ID type atom1 atom2 atom3
* type = angle type (1-Nangletype, or type label)
* atom1,atom2,atom3 = IDs of atoms in angle

The IDs for the three atoms in each angle should be values from 1 to
Natoms, where Natoms = # of atoms in the molecule.  The three atoms are
ordered linearly within the angle.  Thus the central atom (around
which the angle is computed) is the atom2 in the list.

----------

*Dihedrals* section:

* one line per dihedral
* line syntax: ID type atom1 atom2 atom3 atom4
* type = dihedral type (1-Ndihedraltype, or type label)
* atom1,atom2,atom3,atom4 = IDs of atoms in dihedral

The IDs for the four atoms in each dihedral should be values from 1 to
Natoms, where Natoms = # of atoms in the molecule.  The 4 atoms are
ordered linearly within the dihedral.

----------

*Impropers* section:

* one line per improper
* line syntax: ID type atom1 atom2 atom3 atom4
* type = improper type (1-Nimpropertype, or type label)
* atom1,atom2,atom3,atom4 = IDs of atoms in improper

The IDs for the four atoms in each improper should be values from 1 to
Natoms, where Natoms = # of atoms in the molecule.  The ordering of
the 4 atoms determines the definition of the improper angle used in
the formula for the defined :doc:`improper style <improper_style>`.  See
the doc pages for individual styles for details.

----------

*Special Bond Counts* section:

* one line per atom
* line syntax: ID N1 N2 N3
* N1 = # of 1-2 bonds
* N2 = # of 1-3 bonds
* N3 = # of 1-4 bonds

N1, N2, N3 are the number of 1-2, 1-3, 1-4 neighbors respectively of
this atom within the topology of the molecule.  See the
:doc:`special_bonds <special_bonds>` page for more discussion of
1-2, 1-3, 1-4 neighbors.  If this section appears, the Special Bonds
section must also appear.

As explained above, LAMMPS will auto-generate this information if this
section is not specified.  If specified, this section will
override what would be auto-generated.

----------

*Special Bonds* section:

* one line per atom
* line syntax: ID a b c d ...
* a,b,c,d,... = IDs of atoms in N1+N2+N3 special bonds

A, b, c, d, etc are the IDs of the n1+n2+n3 atoms that are 1-2, 1-3,
1-4 neighbors of this atom.  The IDs should be values from 1 to
Natoms, where Natoms = # of atoms in the molecule.  The first N1
values should be the 1-2 neighbors, the next N2 should be the 1-3
neighbors, the last N3 should be the 1-4 neighbors.  No atom ID should
appear more than once.  See the :doc:`special_bonds <special_bonds>` doc
page for more discussion of 1-2, 1-3, 1-4 neighbors.  If this section
appears, the Special Bond Counts section must also appear.

As explained above, LAMMPS will auto-generate this information if this
section is not specified.  If specified, this section will override
what would be auto-generated.

----------

*Shake Flags* section:

* one line per atom
* line syntax: ID flag
* flag = 0,1,2,3,4

This section is only needed when molecules created using the template
will be constrained by SHAKE via the "fix shake" command.  The other
two Shake sections must also appear in the file, following this one.

The meaning of the flag for each atom is as follows.  See the :doc:`fix shake <fix_shake>` page for a further description of SHAKE
clusters.

* 0 = not part of a SHAKE cluster
* 1 = part of a SHAKE angle cluster (two bonds and the angle they form)
* 2 = part of a 2-atom SHAKE cluster with a single bond
* 3 = part of a 3-atom SHAKE cluster with two bonds
* 4 = part of a 4-atom SHAKE cluster with three bonds

----------

*Shake Atoms* section:

* one line per atom
* line syntax: ID a b c d
* a,b,c,d = IDs of atoms in cluster

This section is only needed when molecules created using the template
will be constrained by SHAKE via the "fix shake" command.  The other
two Shake sections must also appear in the file.

The a,b,c,d values are atom IDs (from 1 to Natoms) for all the atoms
in the SHAKE cluster that this atom belongs to.  The number of values
that must appear is determined by the shake flag for the atom (see the
Shake Flags section above).  All atoms in a particular cluster should
list their a,b,c,d values identically.

If flag = 0, no a,b,c,d values are listed on the line, just the
(ignored) ID.

If flag = 1, a,b,c are listed, where a = ID of central atom in the
angle, and b,c the other two atoms in the angle.

If flag = 2, a,b are listed, where a = ID of atom in bond with the
lowest ID, and b = ID of atom in bond with the highest ID.

If flag = 3, a,b,c are listed, where a = ID of central atom,
and b,c = IDs of other two atoms bonded to the central atom.

If flag = 4, a,b,c,d are listed, where a = ID of central atom,
and b,c,d = IDs of other three atoms bonded to the central atom.

See the :doc:`fix shake <fix_shake>` page for a further description
of SHAKE clusters.

----------

*Shake Bond Types* section:

* one line per atom
* line syntax: ID a b c
* a,b,c = bond types (or angle type) of bonds (or angle) in cluster

This section is only needed when molecules created using the template
will be constrained by SHAKE via the "fix shake" command.  The other
two Shake sections must also appear in the file.

The a,b,c values are bond types for all bonds in the SHAKE cluster that
this atom belongs to.  Bond types may be either numbers (from 1 to Nbondtypes)
or bond type labels as defined by the :doc:`labelmap <labelmap>` command
or a "Bond Type Labels" section of a data file.


The number of values that must appear is determined by the shake flag
for the atom (see the Shake Flags section above).  All atoms in a
particular cluster should list their a,b,c values identically.

If flag = 0, no a,b,c values are listed on the line, just the
(ignored) ID.

If flag = 1, a,b,c are listed, where a = bondtype of the bond between
the central atom and the first non-central atom (value b in the Shake
Atoms section), b = bondtype of the bond between the central atom and
the second non-central atom (value c in the Shake Atoms section), and c
= the angle type (1 to Nangletypes, or angle type label) of the angle
between the three atoms.

If flag = 2, only a is listed, where a = bondtype of the bond between
the two atoms in the cluster.

If flag = 3, a,b are listed, where a = bondtype of the bond between
the central atom and the first non-central atom (value b in the Shake
Atoms section), and b = bondtype of the bond between the central atom
and the second non-central atom (value c in the Shake Atoms section).

If flag = 4, a,b,c are listed, where a = bondtype of the bond between
the central atom and the first non-central atom (value b in the Shake
Atoms section), b = bondtype of the bond between the central atom and
the second non-central atom (value c in the Shake Atoms section), and c
= bondtype of the bond between the central atom and the third
non-central atom (value d in the Shake Atoms section).

See the :doc:`fix shake <fix_shake>` page for a further description
of SHAKE clusters.

----------

*Body Integers* section:

* one line
* line syntax: N E F
* N = number of sub-particles or number or vertices
* E,F = number of edges and faces

This section is only needed when the molecule is a body particle. the other
Body section must also appear in the file.

The total number of values that must appear is determined by the body style, and
must be equal to the Ninteger value given in the *body* header.

For *nparticle* and *rounded/polygon*, only the number of sub-particles or
vertices N is required, and Ninteger should have a value of 1.

For *rounded/polyhedron*, the number of edges E and faces F is required, and
Ninteger should have a value of 3.

See the :doc:`Howto body <Howto_body>` page for a further description of
the file format.

----------

*Body Doubles* section:

* first line
* line syntax: Ixx Iyy Izz Ixy Ixz Iyz
* Ixx Iyy Izz Ixy Ixz Iyz = 6 components of inertia tensor of body particle
* one line per sub-particle or vertex
* line syntax: x y z
* x, y, z = coordinates of sub-particle or vertex
* one line per edge
* line syntax: N1 N2
* N1, N2 = vertex indices
* one line per face
* line syntax: N1 N2 N3 N4
* N1, N2, N3, N4 = vertex indices
* last line
* line syntax: diam
* diam = rounded diameter that surrounds each vertex

This section is only needed when the molecule is a body particle. the other
Body section must also appear in the file.

The total number of values that must appear is determined by the body style, and
must be equal to the Ndouble value given in the *body* header. The 6 moments of
inertia and the 3N coordinates of the sub-particles or vertices are required
for all body styles.

For *rounded/polygon*, the E = 6 + 3*N + 1 edges are automatically determined
from the vertices.

For *rounded/polyhedron*, the 2E vertex indices for the end points of the edges
and 4F vertex indices defining the faces are required.

See the :doc:`Howto body <Howto_body>` page for a further description of
the file format.

----------

Format of a JSON molecule file
""""""""""""""""""""""""""""""

.. versionadded:: 12Jun2025

The format of a JSON format individual molecule file must follow the
`JSON format <https://www.json.org/>`_, which evolved from the
JavaScript programming language as a programming-language-neutral data
interchange language.  The JSON syntax is independent of its content,
and thus the data in the file must follow suitable conventions to be
correctly processed.  LAMMPS provides a `JSON schema file
<https://json-schema.org/>`_ for JSON format molecule files in the
:ref:`tools/json folder <json>` to represent those conventions.  Using
the schema file any JSON format molecule files can be validated.  Please
note that the format requirement for JSON are very strict and the JSON
reader in LAMMPS does not accept files with extensions like comments.
Validating a particular JSON format molecule file against this schema
ensures that both, the JSON syntax requirement *and* the LAMMPS
conventions for molecule template files are followed.  LAMMPS should be
able to read and parse any JSON file that passes the schema check.  This
is a formal check only and thus it **cannot** check whether the file
contents are consistent or physically meaningful.

Here is a simple example for the same TIP3P water molecule from above in
JSON format and also using :doc:`type labels <labelmap>` instead of
numeric types:

.. code-block:: json

   {
       "application": "LAMMPS",
       "format": "molecule",
       "revision": 1,
       "title": "Water molecule. TIP3P geometry",
       "schema": "https://download.lammps.org/json/molecule-schema.json",
       "units": "real",
       "coords": {
           "format": ["atom-id", "x", "y", "z"],
           "data": [
               [1,  0.00000, -0.06556,  0.00000],
               [2,  0.75695,  0.52032,  0.00000],
               [3, -0.75695,  0.52032,  0.00000]
           ]
       },
       "types": {
           "format": ["atom-id", "type"],
           "data": [
               [1,  "OW"],
               [2,  "HO1"],
               [3,  "HO1"]
           ]
       },
       "charges": {
           "format": ["atom-id", "charge"],
           "data": [
               [1, -0.834],
               [2,  0.417],
               [3,  0.417]
           ]
       },
       "bonds": {
           "format": ["bond-type", "atom1", "atom2"],
           "data": [
               ["OW-HO1",  1,  2],
               ["OW-HO1",  1,  3]
           ]
       },
       "angles": {
           "format": ["angle-type", "atom1", "atom2", "atom3"],
           "data": [
               ["HO1-OW-HO1",  2,  1,  3]
           ]
       }
   }

Unlike with the native molecule file format, there are no header or body
sections, just a list of keywords with associated data.  JSON format
data is read, parsed, and stored in an internal dictionary data
structure in one step and thus the order of keywords is not relevant.

Data for keywords is either provided directly following the keyword or
as a *data block*.  A *data block* is a list that has to include two
keywords, "format" and "data", where the former lists keywords of the
properties that are stored in the columns of the "data" lists.  The
names and order of entries in the "format" list (and thus how the data
is interpreted) are currently fixed.

Since the length of the various lists can be easily obtained from the
internal data structure, several header keywords of the "native" molecule
file are not needed.  On the other hand, some additional keywords are
required to identify the conventions applied to the generic JSON file
format.  The structure of the data itself mostly follows what is used
for the "native" molecule file format.

.. list-table::
   :header-rows: 1

   * - Keyword
     - Argument(s)
     - Required
     - Description
   * - application
     - "LAMMPS"
     - yes
     - indicates a LAMMPS JSON file; files from other applications may be accepted in the future
   * - format
     - "molecule"
     - yes
     - indicates a molecule template file
   * - revision
     - an integer
     - yes
     - currently 1, to facility backward compatibility on changes to the conventions
   * - title
     - a string
     - no
     - information about the template which will echoed to the screen and log
   * - schema
     - URL as string
     - no
     - location of a JSON schema file for validating the molecule file format
   * - units
     - a string
     - no
     - indicates :doc:`units settings <units>` for this molecule template
   * - com
     - list with 3 doubles
     - no
     - overrides the auto-computed center-of-mass for the template
   * - masstotal
     - double
     - no
     - overrides the auto-computed total mass for the template
   * - inertia
     - list with 6 doubles
     - no
     - overrides the auto-computed moments of inertia
   * - coords
     - a data block
     - no
     - contains atom positions with the format "atom-id", "x", "y", "z" (same as Coords)
   * - types
     - a data block
     - yes
     - assigns atom types to atoms with the format "atom-id", "type" (same as Types)
   * - molecule
     - a data block
     - no
     - assigns molecule-IDs to atoms with the format "atom-id", "molecule-id" (same as Molecules)
   * - fragments
     - a data block
     - no
     - assigns atom-ids to fragment-IDs with the format "fragment-id", "atom-id-list" (same as Fragments)
   * - charges
     - a data block
     - no
     - assigns charges to atoms with the format "atom-id", "charge" (same as Charges)
   * - dipoles
     - a data block
     - no
     - assigns point dipoles to atoms with the format "atom-id", "mux", "muy", "muz" (same as Dipoles)
   * - diameters
     - a data block
     - no
     - assigns diameters to atoms with the format "atom-id", "diameter" (same as Diameters)
   * - masses
     - a data block
     - no
     - assigns per-atom masses to atoms with the format "atom-id", "mass" (same as Masses)
   * - bonds
     - a data block
     - no
     - defines bonds in the molecule template with the format "bond-type", "atom1", "atom2" (same as Bonds without bond-ID)
   * - angles
     - a data block
     - no
     - defines angles in the molecule template with the format "angle-type", "atom1", "atom2", "atom3" (same as Angles without angle-ID)
   * - dihedrals
     - a data block
     - no
     - defines dihedrals in the molecule template with the format "dihedral-type", "atom1", "atom2", "atom3", "atom4" (same as Dihedrals without dihedral-ID)
   * - impropers
     - a data block
     - no
     - defines impropers in the molecule template with the format "improper-type", "atom1", "atom2", "atom3", "atom4" (same as Impropers without improper-ID)
   * - shake
     - 3 JSON objects
     - no
     - contains the sub-sections "flags", "atoms", "bonds" described below
   * - special
     - 2 JSON objects
     - no
     - contains the sub-sections "counts" and "bonds" described below
   * - body
     - 2 JSON objects
     - no
     - contains the "integers" and "doubles" sub-sections with arrays with the same data as Body Integers and Body Doubles, respectively

The following table describes the sub-sections for the "special" entry from above:

.. list-table::
   :header-rows: 1

   * - Subsection
     - Argument(s)
     - Required
     - Description
   * - counts
     - a data block
     - yes
     - contains the counts of 1-2, 1-3, and 1-4 special neighbors with the format "atom-id", "n12", "n13", "n14" (same as Special Bond Counts)
   * - bonds
     - a data block
     - yes
     - contains the lists of special neighbors to atoms with the format "atom-id", "atom-id-list" (same as Special Bonds)

The following table describes the sub-sections for the "shake" entry from above:

.. list-table::
   :header-rows: 1

   * - Subsection
     - Argument(s)
     - Required
     - Description
   * - flags
     - a data block
     - yes
     - contains the counts shake flags for atoms with the format "atom-id", "flag" (same as Shake Flags)
   * - atoms
     - a data block
     - yes
     - contains the lists of shake cluster atom-ids for atoms with the format "atom-id", "atom-id-list" (same as Shake Atoms)
   * - bonds
     - a data block
     - yes
     - contains the lists of shake bond or angle types for atoms with the format "atom-id", "type-list" (same as Shake Bonds)

The "special" and "shake" sections are usually not needed, since the
data can be auto-generated as soon as the simulation box is defined.
Below is an example for what would have to be *added* to the example
JSON file above in case the molecule command needs to be issued earlier.

.. code-block:: json

   "special": {
       "counts": {
           "format": ["atom-id", "n12", "n13", "n14"],
           "data": [
               [1, 2, 0, 0],
               [2, 1, 1, 0],
               [3, 1, 1, 0]
           ]
       },
       "bonds": {
           "format": ["atom-id", "atom-id-list"],
           "data": [
               [1,  [2, 3]],
               [2,  [1, 3]],
               [3,  [1, 2]]
           ]
       }
   },
   "shake": {
       "flags": {
           "format": ["atom-id", "flag"],
           "data": [
               [1, 1],
               [2, 1],
               [3, 1]
           ]
       },
       "atoms": {
           "format": ["atom-id", "atom-id-list"],
           "data": [
               [1,  [1, 2, 3]],
               [2,  [1, 2, 3]],
               [3,  [1, 2, 3]]
           ]
       },
       "types": {
           "format": ["atom-id", "type-list"],
           "data": [
               [1,  ["OW-HO1",  "OW-HO1", "HO1-OW-HO1"]],
               [2,  ["OW-HO1",  "OW-HO1", "HO1-OW-HO1"]],
               [3,  ["OW-HO1",  "OW-HO1", "HO1-OW-HO1"]]
           ]
       }
   }


Below is a minimal example of a JSON format molecule template for a body
particle for :doc:`pair style body/nparticle
<pair_body_nparticle>`. Molecule templates for body particles must
contain only one atom:

.. code-block:: json

   {
       "application": "LAMMPS",
       "format": "molecule",
       "revision": 1,
       "title": "Square body for body/nparticles",
       "schema": "https://download.lammps.org/json/molecule-schema.json",
       "units": "real",
       "coords": {
           "format": ["atom-id", "x", "y", "z"],
           "data": [
               [1,  0.00000,  0.00000,  0.00000]
           ]
       },
       "types": {
           "format": ["atom-id", "type"],
           "data": [
               [1,  1]
           ]
       },
       "masses": {
           "format": ["atom-id", "mass"],
           "data": [
               [1,  1.0]
           ]
       },
       "body": {
           "integers": [4],
           "doubles": [
               1.0, 1.0, 4.0, 0.0, 0.0, 0.0,
               -0.70710678118654752440, -0.70710678118654752440, 0.0,
               -0.70710678118654752440,  0.70710678118654752440, 0.0,
               0.70710678118654752440,  0.70710678118654752440, 0.0,
               0.70710678118654752440, -0.70710678118654752440, 0.0
           ]
       }
   }

----------

Restrictions
""""""""""""

None

Related commands
""""""""""""""""

:doc:`fix deposit <fix_deposit>`, :doc:`fix pour <fix_pour>`,
:doc:`fix gcmc <fix_gcmc>`

Default
"""""""

The default keywords values are offset 0 0 0 0 0 and scale = 1.0.
