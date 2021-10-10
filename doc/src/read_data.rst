.. index:: read_data

read_data command
=================

Syntax
""""""

.. code-block:: LAMMPS

   read_data file keyword args ...

* file = name of data file to read in
* zero or more keyword/arg pairs may be appended
* keyword = *add* or *offset* or *shift* or *extra/atom/types* or *extra/bond/types* or *extra/angle/types* or *extra/dihedral/types* or *extra/improper/types* or *extra/bond/per/atom* or *extra/angle/per/atom* or *extra/dihedral/per/atom* or *extra/improper/per/atom* or *group* or *nocoeff* or *fix*

  .. parsed-literal::

       *add* arg = *append* or *IDoffset* or *IDoffset MOLoffset* or *merge*
         append = add new atoms with atom IDs appended to current IDs
         IDoffset = add new atoms with atom IDs having IDoffset added
         MOLoffset = add new atoms with molecule IDs having MOLoffset added (only when molecule IDs are enabled)
         merge = add new atoms with their atom IDs (and molecule IDs) unchanged
       *offset* args = toff boff aoff doff ioff
         toff = offset to add to atom types
         boff = offset to add to bond types
         aoff = offset to add to angle types
         doff = offset to add to dihedral types
         ioff = offset to add to improper types
       *shift* args = Sx Sy Sz
         Sx,Sy,Sz = distance to shift atoms when adding to system (distance units)
       *extra/atom/types* arg = # of extra atom types
       *extra/bond/types* arg = # of extra bond types
       *extra/angle/types* arg = # of extra angle types
       *extra/dihedral/types* arg = # of extra dihedral types
       *extra/improper/types* arg = # of extra improper types
       *extra/bond/per/atom* arg = leave space for this many new bonds per atom
       *extra/angle/per/atom* arg = leave space for this many new angles per atom
       *extra/dihedral/per/atom* arg = leave space for this many new dihedrals per atom
       *extra/improper/per/atom* arg = leave space for this many new impropers per atom
       *extra/special/per/atom* arg = leave space for extra 1-2,1-3,1-4 interactions per atom
       *group* args = groupID
         groupID = add atoms in data file to this group
       *nocoeff* = ignore force field parameters
       *fix* args = fix-ID header-string section-string
         fix-ID = ID of fix to process header lines and sections of data file
         header-string = header lines containing this string will be passed to fix
         section-string = section names with this string will be passed to fix

Examples
""""""""

.. code-block:: LAMMPS

   read_data data.lj
   read_data ../run7/data.polymer.gz
   read_data data.protein fix mycmap crossterm CMAP
   read_data data.water add append offset 3 1 1 1 1 shift 0.0 0.0 50.0
   read_data data.water add merge 1 group solvent

Description
"""""""""""

Read in a data file containing information LAMMPS needs to run a
simulation.  The file can be ASCII text or a gzipped text file
(detected by a .gz suffix).  This is one of 3 ways to specify initial
atom coordinates; see the :doc:`read_restart <read_restart>` and
:doc:`create_atoms <create_atoms>` commands for alternative methods.
Also see the explanation of the :doc:`-restart command-line switch
<Run_options>` which can convert a restart file to a data file.

This command can be used multiple times to add new atoms and their
properties to an existing system by using the *add*, *offset*, and
*shift* keywords.  See more details below, which includes the use case
for the *extra* keywords.

The *group* keyword adds all the atoms in the data file to the
specified group-ID.  The group will be created if it does not already
exist.  This is useful if you are reading multiple data files and wish
to put sets of atoms into different groups so they can be operated on
later.  E.g. a group of added atoms can be moved to new positions via
the :doc:`displace_atoms <displace_atoms>` command.  Note that atoms
read from the data file are also always added to the "all" group.  The
:doc:`group <group>` command discusses atom groups, as used in LAMMPS.

The *nocoeff* keyword tells read_data to ignore force field parameters.
The various Coeff sections are still read and have to have the correct
number of lines, but they are not applied. This also allows to read a
data file without having any pair, bond, angle, dihedral or improper
styles defined, or to read a data file for a different force field.

The use of the *fix* keyword is discussed below.

----------

Reading multiple data files
"""""""""""""""""""""""""""

The read_data command can be used multiple times with the same or
different data files to build up a complex system from components
contained in individual data files.  For example one data file could
contain fluid in a confined domain; a second could contain wall atoms,
and the second file could be read a third time to create a wall on the
other side of the fluid.  The third set of atoms could be rotated to
an opposing direction using the :doc:`displace_atoms <displace_atoms>`
command, after the third read_data command is used.

The *add*, *offset*, *shift*, *extra*, and *group* keywords are
useful in this context.

If a simulation box does not yet exist, the *add* keyword cannot be
used; the read_data command is being used for the first time.  If a
simulation box does exist, due to using the :doc:`create_box
<create_box>` command, or a previous read_data command, then the *add*
keyword must be used.

.. note::

   The simulation box size (xlo to xhi, ylo to yhi, zlo to zhi) in
   the new data file will be merged with the existing simulation box to
   create a large enough box in each dimension to contain both the
   existing and new atoms.  Each box dimension never shrinks due to this
   merge operation, it only stays the same or grows. Care must be used if
   you are growing the existing simulation box in a periodic dimension.
   If there are existing atoms with bonds that straddle that periodic
   boundary, then the atoms may become far apart if the box size grows.
   This will separate the atoms in the bond, which can lead to "lost"
   bond atoms or bad dynamics.

The three choices for the *add* argument affect how the atom IDs and
molecule IDs of atoms in the data file are treated.  If *append* is
specified, atoms in the data file are added to the current system,
with their atom IDs reset so that an atom-ID = M in the data file
becomes atom-ID = N+M, where N is the largest atom ID in the current
system.  This rule is applied to all occurrences of atom IDs in the
data file, e.g. in the Velocity or Bonds section. This is also done
for molecule IDs, if the atom style does support molecule IDs or
they are enabled via fix property/atom. If *IDoffset* is specified,
then *IDoffset* is a numeric value is given, e.g. 1000, so that an
atom-ID = M in the data file becomes atom-ID = 1000+M. For systems
with enabled molecule IDs, another numerical argument *MOLoffset*
is required representing the equivalent offset for molecule IDs.
If *merge* is specified, the data file atoms
are added to the current system without changing their IDs.  They are
assumed to merge (without duplication) with the currently defined
atoms.  It is up to you to insure there are no multiply defined atom
IDs, as LAMMPS only performs an incomplete check that this is the case
by insuring the resulting max atom-ID >= the number of atoms. For
molecule IDs, there is no check done at all.

The *offset* and *shift* keywords can only be used if the *add*
keyword is also specified.

The *offset* keyword adds the specified offset values to the atom
types, bond types, angle types, dihedral types, and improper types as
they are read from the data file.  E.g. if *toff* = 2, and the file
uses atom types 1,2,3, then the added atoms will have atom types
3,4,5.  These offsets apply to all occurrences of types in the data
file, e.g. for the Atoms or Masses or Pair Coeffs or Bond Coeffs
sections.  This makes it easy to use atoms and molecules and their
attributes from a data file in different simulations, where you want
their types (atom, bond, angle, etc) to be different depending on what
other types already exist.  All five offset values must be specified,
but individual values will be ignored if the data file does not use
that attribute (e.g. no bonds).

The *shift* keyword can be used to specify an (Sx, Sy, Sz)
displacement applied to the coordinates of each atom.  Sz must be 0.0
for a 2d simulation.  This is a mechanism for adding structured
collections of atoms at different locations within the simulation box,
to build up a complex geometry.  It is up to you to insure atoms do
not end up overlapping unphysically which would lead to bad dynamics.
Note that the :doc:`displace_atoms <displace_atoms>` command can be used
to move a subset of atoms after they have been read from a data file.
Likewise, the :doc:`delete_atoms <delete_atoms>` command can be used to
remove overlapping atoms.  Note that the shift values (Sx, Sy, Sz) are
also added to the simulation box information (xlo, xhi, ylo, yhi, zlo,
zhi) in the data file to shift its boundaries.  E.g. xlo_new = xlo +
Sx, xhi_new = xhi + Sx.

The *extra* keywords can only be used the first time the read_data
command is used.  They are useful if you intend to add new atom, bond,
angle, etc types later with additional read_data commands.  This is
because the maximum number of allowed atom, bond, angle, etc types is
set by LAMMPS when the system is first initialized.  If you do not use
the *extra* keywords, then the number of these types will be limited
to what appears in the first data file you read.  For example, if the
first data file is a solid substrate of Si, it will likely specify a
single atom type.  If you read a second data file with a different
material (water molecules) that sit on top of the substrate, you will
want to use different atom types for those atoms.  You can only do
this if you set the *extra/atom/types* keyword to a sufficiently large
value when reading the substrate data file.  Note that use of the
*extra* keywords also allows each data file to contain sections like
Masses or Pair Coeffs or Bond Coeffs which are sized appropriately for
the number of types in that data file.  If the *offset* keyword is
used appropriately when each data file is read, the values in those
sections will be stored correctly in the larger data structures
allocated by the use of the *extra* keywords.  E.g. the substrate file
can list mass and pair coefficients for type 1 silicon atoms.  The
water file can list mass and pair coefficients for type 1 and type 2
hydrogen and oxygen atoms.  Use of the *extra* and *offset* keywords
will store those mass and pair coefficient values appropriately in
data structures that allow for 3 atom types (Si, H, O).  Of course,
you would still need to specify coefficients for H/Si and O/Si
interactions in your input script to have a complete pairwise
interaction model.

An alternative to using the *extra* keywords with the read_data
command, is to use the :doc:`create_box <create_box>` command to
initialize the simulation box and all the various type limits you need
via its *extra* keywords.  Then use the read_data command one or more
times to populate the system with atoms, bonds, angles, etc, using the
*offset* keyword if desired to alter types used in the various data
files you read.

----------

Format of a data file
"""""""""""""""""""""

The structure of the data file is important, though many settings and
sections are optional or can come in any order.  See the examples
directory for sample data files for different problems.

The file will be read line by line, but there is a limit of 254
characters per line and characters beyond that limit will be ignored.

A data file has a header and a body.  The header appears first.  The
first line of the header is always skipped; it typically contains a
description of the file.  Then lines are read one at a time.  Lines
can have a trailing comment starting with '#' that is ignored.  If the
line is blank (only white-space after comment is deleted), it is
skipped.  If the line contains a header keyword, the corresponding
value(s) is read from the line.  If it does not contain a header
keyword, the line begins the body of the file.

The body of the file contains zero or more sections.  The first line
of a section has only a keyword.  This line can have a trailing
comment starting with '#' that is either ignored or can be used to
check for a style match, as described below.  The next line is
skipped.  The remaining lines of the section contain values.  The
number of lines depends on the section keyword as described below.
Zero or more blank lines can be used between sections.  Sections can
appear in any order, with a few exceptions as noted below.

The keyword *fix* can be used one or more times.  Each usage specifies
a fix that will be used to process a specific portion of the data
file.  Any header line containing *header-string* and any section that
is an exact match to *section-string* will be passed to the specified
fix.  See the :doc:`fix property/atom <fix_property_atom>` command for
an example of a fix that operates in this manner.  The doc page for
the fix defines the syntax of the header line(s) and section that it
reads from the data file.  Note that the *header-string* can be
specified as NULL, in which case no header lines are passed to the
fix.  This means the fix can infer the length of its Section from
standard header settings, such as the number of atoms.

The formatting of individual lines in the data file (indentation,
spacing between words and numbers) is not important except that header
and section keywords (e.g. atoms, xlo xhi, Masses, Bond Coeffs) must
be capitalized as shown and cannot have extra white-space between
their words - e.g. two spaces or a tab between the 2 words in "xlo
xhi" or the 2 words in "Bond Coeffs", is not valid.

----------

Format of the header of a data file
"""""""""""""""""""""""""""""""""""

These are the recognized header keywords.  Header lines can come in
any order.  The value(s) are read from the beginning of the line.
Thus the keyword *atoms* should be in a line like "1000 atoms"; the
keyword *ylo yhi* should be in a line like "-10.0 10.0 ylo yhi"; the
keyword *xy xz yz* should be in a line like "0.0 5.0 6.0 xy xz yz".
All these settings have a default value of 0, except the lo/hi box
size defaults are -0.5 and 0.5.  A line need only appear if the value
is different than the default.

* *atoms* = # of atoms in system
* *bonds* = # of bonds in system
* *angles* = # of angles in system
* *dihedrals* = # of dihedrals in system
* *impropers* = # of impropers in system
* *atom types* = # of atom types in system
* *bond types* = # of bond types in system
* *angle types* = # of angle types in system
* *dihedral types* = # of dihedral types in system
* *improper types* = # of improper types in system
* *extra bond per atom* = leave space for this many new bonds per atom (deprecated, use extra/bond/per/atom keyword)
* *extra angle per atom* = leave space for this many new angles per atom (deprecated, use extra/angle/per/atom keyword)
* *extra dihedral per atom* = leave space for this many new dihedrals per atom (deprecated, use extra/dihedral/per/atom keyword)
* *extra improper per atom* = leave space for this many new impropers per atom (deprecated, use extra/improper/per/atom keyword)
* *extra special per atom* = leave space for this many new special bonds per atom (deprecated, use extra/special/per/atom keyword)
* *ellipsoids* = # of ellipsoids in system
* *lines* = # of line segments in system
* *triangles* = # of triangles in system
* *bodies* = # of bodies in system
* *xlo xhi* = simulation box boundaries in x dimension
* *ylo yhi* = simulation box boundaries in y dimension
* *zlo zhi* = simulation box boundaries in z dimension
* *xy xz yz* = simulation box tilt factors for triclinic system

The initial simulation box size is determined by the lo/hi settings.
In any dimension, the system may be periodic or non-periodic; see the
:doc:`boundary <boundary>` command.  When the simulation box is created
it is also partitioned into a regular 3d grid of rectangular bricks,
one per processor, based on the number of processors being used and
the settings of the :doc:`processors <processors>` command.  The
partitioning can later be changed by the :doc:`balance <balance>` or
:doc:`fix balance <fix_balance>` commands.

If the *xy xz yz* line does not appear, LAMMPS will set up an
axis-aligned (orthogonal) simulation box.  If the line does appear,
LAMMPS creates a non-orthogonal simulation domain shaped as a
parallelepiped with triclinic symmetry.  The parallelepiped has its
"origin" at (xlo,ylo,zlo) and is defined by 3 edge vectors starting
from the origin given by A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C =
(xz,yz,zhi-zlo).  *Xy,xz,yz* can be 0.0 or positive or negative values
and are called "tilt factors" because they are the amount of
displacement applied to faces of an originally orthogonal box to
transform it into the parallelepiped.

By default, the tilt factors (xy,xz,yz) can not skew the box more than
half the distance of the corresponding parallel box length.  For
example, if xlo = 2 and xhi = 12, then the x box length is 10 and the
xy tilt factor must be between -5 and 5.  Similarly, both xz and yz
must be between -(xhi-xlo)/2 and +(yhi-ylo)/2.  Note that this is not
a limitation, since if the maximum tilt factor is 5 (as in this
example), then configurations with tilt = ..., -15, -5, 5, 15, 25,
... are all geometrically equivalent.  If you wish to define a box
with tilt factors that exceed these limits, you can use the :doc:`box tilt <box>` command, with a setting of *large*\ ; a setting of
*small* is the default.

See the :doc:`Howto triclinic <Howto_triclinic>` page for a
geometric description of triclinic boxes, as defined by LAMMPS, and
how to transform these parameters to and from other commonly used
triclinic representations.

When a triclinic system is used, the simulation domain should normally
be periodic in the dimension that the tilt is applied to, which is
given by the second dimension of the tilt factor (e.g. y for xy tilt).
This is so that pairs of atoms interacting across that boundary will
have one of them shifted by the tilt factor.  Periodicity is set by
the :doc:`boundary <boundary>` command.  For example, if the xy tilt
factor is non-zero, then the y dimension should be periodic.
Similarly, the z dimension should be periodic if xz or yz is non-zero.
LAMMPS does not require this periodicity, but you may lose atoms if
this is not the case.

Also note that if your simulation will tilt the box, e.g. via the
:doc:`fix deform <fix_deform>` command, the simulation box must be setup
to be triclinic, even if the tilt factors are initially 0.0.  You can
also change an orthogonal box to a triclinic box or vice versa by using
the :doc:`change box <change_box>` command with its *ortho* and
*triclinic* options.

For 2d simulations, the *zlo zhi* values should be set to bound the z
coords for atoms that appear in the file; the default of -0.5 0.5 is
valid if all z coords are 0.0.  For 2d triclinic simulations, the xz
and yz tilt factors must be 0.0.

If the system is periodic (in a dimension), then atom coordinates can
be outside the bounds (in that dimension); they will be remapped (in a
periodic sense) back inside the box.  Note that if the *add* option is
being used to add atoms to a simulation box that already exists, this
periodic remapping will be performed using simulation box bounds that
are the union of the existing box and the box boundaries in the new
data file.

If the system is non-periodic (in a dimension), then an image flag for
that direction has no meaning, since there cannot be periodic images
without periodicity and the data file is therefore - technically speaking
- invalid.  This situation would happen when a data file was written
with periodic boundaries and then read back for non-periodic boundaries.
Accepting a non-zero image flag can lead to unexpected results for any
operations and computations in LAMMPS that internally use unwrapped
coordinates (for example computing the center of mass of a group of
atoms). Thus all non-zero image flags for non-periodic dimensions will
be be reset to zero on reading the data file and LAMMPS will print a
warning message, if that happens.  This is equivalent to wrapping atoms
individually back into the principal unit cell in that direction.  This
operation is equivalent to the behavior of the :doc:`change_box command
<change_box>` when used to change periodicity.


If those atoms with non-zero image flags are involved in bonded
interactions, this reset can lead to undesired changes, when the image
flag values differ between the atoms, i.e. the bonded interaction
straddles domain boundaries.  For example a bond can become stretched
across the unit cell if one of its atoms is wrapped to one side of the
cell and the second atom to the other. In those cases the data file
needs to be pre-processed externally to become valid again.  This can be
done by first unwrapping coordinates and then wrapping entire molecules
instead of individual atoms back into the principal simulation cell and
finally expanding the cell dimensions in the non-periodic direction as
needed, so that the image flag would be zero.

.. note::

   If the system is non-periodic (in a dimension), then all atoms in the
   data file must have coordinates (in that dimension) that are "greater
   than or equal to" the lo value and "less than or equal to" the hi
   value.  If the non-periodic dimension is of style "fixed" (see the
   :doc:`boundary <boundary>` command), then the atom coords must be
   strictly "less than" the hi value, due to the way LAMMPS assign atoms
   to processors.  Note that you should not make the lo/hi values
   radically smaller/larger than the extent of the atoms.  For example,
   if your atoms extend from 0 to 50, you should not specify the box
   bounds as -10000 and 10000 unless you also use the :doc:`processors
   command <processors>`.  This is because LAMMPS uses the specified box
   size to layout the 3d grid of processors.  A huge (mostly empty) box
   will be sub-optimal for performance when using "fixed" boundary
   conditions (see the :doc:`boundary <boundary>` command).  When using
   "shrink-wrap" boundary conditions (see the :doc:`boundary <boundary>`
   command), a huge (mostly empty) box may cause a parallel simulation
   to lose atoms when LAMMPS shrink-wraps the box around the atoms.  The
   read_data command will generate an error in this case.

The "extra bond per atom" setting (angle, dihedral, improper) is only
needed if new bonds (angles, dihedrals, impropers) will be added to
the system when a simulation runs, e.g. by using the :doc:`fix bond/create <fix_bond_create>` command. Using this header flag
is deprecated; please use the *extra/bond/per/atom* keyword (and
correspondingly for angles, dihedrals and impropers) in the read_data
command instead. Either will pre-allocate space in LAMMPS data
structures for storing the new bonds (angles, dihedrals, impropers).

The "extra special per atom" setting is typically only needed if new
bonds/angles/etc will be added to the system, e.g. by using the :doc:`fix bond/create <fix_bond_create>` command.  Or if entire new molecules
will be added to the system, e.g. by using the
:doc:`fix deposit <fix_deposit>` or :doc:`fix pour <fix_pour>` commands,
which will have more special 1-2,1-3,1-4 neighbors than any other
molecules defined in the data file.  Using this header flag is
deprecated; please use the *extra/special/per/atom* keyword instead.
Using this setting will pre-allocate space in the LAMMPS data
structures for storing these neighbors.  See the
:doc:`special_bonds <special_bonds>` and :doc:`molecule <molecule>` doc
pages for more discussion of 1-2,1-3,1-4 neighbors.

.. note::

   All of the "extra" settings are only applied in the first data
   file read and when no simulation box has yet been created; as soon as
   the simulation box is created (and read_data implies that), these
   settings are *locked* and cannot be changed anymore. Please see the
   description of the *add* keyword above for reading multiple data files.
   If they appear in later data files, they are ignored.

The "ellipsoids" and "lines" and "triangles" and "bodies" settings are
only used with :doc:`atom_style ellipsoid or line or tri or body <atom_style>` and specify how many of the atoms are
finite-size ellipsoids or lines or triangles or bodies; the remainder
are point particles.  See the discussion of ellipsoidflag and the
*Ellipsoids* section below.  See the discussion of lineflag and the
*Lines* section below.  See the discussion of triangleflag and the
*Triangles* section below.  See the discussion of bodyflag and the
*Bodies* section below.

.. note::

   For :doc:`atom_style template <atom_style>`, the molecular
   topology (bonds,angles,etc) is contained in the molecule templates
   read-in by the :doc:`molecule <molecule>` command.  This means you
   cannot set the *bonds*, *angles*, etc header keywords in the data
   file, nor can you define *Bonds*, *Angles*, etc sections as discussed
   below.  You can set the *bond types*, *angle types*, etc header
   keywords, though it is not necessary.  If specified, they must match
   the maximum values defined in any of the template molecules.

----------

Format of the body of a data file
"""""""""""""""""""""""""""""""""

These are the section keywords for the body of the file.

* *Atoms, Velocities, Masses, Ellipsoids, Lines, Triangles, Bodies* = atom-property sections
* *Bonds, Angles, Dihedrals, Impropers* = molecular topology sections
* *Pair Coeffs, PairIJ Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs,    Improper Coeffs* = force field sections
* *BondBond Coeffs, BondAngle Coeffs, MiddleBondTorsion Coeffs,    EndBondTorsion Coeffs, AngleTorsion Coeffs, AngleAngleTorsion Coeffs,    BondBond13 Coeffs, AngleAngle Coeffs* = class 2 force field sections

These keywords will check an appended comment for a match with the
currently defined style:

* *Atoms, Pair Coeffs, PairIJ Coeffs, Bond Coeffs, Angle Coeffs, Dihedral Coeffs, Improper Coeffs*

For example, these lines:

.. parsed-literal::

   Atoms # sphere
   Pair Coeffs # lj/cut

will check if the currently-defined :doc:`atom_style <atom_style>` is
*sphere*, and the current :doc:`pair_style <pair_style>` is *lj/cut*\ .
If not, LAMMPS will issue a warning to indicate that the data file
section likely does not contain the correct number or type of
parameters expected for the currently-defined style.

Each section is listed below in alphabetic order.  The format of each
section is described including the number of lines it must contain and
rules (if any) for where it can appear in the data file.

Any individual line in the various sections can have a trailing
comment starting with "#" for annotation purposes.  E.g. in the
Atoms section:

.. parsed-literal::

   10 1 17 -1.0 10.0 5.0 6.0   # salt ion

----------

*Angle Coeffs* section:

* one line per angle type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = angle type (1-N)
       coeffs = list of coeffs

* example:

  .. parsed-literal::

       6 70 108.5 0 0

The number and meaning of the coefficients are specific to the defined
angle style.  See the :doc:`angle_style <angle_style>` and
:doc:`angle_coeff <angle_coeff>` commands for details.  Coefficients can
also be set via the :doc:`angle_coeff <angle_coeff>` command in the
input script.

----------

*AngleAngle Coeffs* section:

* one line per improper type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = improper type (1-N)
       coeffs = list of coeffs (see :doc:`improper_coeff <improper_coeff>`)

----------

*AngleAngleTorsion Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs (see :doc:`dihedral_coeff <dihedral_coeff>`)

----------

*Angles* section:

* one line per angle
* line syntax: ID type atom1 atom2 atom3

  .. parsed-literal::

       ID = number of angle (1-Nangles)
       type = angle type (1-Nangletype)
       atom1,atom2,atom3 = IDs of 1st,2nd,3rd atom in angle

example:

  .. parsed-literal::

       2 2 17 29 430

The 3 atoms are ordered linearly within the angle.  Thus the central
atom (around which the angle is computed) is the atom2 in the list.
E.g. H,O,H for a water molecule.  The *Angles* section must appear
after the *Atoms* section.  All values in this section must be
integers (1, not 1.0).

----------

*AngleTorsion Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs (see :doc:`dihedral_coeff <dihedral_coeff>`)

----------

*Atoms* section:

* one line per atom
* line syntax: depends on atom style

An *Atoms* section must appear in the data file if natoms > 0 in the
header section.  The atoms can be listed in any order.  These are the
line formats for each :doc:`atom style <atom_style>` in LAMMPS.  As
discussed below, each line can optionally have 3 flags (nx,ny,nz)
appended to it, which indicate which image of a periodic simulation
box the atom is in.  These may be important to include for some kinds
of analysis.

.. list-table::

   * - angle
     - atom-ID molecule-ID atom-type x y z
   * - atomic
     - atom-ID atom-type x y z
   * - body
     - atom-ID atom-type bodyflag mass x y z
   * - bond
     - atom-ID molecule-ID atom-type x y z
   * - charge
     - atom-ID atom-type q x y z
   * - dipole
     - atom-ID atom-type q x y z mux muy muz
   * - dpd
     - atom-ID atom-type theta x y z
   * - edpd
     - atom-ID atom-type edpd_temp edpd_cv x y z
   * - electron
     - atom-ID atom-type q spin eradius x y z
   * - ellipsoid
     - atom-ID atom-type ellipsoidflag density x y z
   * - full
     - atom-ID molecule-ID atom-type q x y z
   * - line
     - atom-ID molecule-ID atom-type lineflag density x y z
   * - mdpd
     - atom-ID atom-type rho x y z
   * - mesont
     - atom-ID molecule-ID atom-type bond_nt mass mradius mlength buckling x y z
   * - molecular
     - atom-ID molecule-ID atom-type x y z
   * - peri
     - atom-ID atom-type volume density x y z
   * - smd
     - atom-ID atom-type molecule volume mass kradius cradius x0 y0 z0 x y z
   * - sph
     - atom-ID atom-type rho esph cv x y z
   * - sphere
     - atom-ID atom-type diameter density x y z
   * - spin
     - atom-ID atom-type x y z spx spy spz sp
   * - tdpd
     - atom-ID atom-type x y z cc1 cc2 ... ccNspecies
   * - template
     - atom-ID atom-type molecule-ID template-index template-atom x y z
   * - tri
     - atom-ID molecule-ID atom-type triangleflag density x y z
   * - wavepacket
     - atom-ID atom-type charge spin eradius etag cs_re cs_im x y z
   * - hybrid
     - atom-ID atom-type x y z sub-style1 sub-style2 ...

The per-atom values have these meanings and units, listed alphabetically:

* atom-ID = integer ID of atom
* atom-type = type of atom (1-Ntype)
* bodyflag = 1 for body particles, 0 for point particles
* bond_nt = bond NT factor for MESONT particles (?? units)
* buckling = buckling factor for MESONT particles (?? units)
* ccN = chemical concentration for tDPD particles for each species (mole/volume units)
* cradius = contact radius for SMD particles (distance units)
* cs_re,cs_im = real/imaginary parts of wave packet coefficients
* cv = heat capacity (need units) for SPH particles
* density = density of particle (mass/distance\^3 or mass/distance\^2 or mass/distance units, depending on dimensionality of particle)
* diameter = diameter of spherical atom (distance units)
* esph = energy (need units) for SPH particles
* edpd_temp = temperature for eDPD particles (temperature units)
* edpd_cv = volumetric heat capacity for eDPD particles (energy/temperature/volume units)
* ellipsoidflag = 1 for ellipsoidal particles, 0 for point particles
* eradius = electron radius (or fixed-core radius)
* etag = integer ID of electron that each wave packet belongs to
* kradius = kernel radius for SMD particles (distance units)
* lineflag = 1 for line segment particles, 0 for point or spherical particles
* mass = mass of particle (mass units)
* mlength = ?? length for MESONT particles (distance units)
* molecule-ID = integer ID of molecule the atom belongs to
* mradius = ?? radius for MESONT particles (distance units)
* mux,muy,muz = components of dipole moment of atom (dipole units)
* q = charge on atom (charge units)
* rho = density (need units) for SPH particles
* spin = electron spin (+1/-1), 0 = nuclei, 2 = fixed-core, 3 = pseudo-cores (i.e. ECP)
* sp = magnitude of magnetic spin of atom (Bohr magnetons)
* spx,spy,spz = components of magnetic spin of atom (unit vector)
* template-atom = which atom within a template molecule the atom is
* template-index = which molecule within the molecule template the atom is part of
* theta = internal temperature of a DPD particle
* triangleflag = 1 for triangular particles, 0 for point or spherical particles
* volume = volume of Peridynamic particle (distance\^3 units)
* x,y,z = coordinates of atom (distance units)
* x0,y0,z0 = original (strain-free) coordinates of atom (distance units)

The units for these quantities depend on the unit style; see the
:doc:`units <units>` command for details.

For 2d simulations specify z as 0.0, or a value within the *zlo zhi*
setting in the data file header.

The atom-ID is used to identify the atom throughout the simulation and
in dump files.  Normally, it is a unique value from 1 to Natoms for
each atom.  Unique values larger than Natoms can be used, but they
will cause extra memory to be allocated on each processor, if an atom
map array is used, but not if an atom map hash is used; see the
:doc:`atom_modify <atom_modify>` command for details.  If an atom map is
not used (e.g. an atomic system with no bonds), and you don't care if
unique atom IDs appear in dump files, then the atom-IDs can all be set
to 0.

The molecule ID is a second identifier attached to an atom.  Normally, it
is a number from 1 to N, identifying which molecule the atom belongs
to.  It can be 0 if it is a non-bonded atom or if you don't care to
keep track of molecule assignments.

The diameter specifies the size of a finite-size spherical particle.
It can be set to 0.0, which means that atom is a point particle.

The ellipsoidflag, lineflag, triangleflag, and bodyflag determine
whether the particle is a finite-size ellipsoid or line or triangle or
body of finite size, or whether the particle is a point particle.
Additional attributes must be defined for each ellipsoid, line,
triangle, or body in the corresponding *Ellipsoids*, *Lines*,
*Triangles*, or *Bodies* section.

The *template-index* and *template-atom* are only defined used by
:doc:`atom_style template <atom_style>`.  In this case the
:doc:`molecule <molecule>` command is used to define a molecule template
which contains one or more molecules (as separate files).  If an atom
belongs to one of those molecules, its *template-index* and *template-atom*
are both set to positive integers; if not the values are both 0.  The
*template-index* is which molecule (1 to Nmols) the atom belongs to.
The *template-atom* is which atom (1 to Natoms) within the molecule
the atom is.

Some pair styles and fixes and computes that operate on finite-size
particles allow for a mixture of finite-size and point particles.  See
the doc pages of individual commands for details.

For finite-size particles, the density is used in conjunction with the
particle volume to set the mass of each particle as mass = density \*
volume.  In this context, volume can be a 3d quantity (for spheres or
ellipsoids), a 2d quantity (for triangles), or a 1d quantity (for line
segments).  If the volume is 0.0, meaning a point particle, then the
density value is used as the mass.  One exception is for the body atom
style, in which case the mass of each particle (body or point
particle) is specified explicitly.  This is because the volume of the
body is unknown.

Note that for 2d simulations of spheres, this command will treat them
as spheres when converting density to mass.  However, they can also be
modeled as 2d discs (circles) if the :doc:`set density/disc <set>`
command is used to reset their mass after the read_data command is
used.  A *disc* keyword can also be used with time integration fixes,
such as :doc:`fix nve/sphere <fix_nve_sphere>` and :doc:`fix
nvt/sphere <fix_nve_sphere>` to time integrate their motion as 2d
discs (not 3d spheres), by changing their moment of inertia.

For atom\_style hybrid, following the 5 initial values
(ID,type,x,y,z), specific values for each sub-style must be listed.
The order of the sub-styles is the same as they were listed in the
:doc:`atom_style <atom_style>` command.  The specific values for each
sub-style are those that are not the 5 standard ones (ID,type,x,y,z).
For example, for the "charge" sub-style, a "q" value would appear.
For the "full" sub-style, a "molecule-ID" and "q" would appear.  These
are listed in the same order they appear as listed above.  Thus if

.. parsed-literal::

   atom_style hybrid charge sphere

were used in the input script, each atom line would have these fields:

.. parsed-literal::

   atom-ID atom-type x y z q diameter density

Note that if a non-standard value is defined by multiple sub-styles,
it only appears once in the atom line.  E.g. the atom line for
atom_style hybrid dipole full would list "q" only once, with the
dipole sub-style fields; "q" does not appear with the full sub-style
fields.

.. parsed-literal::

   atom-ID atom-type x y z q mux muy myz molecule-ID

Atom lines specify the (x,y,z) coordinates of atoms.  These can be
inside or outside the simulation box.  When the data file is read,
LAMMPS wraps coordinates outside the box back into the box for
dimensions that are periodic.  As discussed above, if an atom is
outside the box in a non-periodic dimension, it will be lost.

LAMMPS always stores atom coordinates as values which are inside the
simulation box.  It also stores 3 flags which indicate which image of
the simulation box (in each dimension) the atom would be in if its
coordinates were unwrapped across periodic boundaries.  An image flag
of 0 means the atom is still inside the box when unwrapped.  A value
of 2 means add 2 box lengths to get the unwrapped coordinate.  A value
of -1 means subtract 1 box length to get the unwrapped coordinate.
LAMMPS updates these flags as atoms cross periodic boundaries during
the simulation.  The :doc:`dump <dump>` command can output atom
coordinates in wrapped or unwrapped form, as well as the 3 image
flags.

In the data file, atom lines (all lines or none of them) can
optionally list 3 trailing integer values (nx,ny,nz), which are used
to initialize the atom's image flags.  If nx,ny,nz values are not
listed in the data file, LAMMPS initializes them to 0.  Note that the
image flags are immediately updated if an atom's coordinates need to
wrapped back into the simulation box.

It is only important to set image flags correctly in a data file if a
simulation model relies on unwrapped coordinates for some calculation;
otherwise they can be left unspecified.  Examples of LAMMPS commands
that use unwrapped coordinates internally are as follows:

* Atoms in a rigid body (see :doc:`fix rigid <fix_rigid>`, :doc:`fix rigid/small <fix_rigid>`) must have consistent image flags, so that
  when the atoms are unwrapped, they are near each other, i.e. as a
  single body.
* If the :doc:`replicate <replicate>` command is used to generate a larger
  system, image flags must be consistent for bonded atoms when the bond
  crosses a periodic boundary.  I.e. the values of the image flags
  should be different by 1 (in the appropriate dimension) for the two
  atoms in such a bond.
* If you plan to :doc:`dump <dump>` image flags and perform post-analysis
  that will unwrap atom coordinates, it may be important that a
  continued run (restarted from a data file) begins with image flags
  that are consistent with the previous run.

.. note::

   If your system is an infinite periodic crystal with bonds then
   it is impossible to have fully consistent image flags.  This is because
   some bonds will cross periodic boundaries and connect two atoms with the
   same image flag.

Atom velocities and other atom quantities not defined above are set to
0.0 when the *Atoms* section is read.  Velocities can be set later by
a *Velocities* section in the data file or by a
:doc:`velocity <velocity>` or :doc:`set <set>` command in the input
script.

----------

*Bodies* section:

* one or more lines per body
* first line syntax: atom-ID Ninteger Ndouble

  .. parsed-literal::

       Ninteger = # of integer quantities for this particle
       Ndouble = # of floating-point quantities for this particle

* 0 or more integer lines with total of Ninteger values
* 0 or more double lines with total of Ndouble values
* example:

  .. parsed-literal::

       12 3 6
       2 3 2
       1.0 2.0 3.0 1.0 2.0 4.0

* example:

  .. parsed-literal::

       12 0 14
       1.0 2.0 3.0 1.0 2.0 4.0 1.0
       2.0 3.0 1.0 2.0 4.0 4.0 2.0

The *Bodies* section must appear if :doc:`atom_style body <atom_style>`
is used and any atoms listed in the *Atoms* section have a bodyflag =
1.  The number of bodies should be specified in the header section via
the "bodies" keyword.

Each body can have a variable number of integer and/or floating-point
values.  The number and meaning of the values is defined by the body
style, as described in the :doc:`Howto body <Howto_body>` doc page.  The
body style is given as an argument to the :doc:`atom_style body <atom_style>` command.

The Ninteger and Ndouble values determine how many integer and
floating-point values are specified for this particle.  Ninteger and
Ndouble can be as large as needed and can be different for every body.
Integer values are then listed next on subsequent lines.  Lines are
read one at a time until Ninteger values are read.  Floating-point
values follow on subsequent lines, Again lines are read one at a time
until Ndouble values are read.  Note that if there are no values of a
particular type, no lines appear for that type.

The *Bodies* section must appear after the *Atoms* section.

----------

*Bond Coeffs* section:

* one line per bond type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = bond type (1-N)
       coeffs = list of coeffs

* example:

  .. parsed-literal::

       4 250 1.49

The number and meaning of the coefficients are specific to the defined
bond style.  See the :doc:`bond_style <bond_style>` and
:doc:`bond_coeff <bond_coeff>` commands for details.  Coefficients can
also be set via the :doc:`bond_coeff <bond_coeff>` command in the input
script.

----------

*BondAngle Coeffs* section:

* one line per angle type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = angle type (1-N)
       coeffs = list of coeffs (see class 2 section of :doc:`angle_coeff <angle_coeff>`)

----------

*BondBond Coeffs* section:

* one line per angle type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = angle type (1-N)
       coeffs = list of coeffs (see class 2 section of :doc:`angle_coeff <angle_coeff>`)

----------

*BondBond13 Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs (see class 2 section of :doc:`dihedral_coeff <dihedral_coeff>`)

----------

*Bonds* section:

* one line per bond
* line syntax: ID type atom1 atom2

  .. parsed-literal::

       ID = bond number (1-Nbonds)
       type = bond type (1-Nbondtype)
       atom1,atom2 = IDs of 1st,2nd atom in bond

* example:

  .. parsed-literal::

       12 3 17 29

The *Bonds* section must appear after the *Atoms* section.  All values
in this section must be integers (1, not 1.0).

----------

*Dihedral Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs

* example:

  .. parsed-literal::

       3 0.6 1 0 1

The number and meaning of the coefficients are specific to the defined
dihedral style.  See the :doc:`dihedral_style <dihedral_style>` and
:doc:`dihedral_coeff <dihedral_coeff>` commands for details.
Coefficients can also be set via the
:doc:`dihedral_coeff <dihedral_coeff>` command in the input script.

----------

*Dihedrals* section:

* one line per dihedral
* line syntax: ID type atom1 atom2 atom3 atom4

  .. parsed-literal::

       ID = number of dihedral (1-Ndihedrals)
       type = dihedral type (1-Ndihedraltype)
       atom1,atom2,atom3,atom4 = IDs of 1st,2nd,3rd,4th atom in dihedral

* example:

  .. parsed-literal::

       12 4 17 29 30 21

The 4 atoms are ordered linearly within the dihedral.  The *Dihedrals*
section must appear after the *Atoms* section.  All values in this
section must be integers (1, not 1.0).

----------

*Ellipsoids* section:

* one line per ellipsoid
* line syntax: atom-ID shapex shapey shapez quatw quati quatj quatk

  .. parsed-literal::

       atom-ID = ID of atom which is an ellipsoid
       shapex,shapey,shapez = 3 diameters of ellipsoid (distance units)
       quatw,quati,quatj,quatk = quaternion components for orientation of atom

* example:

  .. parsed-literal::

       12 1 2 1 1 0 0 0

The *Ellipsoids* section must appear if :doc:`atom_style ellipsoid <atom_style>` is used and any atoms are listed in the
*Atoms* section with an ellipsoidflag = 1.  The number of ellipsoids
should be specified in the header section via the "ellipsoids"
keyword.

The 3 shape values specify the 3 diameters or aspect ratios of a
finite-size ellipsoidal particle, when it is oriented along the 3
coordinate axes.  They must all be non-zero values.

The values *quatw*, *quati*, *quatj*, and *quatk* set the orientation
of the atom as a quaternion (4-vector).  Note that the shape
attributes specify the aspect ratios of an ellipsoidal particle, which
is oriented by default with its x-axis along the simulation box's
x-axis, and similarly for y and z.  If this body is rotated (via the
right-hand rule) by an angle theta around a unit vector (a,b,c), then
the quaternion that represents its new orientation is given by
(cos(theta/2), a\*sin(theta/2), b\*sin(theta/2), c\*sin(theta/2)).  These
4 components are quatw, quati, quatj, and quatk as specified above.
LAMMPS normalizes each atom's quaternion in case (a,b,c) is not
specified as a unit vector.

The *Ellipsoids* section must appear after the *Atoms* section.

----------

*EndBondTorsion Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs (see class 2 section of :doc:`dihedral_coeff <dihedral_coeff>`)

----------

*Improper Coeffs* section:

* one line per improper type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = improper type (1-N)
       coeffs = list of coeffs

* example:

  .. parsed-literal::

       2 20 0.0548311

The number and meaning of the coefficients are specific to the defined
improper style.  See the :doc:`improper_style <improper_style>` and
:doc:`improper_coeff <improper_coeff>` commands for details.
Coefficients can also be set via the
:doc:`improper_coeff <improper_coeff>` command in the input script.

----------

*Impropers* section:

* one line per improper
* line syntax: ID type atom1 atom2 atom3 atom4

  .. parsed-literal::

       ID = number of improper (1-Nimpropers)
       type = improper type (1-Nimpropertype)
       atom1,atom2,atom3,atom4 = IDs of 1st,2nd,3rd,4th atom in improper

* example:

  .. parsed-literal::

       12 3 17 29 13 100

The ordering of the 4 atoms determines the definition of the improper
angle used in the formula for each :doc:`improper style <improper_style>`.  See the doc pages for individual styles
for details.

The *Impropers* section must appear after the *Atoms* section.  All
values in this section must be integers (1, not 1.0).

----------

*Lines* section:

* one line per line segment
* line syntax: atom-ID x1 y1 x2 y2

  .. parsed-literal::

       atom-ID = ID of atom which is a line segment
       x1,y1 = 1st end point
       x2,y2 = 2nd end point

* example:

  .. parsed-literal::

       12 1.0 0.0 2.0 0.0

The *Lines* section must appear if :doc:`atom_style line <atom_style>`
is used and any atoms are listed in the *Atoms* section with a
lineflag = 1.  The number of lines should be specified in the header
section via the "lines" keyword.

The 2 end points are the end points of the line segment.  The ordering
of the 2 points should be such that using a right-hand rule to cross
the line segment with a unit vector in the +z direction, gives an
"outward" normal vector perpendicular to the line segment.
I.e. normal = (c2-c1) x (0,0,1).  This orientation may be important
for defining some interactions.

The *Lines* section must appear after the *Atoms* section.

----------

*Masses* section:

* one line per atom type
* line syntax: ID mass

  .. parsed-literal::

       ID = atom type (1-N)
       mass = mass value

* example:

  .. parsed-literal::

       3 1.01

This defines the mass of each atom type.  This can also be set via the
:doc:`mass <mass>` command in the input script.  This section cannot be
used for atom styles that define a mass for individual atoms -
e.g. :doc:`atom_style sphere <atom_style>`.

----------

*MiddleBondTorsion Coeffs* section:

* one line per dihedral type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = dihedral type (1-N)
       coeffs = list of coeffs (see class 2 section of :doc:`dihedral_coeff <dihedral_coeff>`)

----------

*Pair Coeffs* section:

* one line per atom type
* line syntax: ID coeffs

  .. parsed-literal::

       ID = atom type (1-N)
       coeffs = list of coeffs

* example:

  .. parsed-literal::

       3 0.022 2.35197 0.022 2.35197

The number and meaning of the coefficients are specific to the defined
pair style.  See the :doc:`pair_style <pair_style>` and
:doc:`pair_coeff <pair_coeff>` commands for details.  Since pair
coefficients for types I != J are not specified, these will be
generated automatically by the pair style's mixing rule.  See the
individual pair_style doc pages and the :doc:`pair_modify mix
<pair_modify>` command for details.  Pair coefficients can also be set
via the :doc:`pair_coeff <pair_coeff>` command in the input script.

----------

*PairIJ Coeffs* section:

* one line per pair of atom types for all I,J with I <= J
* line syntax: ID1 ID2 coeffs

  .. parsed-literal::

       ID1 = atom type I = 1-N
       ID2 = atom type J = I-N, with I <= J
       coeffs = list of coeffs

* examples:

  .. parsed-literal::

       3 3 0.022 2.35197 0.022 2.35197
       3 5 0.022 2.35197 0.022 2.35197

This section must have N\*(N+1)/2 lines where N = # of atom types.
The number and meaning of the coefficients are specific to the defined
pair style.  See the :doc:`pair_style <pair_style>` and
:doc:`pair_coeff <pair_coeff>` commands for details.  Since pair
coefficients for types I != J are all specified, these values will
turn off the default mixing rule defined by the pair style.  See the
individual pair_style doc pages and the :doc:`pair_modify mix
<pair_modify>` command for details.  Pair coefficients can also be set
via the :doc:`pair_coeff <pair_coeff>` command in the input script.

----------

*Triangles* section:

* one line per triangle
* line syntax: atom-ID x1 y1 z1 x2 y2 z2 x3 y3 z3

  .. parsed-literal::

       atom-ID = ID of atom which is a line segment
       x1,y1,z1 = 1st corner point
       x2,y2,z2 = 2nd corner point
       x3,y3,z3 = 3rd corner point

* example:

  .. parsed-literal::

       12 0.0 0.0 0.0 2.0 0.0 1.0 0.0 2.0 1.0

The *Triangles* section must appear if :doc:`atom_style tri <atom_style>` is used and any atoms are listed in the *Atoms*
section with a triangleflag = 1.  The number of lines should be
specified in the header section via the "triangles" keyword.

The 3 corner points are the corner points of the triangle.  The
ordering of the 3 points should be such that using a right-hand rule
to go from point1 to point2 to point3 gives an "outward" normal vector
to the face of the triangle.  I.e. normal = (c2-c1) x (c3-c1).  This
orientation may be important for defining some interactions.

The *Triangles* section must appear after the *Atoms* section.

----------

*Velocities* section:

* one line per atom
* line syntax: depends on atom style

+--------------------------------+--------------------------------------------+
| all styles except those listed | atom-ID vx vy vz                           |
+--------------------------------+--------------------------------------------+
| electron                       | atom-ID vx vy vz ervel                     |
+--------------------------------+--------------------------------------------+
| ellipsoid                      | atom-ID vx vy vz lx ly lz                  |
+--------------------------------+--------------------------------------------+
| sphere                         | atom-ID vx vy vz wx wy wz                  |
+--------------------------------+--------------------------------------------+
| hybrid                         | atom-ID vx vy vz sub-style1 sub-style2 ... |
+--------------------------------+--------------------------------------------+

where the keywords have these meanings:

vx,vy,vz = translational velocity of atom
lx,ly,lz = angular momentum of aspherical atom
wx,wy,wz = angular velocity of spherical atom
ervel = electron radial velocity (0 for fixed-core):ul

The velocity lines can appear in any order.  This section can only be
used after an *Atoms* section.  This is because the *Atoms* section
must have assigned a unique atom ID to each atom so that velocities
can be assigned to them.

Vx, vy, vz, and ervel are in :doc:`units <units>` of velocity.  Lx, ly,
lz are in units of angular momentum (distance-velocity-mass).  Wx, Wy,
Wz are in units of angular velocity (radians/time).

For atom_style hybrid, following the 4 initial values (ID,vx,vy,vz),
specific values for each sub-style must be listed.  The order of the
sub-styles is the same as they were listed in the
:doc:`atom_style <atom_style>` command.  The sub-style specific values
are those that are not the 5 standard ones (ID,vx,vy,vz).  For
example, for the "sphere" sub-style, "wx", "wy", "wz" values would
appear.  These are listed in the same order they appear as listed
above.  Thus if

.. code-block:: LAMMPS

   atom_style hybrid electron sphere

were used in the input script, each velocity line would have these
fields:

.. parsed-literal::

   atom-ID vx vy vz ervel wx wy wz

Translational velocities can also be set by the
:doc:`velocity <velocity>` command in the input script.

----------

Restrictions
""""""""""""

To read gzipped data files, you must compile LAMMPS with the
-DLAMMPS_GZIP option.  See the :doc:`Build settings <Build_settings>`
doc page for details.

Related commands
""""""""""""""""

:doc:`read_dump <read_dump>`, :doc:`read_restart <read_restart>`,
:doc:`create_atoms <create_atoms>`, :doc:`write_data <write_data>`

Default
"""""""

The default for all the *extra* keywords is 0.
