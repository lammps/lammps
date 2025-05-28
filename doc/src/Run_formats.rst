
File formats used by LAMMPS
===========================

This page provides a general overview of the kinds of files and file
formats that LAMMPS is reading and writing.

.. contents:: On this page
   :depth: 2
   :backlinks: top

-------------------

Character Encoding
^^^^^^^^^^^^^^^^^^

For processing text files, the LAMMPS source code assumes `ASCII
character encoding <https://en.wikipedia.org/wiki/ASCII>`_ which
represents the digits 0 to 9, the lower and upper case letters a to z,
some common punctuation and other symbols and a few whitespace
characters including a regular "space character", "line feed", "carriage
return", "tabulator".  These characters are all represented by single
bytes with a value smaller than 128 and only 95 of those 128 values
represent printable characters.  This list is sufficient to represent
most English text, but misses accented characters or umlauts or Greek
symbols and more.

Modern text often uses `UTF-8 character encoding
<https://en.wikipedia.org/wiki/UTF-8>`_ instead.  This encoding is a way
to represent many more different characters as defined by the Unicode
standard.  UFT-8 is compatible with ASCII, since the first 128 values
are identical with the ASCII encoding.  It is important to note,
however, that there are Unicode characters that *look* similar to ASCII
characters, but have a different binary representation.  As a general
rule, these characters may not be correctly recognized by LAMMPS.  For
some parts of LAMMPS' text processing, translation tables with known
"lookalike" characters are used.  The tables are used to substitute
non-ASCII characters with their ASCII equivalents.  Non-ASCII lookalike
characters are often used by web browsers or PDF viewers to improve the
readability of text.  Thus, when using copy and paste to transfer text
from such an application to your input file, you may unintentionally
create text that is not exclusively using ASCII encoding and may cause
errors when LAMMPS is trying to read it.

Lines with non-printable and non-ASCII characters in text files can be
detected for example with a (Linux) command like the following:

.. code-block:: bash

   env LC_ALL=C grep -n '[^ -~]' some_file.txt

Number Formatting
^^^^^^^^^^^^^^^^^

Different countries and languages have different conventions to format
numbers.  While in some regions commas are used for fractions and points
to indicate thousand, million and so on, this is reversed in other
regions.  Modern operating systems have facilities to adjust input and
output accordingly that are collectively referred to as "native language
support" (NLS).  The exact rules are often applied according to the
value of the ``$LANG`` environment variable (e.g. "en_US.utf8" for
English text in UTF-8 encoding).

For the sake of simplicity of the implementation and transferability of
results, LAMMPS does not support this and instead expects numbers being
formatted in the generic or "C" locale.  The "C" locale has no
punctuation for thousand, million and so on and uses a decimal point for
fractions.  One thousand would be represented as "1000.0" and not as
"1,000.0" nor as "1.000,0".  Having native language support enabled for
a locale other than "C" will result in different behavior when
converting or formatting numbers that can trigger unexpected errors.

LAMMPS also only accepts integer numbers when an integer is required, so
using floating point equivalents like "1.0" are not accepted; you *must*
use "1" instead.

For floating point numbers in scientific notation, the Fortran double
precision notation "1.1d3" is not accepted; you have to use "1100",
"1100.0" or "1.1e3".

Input file
^^^^^^^^^^

A LAMMPS input file is a text file with commands.  It is read
line-by-line and each line is processed *immediately*.  Before looking
for commands and executing them, there is a pre-processing step where
comments (non-quoted text starting with a pound sign '#') are removed,
``${variable}`` and ``$(expression)`` constructs are expanded or
evaluated, and lines that end in the ampersand character '&' are
combined with the next line (similar to Fortran 90 free-format source
code).  After the pre-processing, lines are split into "words" and
evaluated.  The first word must be a :doc:`command <Commands_all>` and
all following words are arguments.  Below are some example lines:

.. code-block:: LAMMPS

   # full line comment

   # some global settings
   units           lj
   atom_style      atomic
   # ^^ command    ^^ argument(s)

   variable        x index 1       # may be overridden from command line with -var x <value>
   variable        xx equal 20*$x  # variable "xx" is always 20 times "x"

   lattice         fcc 0.8442

   # example of a command written across multiple lines
   # the "region" command uses spacing from "lattice" command, unless "units box" is specified
   region          box block 0.0 ${xx} &
                             0.0 40.0  &
                             0.0 30.0
   # create simulation box and fill with atoms according to lattice setting
   create_box      1 box
   create_atoms    1 box

   # set force field and parameters
   mass            1 1.0
   pair_style      lj/cut 2.5
   pair_coeff      1 1 1.0 1.0 2.5

   # run simulation
   fix             1 all nve
   run             1000

The pivotal command in this example input is the :doc:`create_box
command <create_box>`.  It defines the simulation system and many
parameters that go with it: units, atom style, number of atom types (and
other types) and more.  Those settings are *locked in* after the box is
created.  Commands that change these kind of settings are only allowed
**before** a simulation box is created and many other commands are only
allowed **after** the simulation box is defined (e.g. :doc:`pair_coeff
<pair_coeff>`).  Very few commands (e.g. :doc:`pair_style <pair_style>`)
may be used in either part of the input.  The :doc:`read_data
<read_data>` and :doc:`read_restart <read_restart>` commands also create
the system box and thus have a similar pivotal function.

The LAMMPS input syntax has minimal support for conditionals and loops,
but if more complex operations are required, it is recommended to use
the library interface, e.g. :doc:`from Python using the LAMMPS Python
module <Python_run>`.

There is a frequent misconception about the :doc:`if command <if>`:
this is a command for conditional execution **outside** a run or
minimization.  To trigger actions on specific conditions **during**
a run is a non-trivial operation that usually requires adopting one
of the available "fix" commands or creating a new "fix" command.

LAMMPS commands change the internal state and thus the order of commands
matters and reordering them can produce different results.  For example,
the region defined by the :doc:`region command <region>` in the example
above depends on the :doc:`lattice setting <lattice>` and thus its
dimensions will be different depending on the order of the two commands.

Each line must have an "end-of-line" character (line feed or carriage
return plus line feed).  Some text editors do not automatically insert
one which may cause LAMMPS to ignore the last command.  It is thus
recommended to always have an empty line at the end of an input file.

The specific details describing how LAMMPS input is processed and parsed
are explained in :doc:`Commands_parse`.

Data file
^^^^^^^^^

A LAMMPS data file contains a description of a system suitable for
reading with the :doc:`read_data command <read_data>`.  Data files are
commonly used for setting up complex molecular systems that can be
difficult to achieve with the commands :doc:`create_box <create_box>`
and :doc:`create_atoms <create_atoms>` alone.  Also, data files can be
used as a portable alternatives to a :doc:`binary restart file
<restart>`.  A restart file can be converted into a data file from the
:doc:`command line <Run_options>`.

Data files have a header section at the very beginning of the file and
multiple titled sections such as "Atoms", Masses", "Pair Coeffs", and so
on.  Header keywords can only be used *before* the first title section.

The data file **always** starts with a "title" line, which will be
**ignored** by LAMMPS.  Omitting the title line can lead to unexpected
behavior because a line of the header with an actual setting may be
ignored.  In this case, the mistakenly ignored line often contains the
"atoms" keyword, which results in LAMMPS assuming that there are no
atoms in the data file and thus throwing an error on the contents of the
"Atoms" section.  The title line may contain some keywords that can be
used by external programs to convey information about the system
(included as comments), that is not required and not read by LAMMPS.

The line following a section title is also **ignored**.  An error will
occur if an empty line is not placed after a section title.  The number
of lines in titled sections depends on header keywords, like the number
of atom types, the number of atoms, the number of bond types, the number
of bonds, and so on.  The data in those sections has to be complete.  A
special case are the "Pair Coeffs" and "PairIJ Coeffs" sections; the
former is for force fields and pair styles that use mixing of non-bonded
potential parameters, the latter for pair styles and force fields
requiring explicit coefficients.  Thus with *N* being the number of atom
types, the "Pair Coeffs" section has *N* entries while "PairIJ Coeffs"
has :math:`N \cdot (N-1)` entries.  Internally, these sections will be
converted to :doc:`pair_coeff <pair_coeff>` commands.  Thus the
corresponding :doc:`pair style <pair_style>` must have been set *before*
the :doc:`read_data command <read_data>` reads the data file.

Data files may contain comments, which start with the pound sign '#'.
There must be at least one blank between a valid keyword and the pound
sign. Below is a simple example case of a data file for :doc:`atom style
full <atom_style>`.

.. code-block:: bash

   LAMMPS Title line (ignored)
   # full line comment

           10  atoms # comment
            4  atom types

    -36.840194 64.211560 xlo xhi
    -41.013691 68.385058 ylo yhi
    -29.768095 57.139462 zlo zhi

   Masses

     1 12.0110
     2 12.0110
     3 15.9990
     4  1.0080

   Pair Coeffs  # this section is optional

     1    0.110000    3.563595    0.110000    3.563595
     2    0.080000    3.670503    0.010000    3.385415
     3    0.120000    3.029056    0.120000    2.494516
     4    0.022000    2.351973    0.022000    2.351973

   Atoms # full

         1      1       1       0.560   43.99993  58.52678  36.78550   0   0   0
         2      1       2      -0.270   45.10395  58.23499  35.86693   0   0   0
         3      1       3      -0.510   43.81519  59.54928  37.43995   0   0   0
         4      1       4       0.090   45.71714  57.34797  36.13434   0   0   0
         5      1       4       0.090   45.72261  59.13657  35.67007   0   0   0
         6      1       4       0.090   44.66624  58.09539  34.85538   0   0   0
         7      1       3      -0.470   43.28193  57.47427  36.91953   0   0   0
         8      1       4       0.070   42.07157  57.45486  37.62418   0   0   0
         9      1       1       0.510   42.19985  57.57789  39.12163   0   0   0
        10      1       1       0.510   41.88641  58.62251  39.70398   0   0   0
   #  ^^atomID ^^molID ^^type  ^^charge ^^xcoord  ^^ycoord  ^^ycoord  ^^image^^flags (optional)

   Velocities # this section is optional

         1  0.0050731  -0.00398928  0.00391473
         2 -0.0175184   0.0173484  -0.00489207
         3  0.00597225 -0.00202006  0.00166454
         4 -0.010395   -0.0082582   0.00316419
         5 -0.00390877  0.00470331 -0.00226911
         6 -0.00111157 -0.00374545 -0.0169374
         7  0.00209054 -0.00594936 -0.000124563
         8  0.00635002 -0.0120093  -0.0110999
         9 -0.004955   -0.0123375   0.000403422
        10  0.00265028 -0.00189329 -0.00293198

The common problem is processing the "Atoms" section, since its format
depends on the :doc:`atom style <atom_style>` used, and that setting
must be done in the input file *before* reading the data file.  To
assist with detecting incompatible data files, a comment is appended to
the "Atoms" title indicating the atom style used (or intended) when
*writing* the data file.  For example, below is an "Atoms" section for
:doc:`atom style charge <atom_style>`, which omits the molecule ID
column.

.. code-block:: bash

   Atoms # charge

         1      1       0.560   43.99993  58.52678  36.78550
         2      2      -0.270   45.10395  58.23499  35.86693
         3      3      -0.510   43.81519  59.54928  37.43995
         4      4       0.090   45.71714  57.34797  36.13434
         5      4       0.090   45.72261  59.13657  35.67007
         6      4       0.090   44.66624  58.09539  34.85538
         7      3      -0.470   43.28193  57.47427  36.91953
         8      4       0.070   42.07157  57.45486  37.62418
         9      1       0.510   42.19985  57.57789  39.12163
        10      1       0.510   41.88641  58.62251  39.70398
   #  ^^atomID ^^type  ^^charge ^^xcoord  ^^ycoord  ^^ycoord

Another source of confusion about the "Atoms" section format is the
ordering of columns.  The three atom style variants `atom_style full`,
`atom_style hybrid charge molecular`, and `atom_style hybrid molecular
charge` all carry the same per-atom information. However, in data files,
the Atoms section has the columns 'Atom-ID Molecule-ID Atom-type Charge
X Y Z' for atom style full, but for hybrid atom styles the first columns
are always 'Atom-ID Atom-type X Y Z' followed by any *additional* data
added by the hybrid styles, for example, 'Charge Molecule-ID' for the
first hybrid style and 'Molecule-ID Charge' in the second hybrid style
variant.  Finally, an alternative to a hybrid atom style is to use fix
property/atom, e.g. to add molecule IDs to atom style charge.  In this
case the "Atoms" section is formatted according to atom style charge and
a new section, "Molecules" is added that contains lines with 'Atom-ID
Molecule-ID', one for each atom in the system.  For adding charges to
atom style molecular with fix property/atom, the "Atoms" section is now
formatted according to the atom style and a "Charges" section is added.

Molecule file
^^^^^^^^^^^^^

Molecule files for use with the :doc:`molecule command <molecule>` look
quite similar to data files but they do not have a compatible format,
i.e., one cannot use a data file as molecule file and vice versa.  Below
is a simple example for a water molecule (SPC/E model).  Same as a data
file, there is an ignored title line and you can use comments.  However,
there is no information about the number of types or the box dimensions.
These parameters are set when the simulation box is created.  Thus the
header only has the count of atoms, bonds, and so on.

Molecule files have a header followed by sections (just as in data
files), but the section names are different than those of a data file.
There is no "Atoms" section and the section formats in molecule files is
independent of the atom style.  Its information is split across multiple
sections, like "Coords", "Types", and "Charges".  Note that no "Masses"
section is needed here.  The atom masses are by default tied to the atom
type and set with a data file or the :doc:`mass command <mass>`.  A
"Masses" section would only be required for atom styles with per-atom
masses, e.g. atom style sphere, where in data files you would provide
the density and the diameter instead of the mass.

Since the entire file is a 'molecule', LAMMPS will assign a new
molecule-ID (if supported by the atom style) when atoms are instantiated
from a molecule file, e.g. with the :doc:`create_atoms command
<create_atoms>`.  It is possible to include a "Molecules" section to
indicate that the atoms belong to multiple 'molecules'.  Atom-IDs and
molecule-IDs in the molecule file are relative for the file
(i.e. starting from 1) and will be translated into actual atom-IDs also
when the atoms from the molecule are created.

.. code-block:: bash

   # Water molecule. SPC/E model.

   3 atoms
   2 bonds
   1 angles

   Coords

   1    1.12456   0.09298   1.27452
   2    1.53683   0.75606   1.89928
   3    0.49482   0.56390   0.65678

   Types

   1        1
   2        2
   3        2

   Charges

   1       -0.8472
   2        0.4236
   3        0.4236

   Bonds

   1   1      1      2
   2   1      1      3

   Angles

   1   1      2      1      3


There are also optional sections, e.g. about :doc:`SHAKE <fix_shake>`
and :doc:`special bonds <special_bonds>`. Those sections are only needed
if the molecule command is issued *before* the simulation box is
defined.  Otherwise, the molecule command can derive the required
settings internally.

Restart file
^^^^^^^^^^^^

LAMMPS restart files are binary files and not available in text format.
They can be identified by the first few bytes that contain the (C-style)
string ``LammpS RestartT`` as `magic string
<https://en.wikipedia.org/wiki/Magic_string>`_.  This string is followed
by a 16-bit integer of the number 1 used for detecting whether the
computer writing the restart has the same `endianness
<https://en.wikipedia.org/wiki/Endianness>`_ as the computer reading it.
If not, the file cannot be read correctly.  This integer is followed by
a 32-bit integer indicating the file format revision (currently 3),
which can be used to implement backward compatibility for reading older
revisions.

This information has been added to the `Unix "file" command's
<https://www.darwinsys.com/file/>` "magic" file so that restart files
can be identified without opening them.  If you have a fairly recent
version, it should already be included.  If you have an older version,
the LAMMPS source package :ref:`contains a file with the necessary
additions <magic>`.

The rest of the file is organized in sections of a 32-bit signed integer
constant indicating the kind of content and the corresponding value (or
values).  If those values are arrays (including C-style strings), then
the integer constant is followed by a 32-bit integer indicating the
length of the array.  This mechanism will read the data regardless of
the ordering of the sections.  Symbolic names of the section constants
are in the ``lmprestart.h`` header file.

LAMMPS restart files are not expected to be portable between platforms
or LAMMPS versions, but changes to the file format are rare.

.. Native Dump file
.. ^^^^^^^^^^^^^^^^
..
.. Potential files
.. ^^^^^^^^^^^^^^^
