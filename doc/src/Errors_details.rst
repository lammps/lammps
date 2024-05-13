Error and warning details
=========================

Many errors or warnings are self-explanatory and thus straightforward to
resolve.  However, there are also cases, where there is no single cause
and explanation, where LAMMPS can only detect symptoms of an error but
not the exact cause, or where the explanation needs to be more detailed than
what can be fit into a message printed by the program.  The following are
discussions of such cases.

.. _err0001:

Unknown identifier in data file
-------------------------------

This error happens when LAMMPS encounters a line of text with an
unexpected keyword while :doc:`reading a data file <read_data>`.  This
would be either header keywords or section header keywords.  This is
most commonly due to a mistyped keyword or due to a keyword that is
inconsistent with the :doc:`atom style <atom_style>` used.

The header section informs LAMMPS how many entries or lines are expected
in the various sections (like Atoms, Masses, Pair Coeffs, *etc.*\ ) of
the data file.  If there is a mismatch, LAMMPS will either keep reading
beyond the end of a section or stop reading before the section has
ended.  In that case the next line will not contain a recognized keyword.

Such a mismatch can also happen when the first line of the data
is *not* a comment as required by the format, but a line with a valid
header keyword.  That would result in LAMMPS expecting, for instance,
0 atoms because the "atoms" header line is the first line and thus
treated as a comment.

Another possibility to trigger this error is to have a keyword in the
data file that corresponds to a fix (e.g. :doc:`fix cmap <fix_cmap>`)
but the :doc:`read_data <read_data>` command is missing the (optional)
arguments that identify the fix and the header keyword and section
keyword or those arguments are inconsistent with the keywords in the
data file.

.. _err0002:

Incorrect format in ... section of data file
--------------------------------------------

This error happens when LAMMPS reads the contents of a section of a
:doc:`data file <read_data>` and the number of parameters in the line
differs from what is expected.  This most commonly happens, when the
atom style is different from what is expected for a specific data file
since changing the atom style usually changes the format of the line.

This error can also happen when the number of entries indicated in the
header of a data file (e.g. the number of atoms) is larger than the
number of lines provided (e.g. in the corresponding Atoms section)
and then LAMMPS will continue reading into the next section and that
would have a completely different format.
