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

This error happens when LAMMPS encounters a line of text in an unexpected format
while reading a data file. This is most commonly cause by inconsistent header and
section data.  The header section informs LAMMPS how many entries or lines are expected in the
various sections (like Atoms, Masses, Pair Coeffs, *etc.*\ ) of the data file.
If there is a mismatch, LAMMPS will either keep reading beyond the end of a section
or stop reading before the section has ended.

Such a mismatch can happen unexpectedly when the first line of the data
is *not* a comment as required by the format.  That would result in
LAMMPS expecting, for instance, 0 atoms because the "atoms" header line
is treated as a comment.

