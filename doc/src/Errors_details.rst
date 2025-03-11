Error and warning details
=========================

Many errors or warnings are self-explanatory and thus straightforward to
resolve.  However, there are also cases, where there is no single cause
and explanation, where LAMMPS can only detect symptoms of an error but
not the exact cause, or where the explanation needs to be more detailed
than what can be fit into a message printed by the program.  The
following are discussions of such cases.

- :ref:`Unknown identifier in data file <err0001>`
- :ref:`Incorrect format in ... section of data file <err0002>`
- :ref:`Illegal variable command: expected X arguments but found Y <err0003>`
- :ref:`Out of range atoms - cannot compute ... <err0004>`
- :ref:`Too many neighbor bins <err0009>`
- :ref:`Cannot use neighbor bins - box size \<\< cutoff <err0015>`
- :ref:`Domain too large for neighbor bins <err0017>`
- :ref:`Molecule topology/atom exceeds system topology/atom <err0024>`
- :ref:`Molecule topology type exceeds system topology type <err0025>`
- :ref:`Molecule attributes do not match system attributes <err0026>`

------

General troubleshooting advice
------------------------------

Below are suggestions that can help to understand the causes of problems
with simulations leading to errors or unexpected results.

Create a small test system
^^^^^^^^^^^^^^^^^^^^^^^^^^

Debugging problems often requires running a simulation many times with
small modifications, thus it can be a huge time saver to first assemble
a small test system input that has the same issue, but will take much
time until it triggers the error condition.  Also, it will be easier to
see what happens.

Visualize your trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^

To better understand what is causing problems, it is often very useful
to visualize the system close to the point of failure.  It may be
necessary to have LAMMPS output trajectory frames rather frequently.  To
avoid gigantic files, you can use :doc:`dump_modify delay <dump_modify>`
to delay output until the critical section is reached, and you can use a
smaller test system (see above).

Parallel versus serial
^^^^^^^^^^^^^^^^^^^^^^

Issues where something is "lost" or "missing" often exhibit that issue
only when running in parallel.  That doesn't mean there is no problem,
only the symptoms are not triggering an error quickly.  Correspondingly,
errors may be triggered faster with more processors and thus smaller
sub-domains.

Fast moving atoms
^^^^^^^^^^^^^^^^^

Fast moving atoms may be "lost" or "missing" when their velocity becomes
so large that they can cross a sub-domain within one timestep.  This
often happens when atoms are too close, but atoms may also "move" too
fast from sub-domain to sub-domain if the box changes rapidly, e.g. when
setting a large an initial box with :doc:`shrink-wrap boundary
conditions <boundary>` that collapses on the first step (in this case
the solution is often using 'm' instead of 's' as boundary condition).

To reduce the impact of "close contacts", one can remove those atoms or
molecules with something like :doc:`delete_atoms overlap 0.1 all all
<delete_atoms>`.  With periodic boundaries, a close contact pair of atoms
may be on opposite sides of the simulation box.  Another option would be
to first run a minimization (aka quench) before starting the MD.  Reducing
the time step can also help.  Many times, one just needs to "ease" the
system into a balanced state and can then switch to more aggressive settings.

The speed of atoms during an MD depends on the steepness of the
potential function and their mass.  Since the positions and velocities
of atoms are computed with finite timesteps, they choice of timestep can
be too large for a stable numeric integration of the trajectory.  In
those cases using (temporarily) :doc:`fix nve/limit <fix_nve_limit>` or
:doc:`fix dt/reset <fix_dt_reset>` can help to avoid too large updates
or adapt the timestep according to the displacements.


Pressure, forces, positions becoming NaN of Inf
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some potentials can overflow or have a division by zero with close contacts
or bad geometries (for the given force styles in use) leading to forces
that can no longer be represented as numbers.  Those will show as "NaN" or
"Inf".  On most machines, the program will continue, but there is no way
to recover from it and those NaN or Inf values will propagate.  So-called
:doc:`"soft-core" potentials <pair_fep_soft>` or the :doc:`"soft" repulsive-only
pair style <pair_soft>` are less prone for this behavior (depending on the
settings in use) and can be used at the beginning of a simulation.  Also,
single precision numbers can overflow much faster, so for the GPU or INTEL
package it may be beneficial to run with double precision initially before
switching to mixed or single precision for faster execution when the system
has relaxed.

Communication cutoff
^^^^^^^^^^^^^^^^^^^^

The communication cutoff determines the "overlap" between sub-domains
and atoms in these regions are referred to in LAMMPS as "ghost atoms".
This region has to be large enough to contain all atoms of a bond,
angle, dihedral or improper with just one atom in the actual sub-domain.
Typically, this cutoff is set to the largest cutoff from the :doc:`pair
style(s) <pair_style>` plus the :doc:`neighbor list skin distance
<neighbor>` and will be more than sufficient for all bonded
interactions.  But if the pair style cutoff is small this may bot be
enough.  LAMMPS will print a warning in this case using some heuristic
based on the equilibrium bond length, but that may not be sufficient for
cases where the force constants are small and thus bonds may be
stretched very far.  The communication cutoff can be adjusted with
:doc:`comm_modify cutoff \<value\> <comm_modify>`, but setting this too
large will waste CPU time and memory.

Neighbor list settings
^^^^^^^^^^^^^^^^^^^^^^

Every time LAMMPS rebuilds the neighbor lists, LAMMPS will also check
for "lost" or "missing" atoms.  Thus it can help to use very
conservative :doc:`neighbor list settings <neigh_modify>` and then
examine the neighbor list statistics if the neighbor list rebuild can be
safely delayed.  Rebuilding the neighbor list less frequently
(i.e. through increasing the *delay* or *every* setting has diminishing
returns and increasing risks).

Ignoring lost atoms
^^^^^^^^^^^^^^^^^^^

It is tempting to use the :doc:`thermo_modify lost ignore <thermo_modify>`
to avoid that LAMMPS stops with an error.  This setting should, however,
*only* be used when atoms *should* leave the system.  In general, ignoring
a problem does not solve it.

Units
^^^^^

A frequent cause for a variety of problems is due to using the wrong
:doc:`units <units>` settings for a particular potentials, especially
when reading them from a potential file.  Most of the (example)
potentials bundled with LAMMPS have a "UNITS:" tag that allows LAMMPS to
check of the units are consistent with what is intended, but potential
files from publications or potential parameter databases may lack this
metadata information and thus will not error out or warn when using the
wrong setting.  Most potential files usually use "metal" units, but some
are parameterized for other settings, most notably :doc:`ReaxFF
potentials <pair_reaxff>` that use "real" units.

Also, individual parameters for :doc:`pair_coeff <pair_coeff>` commands
taken from publications or other MD software, may need to be converted
and sometimes in unexpected ways.  Thus some careful checking is
recommended.

No error message printed
^^^^^^^^^^^^^^^^^^^^^^^^

In some cases - especially when running in parallel with MPI - LAMMPS
may stop without displaying an error.  But that does not mean, that
there was no error message, instead it is highly likely that the message
was written to a buffer and LAMMPS was aborted before the buffer was
output.  Usually, output buffers are output for every line of output,
but sometimes, this is delayed until 4096 or 8192 bytes of output have
been accumulated.  This buffering for screen and logfile output can be
disabled by using the :ref:`-nb or -nonbuf <nonbuf>` command-line flag.
This is most often needed when debugging crashing multi-replica
calculations.

------

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

.. _err0003:

Illegal variable command: expected X arguments but found Y
----------------------------------------------------------

This error indicates that there are the wrong number of arguments for a
specific variable command, but a common reason for that is a variable
expression that has whitespace but is not enclosed in single or double
quotes.

To explain, the LAMMPS input parser reads and processes lines.  The
resulting line is broken down into "words".  Those are usually
individual commands, labels, names, values separated by whitespace (a
space or tab character).  For "words" that may contain whitespace, they
have to be enclosed in single (') or double (") quotes.  The parser will
then remove the outermost pair of quotes and then pass that string as
"word" to the variable command.

Thus missing quotes or accidental extra whitespace will lead to the
error shown in the header because the unquoted whitespace will result
in the text being broken into more "words", i.e. the variable expression
being split.

.. _err0004:

Out of range atoms - cannot compute ...
---------------------------------------

The PPPM (and also PPPMDisp and MSM) methods require to assemble a grid
of electron density data derived from the (partial) charges assigned to
the atoms.  This charges are smeared out across multiple grid points
(see :doc:`kspace_modify order <kspace_modify>`).  When running in
parallel with MPI, LAMMPS uses a :doc:`domain decomposition scheme
<Developer_par_part>` where each processor manages a subset of atoms and
thus also a grid representing the density, which covers the actual
volume of the sub-domain and some extra space corresponding to the
:doc:`neighbor list skin <neighbor>`.  These are then :doc:`combined and
redistributed <Developer_par_long>` for parallel processing of the
long-range component of the Coulomb interaction.

The ``Out of range atoms`` error can happen, when atoms move too fast or
the neighbor list skin is too small or the neighbor lists are not
updated frequently enough.  Then the smeared charges cannot be fully
assigned to the density grid for all atoms.  LAMMPS checks for this
condition and stops with an error.  Most of the time, this is an
indication of a system with very high forces, most often at the
beginning of a simulation or when boundary conditions are changed.  The
error becomes more likely with more MPI processes.

There are multiple options to explore for avoiding the error.  The best
choice depends strongly on the individual system, and often a
combination of changes is required.  For example, more conservative MD
parameter settings can be used (larger neighbor skin, shorter time step,
more frequent neighbor list updates).  Sometimes, it helps to revisit
the system generation and avoid close contacts when building it, or use
the :doc:`delete_atoms overlap<delete_atoms>` command to delete those
close contact atoms, or run a minimization before the MD.  It can also
help to temporarily use a cutoff-Coulomb pair style and no kspace style
until the system has somewhat equilibrated and then switch to the
long-range solver.

.. _err0009:

Too many neighbor bins
----------------------

The simulation box has become too large relative to the size of a
neighbor bin and LAMMPS is unable to store the needed number of
bins. This typically implies the simulation box has expanded too far.
This can happen when some atoms move rapidly apart with shrink-wrap
boundaries or when a fix (like fix deform or a barostat) excessively
grows the simulation box.

.. _err0015:

Cannot use neighbor bins - box size \<\< cutoff
-----------------------------------------------

LAMMPS is unable to build neighbor bins since the size of the box is
much smaller than an interaction cutoff in at least one of its dimensions.
Typically, this error is triggered when the simulation box has one very
thin dimension. If a cubic neighbor bin had to fit exactly within
the thin dimension, then an inordinate amount of bins would be created to
fill space. This error can be avoided using the generally slower
:doc:`nsq neighbor style <neighbor>` or by increasing the size of the
smallest box lengths.

.. _err0017:

Domain too large for neighbor bins
----------------------------------

The domain has become extremely large so that neighbor bins cannot
be used. Too many neighbor bins would need to be created to fill space
Most likely, one or more atoms have been blown out of the simulation
box to a great distance or a fix (like fix deform or a barostat) has
excessively grown the simulation box.

.. _err0024:

Molecule topology/atom exceeds system topology/atom
---------------------------------------------------

LAMMPS uses :doc:`domain decomposition <Developer_par_part>` to
distribute data (i.e. atoms) across the MPI processes in parallel runs.
This includes topology data, that is data about bonds, angles,
dihedrals, impropers and :doc:`"special" neighbors <special_bonds>`.
This information is stored with either one or all atoms involved in such
a topology entry (which of the two option applies depends on the
:doc:`newton <newton>` setting for bonds. When reading a data file,
LAMMPS analyzes the requirements for this file and then the values
are "locked in" and cannot be extended.

So loading a molecule file that requires more of the topology per atom
storage or adding a data file with such needs will lead to an error.  To
avoid the error, one or more of the `extra/XXX/per/atom` keywords are
required to extend the corresponding storage.  It is no problem to
choose those numbers generously and have more storage reserved than
actually needed, but having these numbers set too small will lead to an
error.

.. _err0025:

Molecule topology type exceeds system topology type
---------------------------------------------------

The total number of atom, bond, angle, dihedral, and improper types is
"locked in" when LAMMPS creates the simulation box. This can happen
through either the :doc:`create_box <create_box>`, the :doc:`read_data
<read_data>`, or the :doc:`read_restart <read_restart>` command.  After
this it is not possible to refer to an additional type. So loading a
molecule file that uses additional types or adding a data file that
would require additional types will lead to an error.  To avoid the
error, one or more of the `extra/XXX/types` keywords are required to
extend the maximum number of the individual types.

.. _err0026:

Molecule attributes do not match system attributes
--------------------------------------------------

Choosing an :doc:`atom_style <atom_style>` in LAMMPS determines which
per-atom properties are available.  In a :doc:`molecule file
<molecule>`, however, it is possible to add sections (for example Masses
or Charges) that are not supported by the atom style.  Masses for
example, are usually not a per-atom property, but defined through the
atom type.  Thus it would not be required to have a Masses section and
the included data would be ignored.  LAMMPS prints this warning to
inform about this case.
