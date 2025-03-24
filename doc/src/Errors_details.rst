Error and warning details
=========================

Many errors or warnings are self-explanatory and thus straightforward to
resolve.  However, there are also cases, where there is no single cause
and explanation, where LAMMPS can only detect symptoms of an error but
not the exact cause, or where the explanation needs to be more detailed
than what can be fit into a message printed by the program.  The
following are discussions of such cases.

.. contents::

------

General troubleshooting advice
------------------------------

Below are suggestions that can help to understand the causes of problems
with simulations leading to errors or unexpected results.

.. _hint01:

Create a small test system
^^^^^^^^^^^^^^^^^^^^^^^^^^

Debugging problems often requires running a simulation many times with
small modifications, thus it can be a huge time saver to first assemble
a small test system input that has the same issue, but will take much
time until it triggers the error condition.  Also, it will be easier to
see what happens.

.. _hint02:

Visualize your trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^

To better understand what is causing problems, it is often very useful
to visualize the system close to the point of failure.  It may be
necessary to have LAMMPS output trajectory frames rather frequently.  To
avoid gigantic files, you can use :doc:`dump_modify delay <dump_modify>`
to delay output until the critical section is reached, and you can use a
smaller test system (see above).

.. _hint03:

Parallel versus serial
^^^^^^^^^^^^^^^^^^^^^^

Issues where something is "lost" or "missing" often exhibit that issue
only when running in parallel.  That doesn't mean there is no problem,
only the symptoms are not triggering an error quickly.  Correspondingly,
errors may be triggered faster with more processors and thus smaller
sub-domains.

.. _hint04:

Segmentation Fault
^^^^^^^^^^^^^^^^^^

A segmentation fault is an error reported by the **operating system**
and not LAMMPS itself.  It happens when a process tries to access a
memory address that is not available.  This can have **many** reasons:
memory has not been allocated, a memory buffer is not large enough, a
memory address is computed from an incorrect index, a memory buffer is
used after it has been freed, some general memory corruption.  When
investigating a segmentation fault (aka segfault), it is important to
determine which process is causing it; it may not always be LAMMPS.  For
example, some MPI library implementations report a segmentation fault
from their "mpirun" or "mpiexec" command when the application has been
terminated unexpectedly.

While a segmentation fault is likely an indication of a bug in LAMMPS,
it need not always be; it can also be the consequence of too aggressive
simulation settings.  For time critical code paths, LAMMPS will assume
the user has chosen the settings carefully and will not make any checks
to avoid to avoid performance penalties.

A crucial step in resolving a segmentation fault is to identify the exact
location in the code where it happens.  Please see `Errors_debug` for
a couple of examples showing how to do this on a Linux machine.  With
this information -- a simple way to reproduce the segmentation fault and
the exact :doc:`LAMMPS version <Manual_version>` and platform you are
running on -- you can contact the LAMMPS developers or post in the LAMMPS
forum to get assistance.

.. _hint05:

Fast moving atoms
^^^^^^^^^^^^^^^^^

Fast moving atoms may be "lost" or "missing" when their velocity becomes
so large that they can cross a sub-domain within one timestep.  This
often happens when atoms are too close, but atoms may also "move" too
fast from sub-domain to sub-domain if the box changes rapidly. E.g. when
setting a large an initial box with :doc:`shrink-wrap boundary
conditions <boundary>` that collapses on the first step (in this case
the solution is often using 'm' instead of 's' as a boundary condition).

To reduce the impact of "close contacts", one can remove those atoms or
molecules with something like :doc:`delete_atoms overlap 0.1 all all
<delete_atoms>`.  With periodic boundaries, a close contact pair of atoms
may be on opposite sides of the simulation box.  Another option would be
to first run a minimization (aka quench) before starting the MD.  Reducing
the time step can also help.  Many times, one just needs to "ease" the
system into a balanced state and can then switch to more aggressive settings.

The speed of atoms during an MD run depends on the steepness of the
potential function and their mass.  Since the positions and velocities
of atoms are computed with finite timesteps, the timestep needs to be
small enough for stable numeric integration of the trajectory.  If the timestep
is too large during initialization (or other instances of extreme dynamics),
using :doc:`fix nve/limit <fix_nve_limit>` or :doc:`fix dt/reset <fix_dt_reset>`
temporarily can help to avoid too large updates or adapt the timestep according
to the displacements.

.. _hint06:

Ignoring lost atoms
^^^^^^^^^^^^^^^^^^^

It is tempting to use the :doc:`thermo_modify lost ignore <thermo_modify>`
to avoid LAMMPS aborting with an error on lost atoms.  This setting should,
however, *only* be used when atoms *should* leave the system.  In general,
ignoring a problem does not solve it.

.. _hint07:

Pressure, forces, positions becoming NaN or Inf
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

.. _hint08:

Communication cutoff
^^^^^^^^^^^^^^^^^^^^

The communication cutoff determines the "overlap" between sub-domains
and atoms in these regions are referred to in LAMMPS as "ghost atoms".
This region has to be large enough to contain all atoms of a bond,
angle, dihedral, or improper with just one atom in the actual sub-domain.
Typically, this cutoff is set to the largest cutoff from the :doc:`pair
style(s) <pair_style>` plus the :doc:`neighbor list skin distance
<neighbor>` and will typically be sufficient for all bonded
interactions.  But if the pair style cutoff is small, this may not be
enough.  LAMMPS will print a warning in this case using some heuristic
based on the equilibrium bond length, but that still may not be sufficient
for cases where the force constants are small and thus bonds may be
stretched very far.  The communication cutoff can be adjusted with
:doc:`comm_modify cutoff \<value\> <comm_modify>`, but setting this too
large will waste CPU time and memory.

.. _hint09:

Neighbor list settings
^^^^^^^^^^^^^^^^^^^^^^

Every time LAMMPS rebuilds the neighbor lists, LAMMPS will also check
for "lost" or "missing" atoms.  Thus it can help to use very
conservative :doc:`neighbor list settings <neigh_modify>` and then
examine the neighbor list statistics if the neighbor list rebuild can be
safely delayed.  Rebuilding the neighbor list less frequently
(i.e. through increasing the *delay* or *every*) setting has diminishing
returns and increasing risks.

.. _hint10:

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
taken from publications or other MD software may need to be converted
and sometimes in unexpected ways.  Thus some careful checking is
recommended.

.. _hint11:

No error message printed
^^^^^^^^^^^^^^^^^^^^^^^^

In some cases -- especially when running in parallel with MPI -- LAMMPS
may stop without displaying an error.  But the fact that nothing was
displayed does not mean there was not an error message. Instead it is
highly likely that the message was written to a buffer and LAMMPS was
aborted before the buffer was output.  Usually, output buffers are output
for every line of output, but sometimes this is delayed until 4096 or
8192 bytes of output have been accumulated.  This buffering for screen
and logfile output can be disabled by using the :ref:`-nb or -nonbuf
<nonbuf>` command-line flag. This is most often needed when debugging
crashing multi-replica calculations.

.. _hint12:

Errors before or after the simulation box is created
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As critical step in a LAMMPS input is when the simulation box is
defined, either with a :doc:`create_box command <create_box>`, a
:doc:`read_data command <read_data>`, or a :doc:`read_restart command
<read_restart>`.  After this step, certain settings are locked in (e.g.
units, or number of atom, bond, angle, dihedral, improper types) and
cannot be changed after that.  Consequently, commands that change such
settings (e.g. :doc:`units <units>`) are only allowed before the box is
defined.  Very few commands can be used before and after, like
:doc:`pair_style <pair_style>` (but not :doc:`pair_coeff <pair_coeff>`).
Most LAMMPS commands must be used after the simulation box is created.

Consequently, LAMMPS will stop with an error, if a command is used in
the wrong place.  This is not always obvious. So index or string style
:doc:`variables <variable>` can be expanded anywhere in the input, but
equal style (or similar) variables can only be expanded before the box
is defined if they do not reference anything that cannot be defined
before the box (e.g. a compute or fix reference or a thermo keyword).

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
arguments that identify the fix and its header and section keywords.
Alternatively, those arguments are inconsistent with the keywords in the
data file.

.. _err0002:

Incorrect format in ... section of data file
--------------------------------------------

This error happens when LAMMPS reads the contents of a section of a
:doc:`data file <read_data>` and the number of parameters in the line
differs from what is expected.  This most commonly happens when the
atom style is different from what is expected for a specific data file
since changing the atom style usually changes the format of the line.

This error can also occur when the number of entries indicated in the
header of a data file (e.g. the number of atoms) is larger than the
number of lines provided (e.g. in the corresponding Atoms section)
causing LAMMPS to continue reading into the next section which has
a completely different format.

.. _err0003:

Illegal variable command: expected X arguments but found Y
----------------------------------------------------------

This error indicates that a variable command has the wrong number of
arguments. A common reason for this is that the variable expression
has whitespace, but is not enclosed in single or double quotes.

To explain, the LAMMPS input parser reads and processes lines.  The
resulting line is broken down into "words".  Those are usually
individual commands, labels, names, and values separated by whitespace (a
space or tab character).  For "words" that may contain whitespace, they
have to be enclosed in single (') or double (") quotes.  The parser will
then remove the outermost pair of quotes and pass that string as
"word" to the variable command.

Thus missing quotes or accidental extra whitespace will trigger this
error because the unquoted whitespace will result in the text being broken
into more "words", i.e. the variable expression being split.

.. _err0004:

Out of range atoms - cannot compute ...
---------------------------------------

The PPPM (and also PPPMDisp and MSM) methods need to assemble a grid
of electron density data derived from the (partial) charges assigned to
the atoms.  These charges are smeared out across multiple grid points
(see :doc:`kspace_modify order <kspace_modify>`).  When running in
parallel with MPI, LAMMPS uses a :doc:`domain decomposition scheme
<Developer_par_part>` where each processor manages a subset of atoms and
thus also a grid representing the density. The processor's grid covers the
actual volume of the sub-domain and some extra space corresponding to the
:doc:`neighbor list skin <neighbor>`.  These are then :doc:`combined and
redistributed <Developer_par_long>` for parallel processing of the
long-range component of the Coulomb interaction.

The ``Out of range atoms`` error can happen when atoms move too fast,
the neighbor list skin is too small, or the neighbor lists are not
updated frequently enough.  The smeared charges cannot then be fully
assigned to the density grid for all atoms.  LAMMPS checks for this
condition and stops with an error.  Most of the time, this is an
indication of a system with very high forces, often at the beginning
of a simulation or when boundary conditions are changed.  The
error becomes more likely with more MPI processes.

There are multiple options to explore for avoiding the error.  The best
choice depends strongly on the individual system, and often a
combination of changes is required.  For example, more conservative MD
parameter settings can be used (larger neighbor skin, shorter time step,
more frequent neighbor list updates).  Sometimes, it helps to revisit
the system generation and avoid close contacts when building it. Otherwise
one can use the :doc:`delete_atoms overlap<delete_atoms>` command to delete
those close contact atoms or run a minimization before the MD.  It can also
help to temporarily use a cutoff-Coulomb pair style and no kspace style
until the system has somewhat equilibrated and then switch to the
long-range solver.

.. _err0005:

Bond (or angle, dihedral, or improper) atoms missing
----------------------------------------------------

The second atom needed to compute a particular bond (or the third or fourth
atom for angle, dihedral, or improper) is missing on the indicated timestep
and processor. Typically, this is because the two bonded atoms have become
too far apart relative to the communication cutoff distance for ghost atoms.
By default, the communication cutoff is set by the pair cutoff. However, to
accommodate larger distances between topologically connected atoms, it can
be manually adjusted using :doc:`comm_modify <comm_modify>` at the cost of
increased communication and more ghost atoms. However, missing bond atoms
may also indicate that there are unstable dynamics which caused the atoms
to blow apart. In this scenario, increasing the communication distance will
not solve the underlying issue. Rather, see :ref:`Fast moving atoms <hint05>`
and :ref:`Neighbor list settings <hint09>` in the general troubleshooting
section above for ideas to fix unstable dynamics.

If atoms are intended to be lost during a simulation (e.g. due to open boundary
conditions or :doc:`fix evaporate <fix_evaporate>`) such that two bonded atoms
may be lost at different times from each other, this error can be converted to a
warning or turned off using the *lost/bond* keyword in the :doc:`thermo_modify
<thermo_modify>` command.

.. _err0006:

Non-numeric atom coords - simulation unstable
---------------------------------------------
This error usually occurs due to issues with system geometry or the potential in
use. See :ref:`Pressure, forces, positions becoming NaN or Inf <hint07>` above in the
general troubleshooting section.

.. _err0007:

Non-numeric pressure - simulation unstable
------------------------------------------
This error usually occurs due to issues with system geometry or the potential in
use. See :ref:`Pressure, forces, positions becoming NaN or Inf <hint07>` above in the
general troubleshooting section.


.. _err0008:

Lost atoms ...
--------------

A simulation stopping with an error due to lost atoms can have multiple
causes. In the majority of cases, lost atoms are unexpected and a result
of extremely high velocities causing instabilities in the system, and
those velocities can result from a variety of issues. For ideas on how
to track down issues with unexpected lost atoms, see :ref:`Fast moving
atoms <hint05>` and :ref:`Neighbor list settings <hint09>` in the
general troubleshooting section above. In specific situations however,
losing atoms is expected material behavior (e.g. with sputtering and
surface evaporation simulations) and an unwanted crash can be resolved
by changing the :doc:`thermo_modify lost <thermo_modify>` keyword from
the default 'error' to 'warn' or 'ignore' (though heed the advice in
:ref:`Ignoring lost atoms <hint06>` above!).

.. _err0009:

Too many neighbor bins
----------------------

The simulation box has become too large relative to the size of a
neighbor bin and LAMMPS is unable to store the needed number of
bins. This typically implies the simulation box has expanded too far.
This can happen when some atoms move rapidly apart with shrink-wrap boundaries
or when a fix (like fix deform or a barostat) excessively grows the simulation
box.

.. _err0010:

Unrecognized pair style ... is part of ... package which is not enabled in this LAMMPS binary
---------------------------------------------------------------------------------------------

The LAMMPS executable (binary) being used was not compiled with a package
containing the specified pair style. This indicates that the executable needs to
be re-built after enabling the correct package in the relevant Makefile or CMake
build directory. See :doc:`Section 3. Build LAMMPS <Build>` for more details.
One can check if the expected package and pair style is present in the
executable by running it with the ``-help`` (or ``-h``) flag on the command
line. One common oversight, especially for beginner LAMMPS users, is to enable
the package, but to forget to run commands to rebuild (e.g., to run the final
``make`` or ``cmake`` command).

If this error is occurring with an executable that the user does not control
(e.g., through a module on HPC clusters), the user will need to get in contact
with the relevant person or people who can update the executable.

.. _err0012:

fmt::format_error
-----------------

LAMMPS uses the `{fmt} library <https://fmt.dev>`_ for advanced string
formatting tasks.  This is similar to the ``printf()`` family of
functions from the standard C library, but more flexible.  If there is a
bug in the LAMMPS code and the format string does not match the list of
arguments or has some other error, this error message will be shown.
You should contact the LAMMPS developers and report the bug as a `GitHub
Bug Report Issue <https://github.com/lammps/lammps/issues>`_ along with
sufficient information to easily reproduce it.


.. _err0013:

Substitution for illegal variable
---------------------------------

A variable in an input script or a variable expression was not found in
the list of valid variables.  The most common reason for this is a typo
somewhere in the input file such that the expression uses an invalid variable
name.  The second most common reason is omitting the curly braces for a
direct variable with a name that is not a single letter.  For example:

.. code-block:: LAMMPS

   variable cutoff index 10.0
   pair_style lj/cut ${cutoff}  # this is correct
   pair_style lj/cut $cutoff    # this is incorrect, LAMMPS looks for 'c' instead of 'cutoff'
   variable c      index 5.0    # if $c is defined, LAMMPS subsitutes only '$c' and reads: 5utoff

Another potential source of this error may be invalid command line
variables (-var or -v argument) used when launching LAMMPS from an
interactive shell or shell scripts.  An uncommon source for this error
is using the :doc:`next command <next>` to advance through a list of values
provided by an index style variable.  If there is no remaining element in
the list, LAMMPS will delete the variable and any following expansion or
reference attempt will trigger the error.

Users with harder-to-track variable errors might also find reading
:doc:`Section 5.2. Parsing rules for input scripts<Commands_parse>`
helpful.

.. _err0014:

Bond atom missing in image check or box size check
--------------------------------------------------

This can be either an error or a warning depending on your
:doc:`thermo_modify settings <thermo_modify>`.  It is flagged in a part
of the LAMMPS code where it updates the domain decomposition and before
it builds the neighbor lists.  It checks that both atoms of a bond are
within the communication cutoff of a subdomain.  It is usually caused by
atoms moving too fast (see the :ref:`paragraph on fast moving atoms
<hint05>`), or by the :doc:`communication cutoff being too
small <comm_modify>`, or by waiting too long between :doc:`sub-domain
and neighbor list updates <neigh_modify>`.

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

.. _err0016:

Did not assign all atoms correctly
----------------------------------

This error happens most commonly when :doc:`reading a data file <read_data>`
under :doc:`non-periodic boundary conditions<boundary>`.  Only atoms with
positions **inside** the simulation box will be read and thus any atoms
outside the box will be skipped and the total atom count will not match,
which triggers the error.  This does not happen with periodic boundary
conditions where atoms outside the principal box will be "wrapped" into
the principal box and their image flags set accordingly.

Similar errors can happen with the :doc:`replicate command<replicate>` or
the :doc:`read_restart command<read_restart>`.  In these cases the cause
may be a problematic geometry, an insufficient communication cutoff, or
a bug in the LAMMPS source code.  In these cases it is advisable to set
up :ref:`small test case <hint01>` for testing and debugging.  This will
be required in case you need to get help from a LAMMPS developer.

.. _err0017:

Domain too large for neighbor bins
----------------------------------

The domain has become extremely large so that neighbor bins cannot
be used. Too many neighbor bins would need to be created to fill space.
Most likely, one or more atoms have been blown a great distance out of
the simulation box or a fix (like fix deform or a barostat) has
excessively grown the simulation box.

.. _err0018:

Step X: (h)bondchk failed
-------------------------

This error is a consequence of the heuristic memory allocations for
buffers of the regular ReaxFF version.  In ReaxFF simulations, the lists
of bonds and hydrogen bonds can change due to chemical reactions.  The
default approach, however, assumes that these changes are not very
large, so it allocates buffers for the current system setup plus a
safety margin.  This can be adjusted with the :doc:`safezone, mincap,
and minhbonds settings of the pair style <pair_reaxff>`, but only to some
extent.  When equilibrating a new system, or simulating a sparse system
in parallel, this can be difficult to control and become wasteful.  A
simple workaround is often to break a simulation down in multiple
chunks.  A better approach, however, is to compile and use the KOKKOS
package version of ReaxFF (you do not need a GPU for that, but can also
compile it in serial or OpenMP mode), which uses a more robust
memory allocation approach.

.. _err0019:

Numeric index X is out of bounds
--------------------------------

This error most commonly happens when setting force field coefficients
with either the :doc:`pair_coeff <pair_coeff>`, the :doc:`bond_coeff
<bond_coeff>`, the :doc:`angle_coeff <angle_coeff>`, the
:doc:`dihedral_coeff <dihedral_coeff>`, or the :doc:`improper_coeff
<improper_coeff>` command.  These commands accept type labels,
explicit numbers, and wildcards for ranges of numbers.  If the numeric
value of any of these is outside the valid range (defined by the number
of corresponding types), LAMMPS will stop with this error.  A few other
commands and styles also allow ranges of numbers and check
using the same method and thus print the same kind of error.

The cause is almost always a typo in the input or a logic error
when defining the values or ranges.  So one needs to carefully
review the input.  Along with the error, LAMMPS will print the
valid range as a hint.

.. _err0020:

Compute, fix, or variable vector or array is accessed out-of-range
------------------------------------------------------------------

When accessing an individual element of a global vector or array or a
per-atom vector or array provided by a compute or fix or atom-style or
vector-style variable or data from a specific atom, an index in square
brackets ("[ ]") (or two indices) must be provided to determine which
element to access and it must be in a valid range or else LAMMPS would
access invalid data or crash with a segmentation fault. In the two most
common cases, where this data is accessed, :doc:`variable expressions
<variable>` and :doc:`thermodynamic output <thermo_style>`, LAMMPS will
check for valid indices and stop with an error otherwise.

While LAMMPS is written in C++ (which uses 0 based indexing) these
indices start at 1 (i.e. similar to Fortran).  Any index smaller than 1
or larger than the maximum allowed value should trigger this error.
Since this kind of error frequently happens with rather complex
expressions, it is recommended to test these with small test systems,
where the values can be tracked with output files for all relevant
properties at every step.

.. _err0021:

Incorrect args for pair coefficients (also bond/angle/dihedral/improper coefficients)
-------------------------------------------------------------------------------------

The parameters in the :doc:`pair_coeff <pair_coeff>` command for a specified
:doc:`pair_style <pair_style>` have a missing or erroneous argument. The same
applies when seeing this error for :doc:`bond_coeff <bond_coeff>`,
:doc:`angle_coeff <angle_coeff>`,  :doc:`dihedral_coeff <dihedral_coeff>`, or
:doc:`improper_coeff <improper_coeff>` and their respective style commands when
using the MOLECULE or EXTRA-MOLECULE packages. The cases below describe
some ways to approach pair coefficient errors, but the same strategies
apply to bonded systems as well.

Outside of normal typos, this error can have several sources. In all cases, the
first step is to compare the command arguments to the expected format found in
the corresponding :doc:`pair_style <pair_style>` page. This can reveal cases
where, for example, a pair style was changed, but the pair coefficients were not
updated. This can happen especially with pair style variants such as
:doc:`pair_style eam <pair_eam>` vs. :doc:`pair_style eam/alloy <pair_style>`
that look very similar but accept different parameters (the latter 'eam/alloy'
variant takes element type names while 'eam' does not).

Another common source of coefficient errors is when using multiple pair styles
with commands such as :doc:`pair_style hybrid <pair_hybrid>`. Using hybrid pair
styles requires adding an extra "label" argument in the coefficient commands
that designates which pair style the command line refers to. Moreover, if
the same pair style is used multiple times, this label must be followed by
an additional numeric argument. Also, different pair styles may require
different arguments.

This error message might also require a close look at other LAMMPS input files
that are read in by the input script, such as data files or restart files.

.. _err0022:

Energy was not tallied on needed timestep (also virial, per-atom energy, per-atom virial)
-----------------------------------------------------------------------------------------

This error is generated when LAMMPS attempts to access an out-of-date or
non-existent energy, pressure, or virial.  For efficiency reasons,
LAMMPS does *not* calculate these quantities when the forces are
calculated on every timestep or iteration.  Global quantities are only
calculated when they are needed for :doc:`thermo <thermo_style>` output
(at the beginning, end, and at regular intervals specified by the
:doc:`thermo <thermo>` command). Similarly, per-atom quantities are only
calculated if they are needed to write per-atom energy or virial to a
dump file.  This system works fine for simple input scripts.  However,
the many user-specified `variable`, `fix`, and `compute` commands that
LAMMPS provides make it difficult to anticipate when a quantity will be
requested. In some use cases, LAMMPS will figure out that a quantity is
needed and arrange for it to be calculated on that timestep e.g. if it
is requested by :doc:`fix ave/time <fix_ave_time>` or similar commands.
If that fails, it can be detected by a mismatch between the current
timestep and when a quantity was last calculated, in which case an error
message of this type is generated.

The most common cause of this type of error is requesting a quantity before
the start of the simulation.

.. code-block:: LAMMPS

   # run 0 post no               # this will fix the error
   variable e equal pe           # requesting energy compute
   print "Potential energy = $e" # this will generate the error
   run 1000                      # start of simulation

This situation can be avoided by adding in a "run 0" command, as explained in
more detail in the "Variable Accuracy" section of the
:doc:`variable <variable>` doc page.

Another cause is requesting a quantity on a timestep that is not
a thermo or dump output timestep. This can often be
remedied by increasing the frequency of thermo or dump output.

.. _err0023:

Molecule auto special bond generation overflow
----------------------------------------------

In order to correctly apply the :doc:`special_bonds <special_bonds>`
settings (also known as "exclusions"), LAMMPS needs to maintain for each
atom a list of atoms that are connected to this atom, either directly with
a bond or indirectly through bonding with an intermediate atom(s). The purpose
is to either remove or tag those pairs of atoms in the neighbor list.  This
information is stored with individual
atoms and thus the maximum number of such "special" neighbors is set
when the simulation box is created.  When reading (relative) geometry
and topology of a 'molecule' from a :doc:`molecule file <molecule>`,
LAMMPS will build the list of such "special" neighbors for the molecule atom
(if not given in the molecule file explicitly).  The error is triggered
when the resulting list is too long for the space reserved when
creating the simulation box.  The solution is to increase the
corresponding setting.  Overestimating this value will only consume
more memory, and is thus a safe choice.

.. _err0024:

Molecule topology/atom exceeds system topology/atom
---------------------------------------------------

LAMMPS uses :doc:`domain decomposition <Developer_par_part>` to
distribute data (i.e. atoms) across the MPI processes in parallel
runs. This includes topology data about bonds, angles,
dihedrals, impropers and :doc:`"special" neighbors <special_bonds>`.
This information is stored with either one or all atoms involved in such
a topology entry (which of the two option applies depends on the
:doc:`newton <newton>` setting for bonds). When reading a data file,
LAMMPS analyzes the requirements for this file and then the values are
"locked in" and cannot be extended.

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

.. _err0027:

Inconsistent image flags
------------------------

This warning happens when the distance between the *unwrapped* x-, y-,
or z-components of the coordinates of a bond is larger than half the box
with periodic boundaries or larger than the box with non-periodic
boundaries.  It means that the positions and image flags have become
inconsistent.  LAMMPS will still compute bonded interactions based on
the closest periodic images of the atoms and thus in most cases the
results will be correct.  Nevertheless, it is good practice to update
the system so that the message does not appear.  It will help with
future manipulations of the system.

There is one case where this warning *must* appear: when you have a
chain of connected bonds that pass through the entire box and connect
back to the first atom in the chain through periodic boundaries,
i.e. some kind of "infinite polymer".  In that case, the bond image
flags *must* be inconsistent for the one bond that reaches back to the
beginning of the chain.


.. _err0028:

No fixes with time integration, atoms won't move
------------------------------------------------

This warning will be issued if LAMMPS encounters a :doc:`run <run>` command that
does not have a preceding :doc:`fix <fix>` command that updates atom/object
positions and velocities per step. In other words, there are no fixes detected
that perform velocity-Verlet time integration, such as :doc:`fix nve <fix_nve>`.
Note that this alert does not mean that there are no active fixes. LAMMPS has a
very wide variety of fixes, many of which do not move objects but also operate
through steps, such as printing outputs (e.g. :doc:`fix print <fix_print>`),
performing calculations (e.g. :doc:`fix ave/time <fix_ave_time>`), or changing
other system parameters (e.g. :doc:`fix dt/reset <fix_dt_reset>`). It is up to
the user to determine whether the lack of a time-integrating fix is intentional
or not.


.. _err0029:

System is not charge neutral, net charge = ...
----------------------------------------------

the sum of charges in the system is not zero. When a system is not
charge-neutral, methods that evolve/manipulate per-atom charges, evaluate
Coulomb interactions, evaluate Coulomb forces, or evaluate/manipulate other
properties relying on per-atom charges may raise this warning. A non-zero
net charge most commonly arises after setting per-atom charges :doc:`set <set>`
such that the sum is non-zero or by reading in a system through :doc:`read_data
<read_data>` where the per-atom charges do not sum to zero. However, a loss of
charge neutrality may occur in other less common ways, like when charge
equilibration methods (e.g., :doc:`fix qeq <fix_qeq>`) fail.

A similar warning/error may be raised when using certain charge equilibration
methods: :doc:`fix qeq <fix_qeq>`, :doc:`fix qeq/comb <fix_qeq_comb>`, :doc:`fix
qeq/reaxff <fix_qeq_reaxff>`, and :doc:`fix qtpie/reaxff <fix_qtpie_reaxff>`. In
such cases, this warning/error will be raised for the fix :doc:`group <group>`
when the group has a non-zero net charge.

When the system is expected to be charge-neutral, this warning often arises due
to an error in the lammps input (e.g., an incorrect :doc:`set <set>` command,
error in the data file read by :doc:`read_data <read_data>`, incorrectly
grouping atoms with charge, etc.). If the system is NOT expected to be
charge-neutral, the user should make sure that the method(s) used are
appropriate for systems with a non-zero net charge. Some commonly used fixes for
charge equilibration :doc:`fix qeq <fix_qeq>`, pair styles that include charge
interactions :doc:`pair_style coul/XXX <pair_coul>`, and kspace methods
:doc:`kspace_style <kspace_style>` can, in theory, support systems with non-zero
net charge. However, non-zero net charge can lead to spurious artifacts. The
severity of these artifacts depends on the magnitude of total charge, system
size, and methods used. Before running simulations or calculations for systems
with non-zero net charge, users should test for artifacts and convergence of
properties.

.. _err0030:

Variable evaluation before simulation box is defined
----------------------------------------------------

This error happens, when trying to expand or use an equal- or atom-style
variable (or an equivalent style), where the expression contains a
reference to something (e.g. a compute reference, a property of an atom,
or a thermo keyword) that is not allowed to be used before the
simulation box is defined.  See the paragraph on :ref:`errors before or
after the simulation box is created <hint12>` for additional
information.

.. _err0031:

Invalid thermo keyword 'X' in variable formula
----------------------------------------------

This error message is often misleading.  It is caused when evaluating a
:doc:`variable command <variable>` expression and LAMMPS comes across a
string that it does not recognize.  LAMMPS first checks if a string is a
reference to a compute, fix, custom property, or another variable by
looking at the first 2-3 characters (and if it is, it checks whether the
referenced item exists).  Next LAMMPS checks if the string matches one
of the available functions or constants.  If that fails, LAMMPS will
assume that this string is a :doc:`thermo keyword <thermo_style>` and
let the code for printing thermodynamic output return the corresponding
value.  However, if this fails too, since the string is not a thermo
keyword, LAMMPS stops with the 'Invalid thermo keyword' error.  But it
is also possible, that there is just a typo in the name of a valid
variable function.  Thus it is recommended to check the failing variable
expression very carefully.
