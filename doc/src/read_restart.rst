.. index:: read_restart

read_restart command
====================

Syntax
""""""

.. code-block:: LAMMPS

   read_restart file

* file = name of binary restart file to read in

Examples
""""""""

.. code-block:: LAMMPS

   read_restart save.10000
   read_restart restart.*

Description
"""""""""""

Read in a previously saved system configuration from a restart file.
This allows continuation of a previous run.  Details about what
information is stored (and not stored) in a restart file is given below.
Basically this operation will re-create the simulation box with all its
atoms and their attributes as well as some related global settings, at
the point in time it was written to the restart file by a previous
simulation.  The simulation box will be partitioned into a regular 3d
grid of rectangular bricks, one per processor, based on the number of
processors in the current simulation and the settings of the
:doc:`processors <processors>` command.  The partitioning can later be
changed by the :doc:`balance <balance>` or :doc:`fix balance
<fix_balance>` commands.

.. deprecated:: 23Jun2022

Atom coordinates that are found to be outside the simulation box when
reading the restart will be remapped back into the box and their image
flags updated accordingly.  This previously required specifying the
*remap* option, but that is no longer required.

Restart files are saved in binary format to enable exact restarts,
meaning that the trajectories of a restarted run will precisely match
those produced by the original run had it continued on.

The binary restart file format was not designed with backward, forward,
or cross-platform compatibility in mind, so the files are only expected
to be read correctly by the same LAMMPS executable on the same platform.
Changes to the architecture, compilation settings, or LAMMPS version can
render a restart file unreadable or it may read the data incorrectly.
If you want a more portable format, you can use the data file format as
created by the :doc:`write_data <write_data>` command.  Binary restart
files can also be converted into a data file from the command line by
the LAMMPS executable that wrote them using the :ref:`-restart2data
<restart2data>` command line flag.

Several things can prevent exact restarts due to round-off effects, in
which case the trajectories in the 2 runs will slowly diverge.  These
include running on a different number of processors or changing
certain settings such as those set by the :doc:`newton <newton>` or
:doc:`processors <processors>` commands.  LAMMPS will issue a warning in
these cases.

Certain fixes will not restart exactly, though they should provide
statistically similar results.  These include :doc:`fix shake
<fix_shake>` and :doc:`fix langevin <fix_langevin>`.

Certain pair styles will not restart exactly, though they should
provide statistically similar results.  This is because the forces
they compute depend on atom velocities, which are used at half-step
values every timestep when forces are computed.  When a run restarts,
forces are initially evaluated with a full-step velocity, which is
different than if the run had continued.  These pair styles include
:doc:`granular pair styles <pair_gran>`, :doc:`pair dpd <pair_dpd>`, and
:doc:`pair lubricate <pair_lubricate>`.

If a restarted run is immediately different than the run which
produced the restart file, it could be a LAMMPS bug, so consider
:doc:`reporting it <Errors_bugs>` if you think the behavior is a bug.

Because restart files are binary, they may not be portable to other
machines.  In this case, you can use the :doc:`-restart command-line
switch <Run_options>` to convert a restart file to a data file.

Similar to how restart files are written (see the :doc:`write_restart
<write_restart>` and :doc:`restart <restart>` commands), the restart
filename can contain two wild-card characters.  If a "\*" appears in the
filename, the directory is searched for all filenames that match the
pattern where "\*" is replaced with a timestep value.  The file with the
largest timestep value is read in.  Thus, this effectively means, read
the latest restart file.  It's useful if you want your script to
continue a run from where it left off.  See the :doc:`run <run>` command
and its "upto" option for how to specify the run command so it does not
need to be changed either.

If a "%" character appears in the restart filename, LAMMPS expects a
set of multiple files to exist.  The :doc:`restart <restart>` and
:doc:`write_restart <write_restart>` commands explain how such sets are
created.  Read_restart will first read a filename where "%" is
replaced by "base".  This file tells LAMMPS how many processors
created the set and how many files are in it.  Read_restart then reads
the additional files.  For example, if the restart file was specified
as save.% when it was written, then read_restart reads the files
save.base, save.0, save.1, ... save.P-1, where P is the number of
processors that created the restart file.

Note that P could be the total number of processors in the previous
simulation, or some subset of those processors, if the *fileper* or
*nfile* options were used when the restart file was written; see the
:doc:`restart <restart>` and :doc:`write_restart <write_restart>` commands
for details.  The processors in the current LAMMPS simulation share
the work of reading these files; each reads a roughly equal subset of
the files.  The number of processors which created the set can be
different the number of processors in the current LAMMPS simulation.
This can be a fast mode of input on parallel machines that support
parallel I/O.

----------

Here is the list of information included in a restart file, which
means these quantities do not need to be re-specified in the input
script that reads the restart file, though you can redefine many of
these settings after the restart file is read.

* :doc:`units <units>`
* :doc:`newton bond <newton>` (see discussion of newton command below)
* :doc:`atom style <atom_style>` and :doc:`atom_modify <atom_modify>` settings id, map, sort
* :doc:`comm style <comm_style>` and :doc:`comm_modify <comm_modify>` settings mode, cutoff, vel
* :doc:`timestep <timestep>`
* simulation box size and shape and :doc:`boundary <boundary>` settings
* atom :doc:`group <group>` definitions
* per-type atom settings such as :doc:`mass <mass>`
* per-atom attributes including their group assignments and molecular topology attributes (bonds, angles, etc)
* force field styles (:doc:`pair <pair_style>`, :doc:`bond <bond_style>`, :doc:`angle <angle_style>`, etc)
* force field coefficients (:doc:`pair <pair_coeff>`, :doc:`bond <bond_coeff>`, :doc:`angle <angle_coeff>`, etc) in some cases (see below)
* :doc:`pair_modify <pair_modify>` settings, except the compute option
* :doc:`special_bonds <special_bonds>` settings

Here is a list of information not stored in a restart file, which
means you must re-issue these commands in your input script, after
reading the restart file.

* :doc:`newton pair <newton>` (see discussion of newton command below)
* :doc:`fix <fix>` commands (see below)
* :doc:`compute <compute>` commands (see below)
* :doc:`variable <variable>` commands
* :doc:`region <region>` commands
* :doc:`neighbor list <neighbor>` criteria including :doc:`neigh_modify <neigh_modify>` settings
* :doc:`kspace_style <kspace_style>` and :doc:`kspace_modify <kspace_modify>` settings
* info for :doc:`thermodynamic <thermo_style>`, :doc:`dump <dump>`, or :doc:`restart <restart>` output

The :doc:`newton <newton>` command has two settings, one for pairwise
interactions, the other for bonded.  Both settings are stored in the
restart file.  For the bond setting, the value in the file will
overwrite the current value (at the time the read_restart command is
issued) and warn if the two values are not the same and the current
value is not the default.  For the pair setting, the value in the file
will not overwrite the current value (so that you can override the
previous run's value), but a warning is issued if the two values are
not the same and the current value is not the default.

Note that some force field styles (pair, bond, angle, etc) do not
store their coefficient info in restart files.  Typically these are
many-body or tabulated potentials which read their parameters from
separate files.  In these cases you will need to re-specify the
:doc:`pair_coeff <pair_coeff>`, :doc:`bond_coeff <bond_coeff>`, etc
commands in your restart input script.  The doc pages for individual
force field styles mention if this is the case.  This is also true of
:doc:`pair_style hybrid <pair_hybrid>` (bond hybrid, angle hybrid, etc)
commands; they do not store coefficient info.

As indicated in the above list, the :doc:`fixes <fix>` used for a
simulation are not stored in the restart file.  This means the new
input script should specify all fixes it will use.  However, note that
some fixes store an internal "state" which is written to the restart
file.  This allows the fix to continue on with its calculations in a
restarted simulation.  To re-enable such a fix, the fix command in the
new input script must be of the same style and use the same fix-ID as
was used in the input script that wrote the restart file.

If a match is found, LAMMPS prints a message indicating that the fix
is being re-enabled.  If no match is found before the first run or
minimization is performed by the new script, the "state" information
for the saved fix is discarded.  At the time the discard occurs,
LAMMPS will also print a list of fixes for which the information is
being discarded.  See the doc pages for individual fixes for info on
which ones can be restarted in this manner.  Note that fixes which are
created internally by other LAMMPS commands (computes, fixes, etc)
will have style names which are all-capitalized, and IDs which are
generated internally.

Likewise, the :doc:`computes <fix>` used for a simulation are not stored
in the restart file.  This means the new input script should specify
all computes it will use.  However, some computes create a fix
internally to store "state" information that persists from timestep to
timestep.  An example is the :doc:`compute msd <compute_msd>` command
which uses a fix to store a reference coordinate for each atom, so
that a displacement can be calculated at any later time.  If the
compute command in the new input script uses the same compute-ID and
group-ID as was used in the input script that wrote the restart file,
then it will create the same fix in the restarted run.  This means the
re-created fix will be re-enabled with the stored state information as
described in the previous paragraph, so that the compute can continue
its calculations in a consistent manner.

.. note::

   There are a handful of commands which can be used before or between
   runs which may require a system initialization.  Examples include the
   "balance", "displace_atoms", "delete_atoms", "set" (some options),
   and "velocity" (some options) commands.  This is because they can
   migrate atoms to new processors.  Thus they will also discard unused
   "state" information from fixes.  You will know the discard has
   occurred because a list of discarded fixes will be printed to the
   screen and log file, as explained above.  This means that if you wish
   to retain that info in a restarted run, you must re-specify the
   relevant fixes and computes (which create fixes) before those
   commands are used.

Some pair styles, like the :doc:`granular pair styles <pair_gran>`, also
use a fix to store "state" information that persists from timestep to
timestep.  In the case of granular potentials, it is contact
information between pairs of touching particles.  This info will also
be re-enabled in the restart script, assuming you re-use the same
granular pair style.

LAMMPS allows bond interactions (angle, etc) to be turned off or
deleted in various ways, which can affect how their info is stored in
a restart file.

If bonds (angles, etc) have been turned off by the :doc:`fix shake
<fix_shake>` or :doc:`delete_bonds <delete_bonds>` command, their info
will be written to a restart file as if they are turned on.  This means
they will need to be turned off again in a new run after the restart
file is read.

Bonds that are broken (e.g. by a bond-breaking potential) are written
to the restart file as broken bonds with a type of 0.  Thus these
bonds will still be broken when the restart file is read.

Bonds that have been broken by the :doc:`fix bond/break
<fix_bond_break>` command have disappeared from the system.  No
information about these bonds is written to the restart file.

----------

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`read_dump <read_dump>`,
:doc:`write_restart <write_restart>`, :doc:`restart <restart>`

Default
"""""""

none
