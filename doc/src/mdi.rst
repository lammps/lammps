.. index:: mdi

mdi command
==================

Syntax
""""""

.. parsed-literal::

   mdi option args

* option = *engine* or *plugin* or *connect* or *exit*

  .. parsed-literal::

     *engine* args = zero or more keyword/args pairs
       keywords = *elements*
         *elements* args = N_1 N_2 ... N_ntypes
           N_1,N_2,...N_ntypes = atomic number for each of ntypes LAMMPS atom types
     *plugin* args = name keyword value keyword value ...
       name = name of plugin library, e.g. lammps means a liblammps.so library will be loaded
       keyword/value pairs in any order, some are required, some are optional
       keywords = *mdi* or *infile* or *extra* or *command*
         *mdi* value = args passed to MDI for driver to operate with plugins (required)
         *infile* value = filename the engine will read at start-up (optional)
         *extra* value = aditional command-line args to pass to engine library when loaded (optional)
         *command* value = a LAMMPS input script command to execute (required)
     *connect* args = none
     *exit* args = none

Examples
""""""""

.. code-block:: LAMMPS

   mdi engine
   mdi engine elements 13 29
   mdi plugin lammps mdi "-role ENGINE -name lammps -method LINK" &
              infile in.aimd.engine extra "-log log.aimd.engine.plugin" &
              command "run 5"
   mdi connect
   mdi exit

Description
"""""""""""

This command implements operations within LAMMPS to use the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`
for coupling to other codes in a client/server protocol.

See the Howto MDI doc page for a discussion of all the different ways
2 or more codes can interact via MDI.

The examples/mdi directory has examples which use LAMMPS in 4
different modes: as a driver using an engine as either a stand-alone
code or as a plugin, and as an engine operating as either a
stand-alone code or as a plugin.  The README file in that directory
shows how to launch and couple codes for all the 4 usage modes, and so
they communicate via the MDI library using either MPI or sockets.

The scripts in that directory illustrate the use of all the options
for this command.

The *engine* option enables LAMMPS to act as an MDI engine (server),
responding to requests from an MDI driver (client) code.

The *plugin* option enables LAMMPS to act as an MDI driver (client),
and load the MDI engine (server) code as a library plugin.  In this
case the MDI engine is a library plugin.  An MDI engine can also be a
stand-alone code, launched separately from LAMMPS, in which case the
mdi plugin command is not used.

The *connect* and *exit* options are only used when LAMMPS is acting
as an MDI driver.  As explained below, these options are normally not
needed, except for a specific kind of use case.

----------

The *mdi engine* command is used to make LAMMPS operate as an MDI
engine.  It is typically used in an input script after LAMMPS has
setup the system it is going to model consistent with what the driver
code expects.  Depending on when the driver code tells the LAMMPS
engine to exit, other commands can be executed after this command, but
typically it is used at the end of a LAMMPS input script.

To act as an MDI engine operating as an MD code (or surrogate QM
code), this is the list of standard MDI commands issued by a driver
code which LAMMPS currently recognizes.  Using standard commands
defined by the MDI library means that a driver code can work
interchangeably with LAMMPS or other MD codes or with QM codes which
support the MDI standard.  See more details about these commands in
the `MDI library documentation
<https://molssi-mdi.github.io/MDI_Library/html/mdi_standard.html>`_

These commands are valid at the @DEFAULT node defined by MDI.
Commands that start with ">" mean the driver is sending information to
LAMMPS.  Commands that start with "<" are requests by the driver for
LAMMPS to send it information.  Commands that start with an alphabetic
letter perform actions.  Commands that start with "@" are MDI "node"
commands, which are described further below.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Command name
     - Action
   * - >CELL or <CELL
     - Send/request 3 simulation box edge vectors (9 values)
   * - >CELL_DISPL or <CELL_DISPL
     - Send/request displacement of the simulation box from the origin (3 values)
   * - >CHARGES or <CHARGES
     - Send/request charge on each atom (N values)
   * - >COORDS or <COORDS
     - Send/request coordinates of each atom (3N values)
   * - >ELEMENTS
     - Send elements (atomic numbers) for each atom (N values)
   * - <ENERGY
     - Request total energy (potential + kinetic) of the system (1 value)
   * - >FORCES or <FORCES
     - Send/request forces on each atom (3N values)
   * - >+FORCES
     - Send forces to add to each atom (3N values)
   * - <LABELS
     - Request string label of each atom (N values)
   * - <MASSES
     - Request mass of each atom (N values)
   * - MD
     - Perform an MD simulation for N timesteps (most recent >NSTEPS value)
   * - OPTG
     - Perform an energy minimization to convergence (most recent >TOLERANCE values)
   * - >NATOMS or <NATOMS
     - Sends/request number of atoms in the system (1 value)
   * - >NSTEPS
     - Send number of timesteps for next MD dynamics run via MD command
   * - <PE
     - Request potential energy of the system (1 value)
   * - <STRESS
     - Request symmetric stress tensor (virial) of the system (9 values)
   * - >TOLERANCE
     - Send 4 tolerance parameters for next MD minimization via OPTG command
   * - >TYPES or <TYPES
     - Send/request the LAMMPS atom type for each atom (N values)
   * - >VELOCITIES or <VELOCITIES
     - Send/request the velocity of each atom (3N values)
   * - @INIT_MD or @INIT_OPTG
     - Driver tells LAMMPS to start single-step dynamics or minimization (see below)
   * - EXIT
     - Driver tells LAMMPS to exit engine mode

.. note::

   The <ENERGY, <FORCES, <PE, and <STRESS commands trigger LAMMPS to
   compute atomic interactions for the current configuration of atoms
   and size/shape of the simulation box.  I.e. LAMMPS invokes its
   pair, bond, angle, ..., kspace styles.  If the driver is updating
   the atom coordinates and/or box incrementally (as in an MD
   simulation which the driver is managing), then the LAMMPS engine
   will do the same, and only occasionally trigger neighbor list
   builds.  If the change in atom positions is large (since the
   previous >COORDS command), then LAMMPS will do a more expensive
   operation to migrate atoms to new processors as needed and
   re-neighbor.  If the >NATOMS or >TYPES or >ELEMENTS commands have
   been sent (since the previous >COORDS command), then LAMMPS assumes
   the system is new and re-initializes an entirely new simulation.

.. note::

   The >TYPES or >ELEMENTS commands are how the MDI driver tells the
   LAMMPS engine which LAMMPS atom type to assign to each atom.  If
   both the MDI driver and the LAMMPS engine are initialized so that
   atom type values are consistent in both codes, then the >TYPES
   command can be used.  If not, the optional *elements* keyword can
   be used to specify what element each LAMMPS atom type corresponds
   to.  This is specified by the atomic number of the element, e.g. 13
   for Al.  An atomic number must be specified for each of the ntypes
   LAMMPS atom types.  Ntypes is typically specified via the
   create_box command or in the data file read by the read_data
   command.  In this has been done, the MDI driver can send an
   >ELEMENTS command to the LAMMPS driver with the atomic number of
   each atom.

The MD and OPTG commands perform an entire MD simulation or energy
minimization (to convergence) with no communication from the driver
until the simulation is complete.  By contrast, the @INIT_MD and
@INIT_OPTG commands allow the driver to communicate with the engine at
each timestep of a dynamics run or iteration of a minimization; see
more info below.

The MD command performs a simulation using the most recent >NSTEPS
value.  The OPTG command performs a minimization using the 4
convergence parameters from the most recent >TOLERANCE command.  The 4
parameters sent are those used by the :doc:`minimize <minimize>`
command in LAMMPS: etol, ftol, maxiter, and maxeval.

The mdi engine command also implements the following custom MDI
commands which are LAMMPS-specific.  These commands are also valid at
the @DEFAULT node defined by MDI:

   * - Command name
     - Action
   * - >NBYTES
     - Send # of datums in a subsequent command (1 value)
   * - >COMMAND
     - Send a LAMMPS input script command as a string (Nbytes in length)
   * - >COMMANDS
     - Send multiple LAMMPS input script commands as a newline-separated string (Nbytes in length)
   * - >INFILE
     - Send filename of an input script to execute (filename Nbytes in length)
   * - <KE
     - Request kinetic energy of the system (1 value)

Note that other custom commands can easily be added if these are not
sufficient to support what a user-written driver code needs.  Code to
support new commands can be added to the MDI package within LAMMPS,
specifically to the src/MDI/mdi_engine.cpp file.

MDI also defines a standard mechanism for the driver to request that
an MD engine (LAMMPS) perform a dynamics simulation one step at a time
or an energy minimization one iteration at a time.  This is so that
the driver can (optionally) communicate with LAMMPS at intermediate
points of the timestep or iteration by issuing MDI node commands which
start with "@".

To tell LAMMPS to run dynamics in single-step mode, the driver sends
as @INIT_MD command followed by the these commands.  The driver
can interact with LAMMPS at 3 node locations within each
timestep: @COORDS, @FORCES, @ENDSTEP:

   * - Command name
     - Action
   * - @COORDS
     - Proceed to next @COORDS node = post-integrate location in LAMMPS timestep
   * - @FORCES
     - Proceed to next @FORCES node = post-force location in LAMMPS timestep
   * - @ENDSTEP
     - Proceed to next @ENDSTEP node = end-of-step location in LAMMPS timestep
   * - @DEFAULT
     - Exit MD simulation, return to @DEFAULT node
   * - EXIT
     - Driver tells LAMMPS to exit the MD simulation and engine mode

To tell LAMMPS to run an energy minimization in single-iteration mode.
The driver can interact with LAMMPS at 2 node locations within each
iteration of the minimizer: @COORDS, @FORCES:

   * - Command name
     - Action
   * - @COORDS
     - Proceed to next @COORDS node = min-pre-force location in LAMMPS min iteration
   * - @FORCES
     - Proceed to next @FORCES node = min-post-force location in LAMMPS min iteration
   * - @DEFAULT
     - Exit minimization, return to @DEFAULT node
   * - EXIT
     - Driver tells LAMMPS to exit the minimization and engine mode

While LAMMPS is at its @COORDS node, the following standard MDI
commands are supported, as documented above: >COORDS or <COORDS,
@COORDS, @FORCES, @ENDSTEP, @DEFAULT, EXIT.

While LAMMPS is at its @FORCES node, the following standard MDI
commands are supported, as documented above: <COORDS, <ENERGY, >FORCES
or >+FORCES or <FORCES, <KE, <PE, <STRESS, @COORDS, @FORCES, @ENDSTEP,
@DEFAULT, EXIT.

While LAMMPS is at its @ENDSTEP node, the following standard MDI
commands are supported, as documented above: <ENERGY, <FORCES, <KE,
<PE, <STRESS, @COORDS, @FORCES, @ENDSTEP, @DEFAULT, EXIT.

----------

The *mdi plugin* command is used to make LAMMPS operate as an MDI
driver which loads an MDI engine as a plugin library.  It is typically
used in an input script after LAMMPS has setup the system it is going
to model consistent with the engine code.

The *name* argument specifies which plugin library to load.  A name
like "lammps" is converted to a filename liblammps.so.  The path for
where this file is located is specified by the -plugin_path switch
within the -mdi command-line switch, which is specified when LAMMPS is
launched.  See the examples/mdi/README files for examples of how this
is done.

The *mdi* keyword is required and is used as the -mdi argument passed
to the library when it is launched.  The -role and -method settings
are required.  The -name setting can be anything you choose.  MDI
drivers and engines can query their names to verify they are values
they expect.

The *infile* keyword is optional.  It sets the name of an input script
which the engine will open and process.  MDI will pass it as a
command-line argument to the library when it is launched.  The file
typically contains settings that an MD or QM code will use for its
calculations.

The *extra* keyword is optional.  It contains additional command-line
arguments which MDI will pass to the library when it is launched.

The *command* keyword is required.  It specifies a LAMMPS input script
command (as a single argument in quotes if it is multiple words).
Once the plugin library is launched, LAMMPS will execute this command.
Other previously-defined commands in the input script, such as the
:doc:`fix mdi/qm <fix_mdi_qm>` command, should perform MDI
communication with the engine, while the specified *command* executes.
Note that if *command* is an :doc:`include <include>` command, then it
could specify a filename with multiple LAMMPS commands.

.. note::

   When the *command* is complete, LAMMPS will send an MDI EXIT
   command to the plugin engine and the plugin will be removed.  The
   "mdi plugin" command will then exit and the next command (if any)
   in the LAMMPS input script will be processed.  A subsequent "mdi
   plugin" command could then load the same or a different MDI
   plugin if desired.

----------

The *mdi connect* and *mdi exit* commands are only used when LAMMPS is
operating as an MDI driver.  And when other LAMMPS command(s) which
send MDI commands and associated data to/from the MDI engine are not
able to initiate and terminate the connection to the engine code.

The only current MDI driver command in LAMMPS is the :doc:`fix mdi/qm
<fix_mdi_qm>` command.  If it is only used once in an input script
then it can initiate and terminate the connection.  But if it is being
issued multiple times, e.g. in a loop that issues a :doc:`clear
<clear>` command, then it cannot initiate or terminate the connection
multiple times.  Instead, the *mdi connect* and *mdi exit* commands
should be used outside the loop to initiate or terminate the connection.

See the examples/mdi/in.series.driver script for an example of how
this is done.  The LOOP in that script is reading a series of data
file configurations and passing them to an MDI engine (e.g. quantum
code) for energy and force evaluation.  A *clear* command inside the
loop wipes out the current system so a new one can be defined.  This
operation also destroys all fixes.  So the :doc:`fix mdi/qm
<fix_mdi_qm>` command is issued once per loop iteration.  Note that it
includes a "connect no" option which disables the initiate/terminate
logic within that fix.


Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To use LAMMPS in conjunction with other MDI-enabled atomistic codes,
the :doc:`units <units>` command should be used to specify *real* or
*metal* units.  This will ensure the correct unit conversions between
LAMMPS and MDI units, which the other codes will also perform in their
preferred units.

LAMMPS can also be used as an MDI engine in other unit choices it
supports, e.g. *lj*, but then no unit conversion is performed.

Related commands
""""""""""""""""

:doc:`fix mdi/qm <fix_mdi_qm>`

Default
"""""""

None
