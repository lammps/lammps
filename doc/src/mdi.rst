.. index:: mdi

mdi command
==================

Syntax
""""""

.. parsed-literal::

   mdi mode args

* mode = *engine* or *plugin*

  .. parsed-literal::

       *engine* args = none
       *plugin* args = keyword value keyword value ...
         keywords = *mdi* or *infile* or *extra* or *command*
           *mdi* value = args passed to MDI for driver to operate with plugins
           *infile* value = filename the engine will read at start-up
           *extra* value = aditional command-line args to pass to engine library when loaded
           *command value = a LAMMPS input script command to execute

Examples
""""""""

.. code-block:: LAMMPS

   mdi engine
   mdi plugin lammps mdi "-role ENGINE -name lammps -method LINK" &
              infile in.aimd.engine extra "-log log.aimd.engine.plugin" &
              command "run 5"

Description
"""""""""""

This command implements two high-level operations within LAMMPS to use
the `MDI Library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>` for
coupling to another code.

The *engine* mode enables LAMMPS to act as a server, responding to
requests from another client code to effectively couple the two codes
together in client/server mode.  In MDI lingo this means LAMMPS
operates as an MDI engine.

The *plugin* mode enables LAMMPS to act as a client, and load the
server code as a library plugin.  In MDI lingo this means LAMMPS
operates as an MDI driver, and the server is an MDI engine.  In this
case the MDI engine is a library plugin.  It can also be a stand-alone
code, lauched separately from LAMMPS.  

See the Howto MDI doc page for a discussion of all the different ways
2 or more codes can interact via MDI.  The examples/mdi/README file
has example use cases and launch commands for using LAMMPS as both a
driver and engine in all these different ways.

The examples/mdi directory contains input scripts for LAMMPS acting as
an MDI engine to operate as a surrogate quantum mechanics (QM) code
for running ab initio MD (AIMD) or within a materials modeling
workflow to perform various MD tasks.  It likewise has example scripts
for LAMMPS acting as an MDI driver for AIMD simulations, which could
use real QM codes which support MDI.  The examples/mdi/README file
explains how to launch both driver and engine codes so that they
communicate using the MDI library via either MPI or sockets.

----------

The *mdi engine* command should typically be used in an input script
after LAMMPS has setup the system it is going to model in
collaboration with the driver code.  Depending on how the driver code
tells the LAMMPS engine to exit, other commands can be executed after
this command, but typically it should it is used at the end of a
LAMMPS input script.



To act as an MDI engine operating as an MD code (or surrogate QM
code), this is the list of standard MDI commands issued by a driver
code which LAMMPS currently recognizes.  Using standard commands
defined by the MDI library means that a driver code can work
interchangeably with LAMMPS or other MD codes which support the MDI
standard.  See more details about these commands in the `MDI library
documentation
<https://molssi-mdi.github.io/MDI_Library/html/mdi_standard.html>`_

These commands are valid at the @DEFAULT node defined by MDI.
Commands that start with ">" mean the driver is sending information to
the engine (LAMMMPS).  Commands that start with "<" are requests by
the driver for LAMMPS to send it information.  Commands that start
with a letter perform actions.  Commands that start with "@" are MDI
"node" commands, which are described further below.

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
       Perform an MD simulation for N timestpes (most recent >NSTEPS value)
   * - OPTG
       Perform an energy minimization to convergence (most recent >TOLERANCE values)
   * - >NATOMS or <NATOMS
     - Sends/request number of atoms in the system (1 value)
   * - >NSTEPS
     - Send number of timesteps for next MD dynamics run via MD command
   * - <PE
     - Request potential energy of the system (1 value)
   * - <STRESS
     - Request stress tensor (virial) of the system (6 values)
   * - >TOLERANCE
     - Send 4 tolerance parameters for next MD minimization via OPTG command
   * - >TYPES or <TYPES
     - Send/request the numeric type of each atom (N values)
   * - >VELOCITIES or <VELOCITIES
     - Send/request the velocity of each atom (3N values)
   * - @INIT_MD or @INIT_OPTG
     - Driver tells LAMMPS to start single-step dynamics or minimization (see below)
   * - EXIT
     - Driver tells LAMMPS to exit engine mode


.. note::

   The <ENERGY, <FORCES, <PE, and <STRESS commands trigger LAMMPS to
   invoke its pair, bond, angle, ..., kspace styles, i.e. to calculate
   the force field for the current configuration of atoms and
   size/shape of the simulation box.  If the driver is updating the
   atom coordinates and/or box incrementally (as in an MD simulation
   which the driver is managing), then the LAMMPS engine will do the
   same, and only occasionally trigger neighbor list builds.  If the
   change in atom positions is large (since the previous >COORDS
   command), then LAMMPS will do a more expensive operation to migrate
   atoms to new processors as needed and re-neighbor.  If the >NATOMS
   or >TYPES commands have been sent (since the previous >COORDS
   command), then LAMMPS assumes the system is new and re-initializes
   an entirely new simulation.

The MD and OPTG commands perform an entire MD simulation or energy
minimization (to convergence) with no communication with the driver
until the simulation is complete.  By contrast, the @INIT_MD and
@INIT_OPTG commands allow the driver to communicate with the engine at
each timestep of a dynamics run or iteration of a minimization; see
more info below.

The MD command performa a simulation using the most recent >NSTEPS
value.  The OPTG command performs a minimization using the 4
convergence paremeters from the most recent >TOLERANCE command.  The 4
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
     - Send filename of an input script to execute (Nbytes in length)
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

The *mdi plugin* command should also typically be used at the end of
an input script after LAMMPS has setup the system it is going to model
in collaboration with the engine code.



Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To use LAMMPS as an MDI engine in conjunction with other MDI-enabled
atomistic codes, the :doc:`units <units>` command should be used to
specify *real* or *metal* units.  This will ensure the correct unit
conversions between LAMMPS and MDI units, which the other codes will
also perform in their preferred units.

LAMMPS can also be used as an MDI engine in other unit choices it
supports, e.g. *lj*, but then no unit conversion is performed.

Related commands
""""""""""""""""

:doc:`fix mdi/aimd <fix_mdi_aimd>`

Default
"""""""

None
