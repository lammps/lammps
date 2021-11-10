Using LAMMPS with the MDI library for code coupling
===================================================

.. note::

  This Howto page will eventually replace the
  :doc:`Howto client/server <Howto_client_server>` doc page.

Client/server coupling of two codes is where one code is the "client"
and sends request messages (data) to a "server" code.  The server
responds to each request with a reply message.  This enables the two
codes to work in tandem to perform a simulation.  LAMMPS can act as
either a client or server code; it does this by using the `MolSSI
Driver Interface (MDI) library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_,
developed by the `Molecular Sciences Software Institute (MolSSI)
<https://molssi.org>`_.

Alternate methods for code coupling with LAMMPS are described on the
:doc:`Howto couple <Howto_couple>` doc page.

Some advantages of client/server coupling are that the two codes can run
as stand-alone executables; they need not be linked together.  Thus
neither code needs to have a library interface.  This also makes it easy
to run the two codes on different numbers of processors.  If a message
protocol (format and content) is defined for a particular kind of
simulation, then in principle any code which implements the client-side
protocol can be used in tandem with any code which implements the
server-side protocol.  Neither code needs to know what specific other
code it is working with.

In MDI nomenclature, a client code is the "driver", and a server code is
an "engine".  One driver code can communicate with one or more instances
of one or more engine codes.  Driver and engine codes can be written in
any language: C, C++, Fortran, Python, etc.

In addition to allowing driver and engine(s) running to run as
stand-alone executables, MDI also enables a server code to be a
"plugin" to the client code.  In this scenario, server code(s) are
compiled as shared libraries, and one (or more) instances of the
server are instantiated by the driver code.  If the driver code runs
in parallel, it can split its MPI communicator into multiple
sub-communicators, and launch each plugin engine instance on a
sub-communicator.  Driver processors in that sub-communicator exchange
messages with that engine instance, and can also send MPI messages to
other processors in the driver.  The driver code can also destroy
engine instances and re-instantiate them.

The way that a driver communicates with an engine is by making
MDI_Send() and MDI_Recv() calls, which are conceptually similar to
MPI_Send() and MPI_Recv() calls.  Each send or receive has a string
which identifies the command name, and optionally some data, which can
be a single value or vector of values of any data type.  Inside the
MDI library, data is exchanged between the driver and engine via MPI
calls or sockets.  This a run-time choice by the user.

-------------

As an example, LAMMPS and the ``pw.x`` command from Quantum Espresso (a
suite of quantum DFT codes), can work together via the MDI library to
perform an ab initio MD (AIMD) simulation, where LAMMPS runs an MD
simulation and sends a message each timestep to ``pw.x`` asking it to
compute quantum forces on the current configuration of atoms.  Here is
how the 2 codes are launched to communicate by MPI:

.. code-block:: bash

   % mpirun -np 2 lmp_mpi -mdi "-role DRIVER -name d -method MPI" \
     -in in.aimd : -np 16 pw.x -in qe.in -mdi "-role ENGINE -name e -method MPI"

In this case LAMMPS runs on 2 processors (MPI tasks), ``pw.x`` runs on 16
processors.

Here is how the 2 codes are launched to communicate by sockets:

.. code-block:: bash

   % mpirun -np 2 lmp_mpi -mdi "-role DRIVER -name d -method TCP -port 8021" -in in.aimd
   % mpirun -np 16 pw.x -in qe.in -mdi "-role ENGINE -name e -method TCP -port 8021 -hostname localhost"

These commands could be issued in different windows on a desktop
machine.  Or in the same window, if the first command is ended with
"&" so as to run in the background.  If "localhost" is replaced by an
IP address, ``pw.x`` could be run on another machine on the same network, or
even on another machine across the country.

After both codes initialize themselves to model the same system, this is
what occurs each timestep:

* LAMMPS send a ">COORDS" message to ``pw.x`` with a 3*N vector of current atom coords
* ``pw.x`` receives the message/coords and computes quantum forces on all the atoms
* LAMMPS send a "<FORCES" message to ``pw.x`` and waits for the result
* ``pw.x`` receives the message (after its computation finishes) and sends a 3*N vector of forces
* LAMMPS receives the forces and time integrates to complete a single timestep

-------------

Examples scripts for using LAMMPS as an MDI engine are in the
examples/mdi directory.  See the README file in that directory for
instructions on how to run the examples.

.. note::

  Work is underway to add commands that allow LAMMPS to be used as an
  MDI driver, e.g. for the AIMD example discussed above.  Example
  scripts for this usage mode will be added the same directory when
  available.

If LAMMPS is used as a stand-alone engine it should set up the system
it will be modeling in its input script, then invoke the
:doc:`mdi/engine <mdi_engine>` command.  This will put LAMMPS into
"engine mode" where it waits for messages and data from the driver.
When the driver sends an "EXIT" command, LAMMPS will exit engine mode
and the input script will continue.

If LAMMPS is used as a plugin engine it operates the same way, except
that the driver will pass LAMMPS an input script to initialize itself.
Upon receiving the "EXIT" command, LAMMPS will exit engine mode and the
input script will continue.  After finishing execution of the input
script, the instance of LAMMPS will be destroyed.

LAMMPS supports the full set of MD-appropriate engine commands defined
by the MDI library.  See the :doc:`mdi/engine <mdi_engine>` page for
a list of these.

If those commands are not sufficient for a user-developed driver to use
LAMMPS as an engine, then new commands can be easily added.  See these
two files which implement the definition of MDI commands and the logic
for responding to them:

* src/MDI/mdi_engine.cpp
* src/MDI/fix_mdi_engine.cpp
