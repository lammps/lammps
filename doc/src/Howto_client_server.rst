Using LAMMPS in client/server mode
==================================

Client/server coupling of two codes is where one code is the "client"
and sends request messages to a "server" code.  The server responds to
each request with a reply message.  This enables the two codes to work
in tandem to perform a simulation.  LAMMPS can act as either a client
or server code.

Some advantages of client/server coupling are that the two codes run
as stand-alone executables; they are not linked together.  Thus
neither code needs to have a library interface.  This often makes it
easier to run the two codes on different numbers of processors.  If a
message protocol (format and content) is defined for a particular kind
of simulation, then in principle any code that implements the
client-side protocol can be used in tandem with any code that
implements the server-side protocol, without the two codes needing to
know anything more specific about each other.

A simple example of client/server coupling is where LAMMPS is the
client code performing MD timestepping.  Each timestep it sends a
message to a server quantum code containing current coords of all the
atoms.  The quantum code computes energy and forces based on the
coords.  It returns them as a message to LAMMPS, which completes the
timestep.

A more complex example is where LAMMPS is the client code and
processes a series of data files, sending each configuration to a
quantum code to compute energy and forces.  Or LAMMPS runs dynamics
with an atomistic force field, but pauses every N steps to ask the
quantum code to compute energy and forces.

Alternate methods for code coupling with LAMMPS are described on
the :doc:`Howto couple <Howto_couple>` doc page.

The protocol for using LAMMPS as a client is to use these 3 commands
in this order (other commands may come in between):

* :doc:`message client <message>`  # initiate client/server interaction
* :doc:`fix client/md <fix_client_md>`   # any client fix which makes specific requests to the server
* :doc:`message quit <message>`    # terminate client/server interaction

In between the two message commands, a client fix command and
:doc:`unfix <unfix>` command can be used multiple times.  Similarly,
this sequence of 3 commands can be repeated multiple times, assuming
the server program operates in a similar fashion, to initiate and
terminate client/server communication.

The protocol for using LAMMPS as a server is to use these 2 commands
in this order (other commands may come in between):

* :doc:`message server <message>`  # initiate client/server interaction
* :doc:`server md <server_md>`    # any server command which responds to specific requests from the client

This sequence of 2 commands can be repeated multiple times, assuming
the client program operates in a similar fashion, to initiate and
terminate client/server communication.

LAMMPS support for client/server coupling is in its :ref:`MESSAGE package <PKG-MESSAGE>` which implements several
commands that enable LAMMPS to act as a client or server, as discussed
below.  The MESSAGE package also wraps a client/server library called
CSlib which enables two codes to exchange messages in different ways,
either via files, sockets, or MPI.  The CSlib is provided with LAMMPS
in the lib/message dir.  The CSlib has its own
`website <http://cslib.sandia.gov>`_ with documentation and test
programs.

.. note::

   For client/server coupling to work between LAMMPS and another
   code, the other code also has to use the CSlib.  This can often be
   done without any modification to the other code by simply wrapping it
   with a Python script that exchanges CSlib messages with LAMMPS and
   prepares input for or processes output from the other code.  The other
   code also has to implement a matching protocol for the format and
   content of messages that LAMMPS exchanges with it.

These are the commands currently in the MESSAGE package for two
protocols, MD and MC (Monte Carlo).  New protocols can easily be
defined and added to this directory, where LAMMPS acts as either the
client or server.

* :doc:`message <message>`
* :doc:`fix client md <fix_client_md>` = LAMMPS is a client for running MD
* :doc:`server md <server_md>` = LAMMPS is a server for computing MD forces
* :doc:`server mc <server_mc>` = LAMMPS is a server for computing a Monte Carlo energy

The server doc files give details of the message protocols
for data that is exchanged between the client and server.

These example directories illustrate how to use LAMMPS as either a
client or server code:

* examples/message
* examples/COUPLE/README
* examples/COUPLE/lammps\_mc
* examples/COUPLE/lammps\_vasp

The examples/message dir couples a client instance of LAMMPS to a
server instance of LAMMPS.

The lammps\_mc dir shows how to couple LAMMPS as a server to a simple
Monte Carlo client code as the driver.

The lammps\_vasp dir shows how to couple LAMMPS as a client code
running MD timestepping to VASP acting as a server providing quantum
DFT forces, through a Python wrapper script on VASP.

Here is how to launch a client and server code together for any of the
4 modes of message exchange that the :doc:`message <message>` command
and the CSlib support.  Here LAMMPS is used as both the client and
server code.  Another code could be substituted for either.

The examples below show launching both codes from the same window (or
batch script), using the "&" character to launch the first code in the
background.  For all modes except *mpi/one*\ , you could also launch the
codes in separate windows on your desktop machine.  It does not
matter whether you launch the client or server first.

In these examples either code can be run on one or more processors.
If running in a non-MPI mode (file or zmq) you can launch a code on a
single processor without using mpirun.

IMPORTANT: If you run in mpi/two mode, you must launch both codes via
mpirun, even if one or both of them runs on a single processor.  This
is so that MPI can figure out how to connect both MPI processes
together to exchange MPI messages between them.

For message exchange in *file*\ , *zmq*\ , or *mpi/two* modes:


.. parsed-literal::

   % mpirun -np 1 lmp_mpi -log log.client < in.client &
   % mpirun -np 2 lmp_mpi -log log.server < in.server

   % mpirun -np 4 lmp_mpi -log log.client < in.client &
   % mpirun -np 1 lmp_mpi -log log.server < in.server

   % mpirun -np 2 lmp_mpi -log log.client < in.client &
   % mpirun -np 4 lmp_mpi -log log.server < in.server

For message exchange in *mpi/one* mode:

Launch both codes in a single mpirun command:


.. parsed-literal::

   mpirun -np 2 lmp_mpi -mpicolor 0 -in in.message.client -log log.client : -np 4 lmp_mpi -mpicolor 1 -in in.message.server -log log.server

The two -np values determine how many procs the client and the server
run on.

A LAMMPS executable run in this manner must use the -mpicolor color
command-line option as their its option, where color is an integer
label that will be used to distinguish one executable from another in
the multiple executables that the mpirun command launches.  In this
example the client was colored with a 0, and the server with a 1.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
