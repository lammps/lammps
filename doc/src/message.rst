.. index:: message

message command
===============

Syntax
""""""


.. parsed-literal::

   message which protocol mode arg

* which = *client* or *server*
* protocol = *md* or *mc*
* mode = *file* or *zmq* or *mpi/one* or *mpi/two*
  
  .. parsed-literal::
  
       *file* arg = filename
         filename = file used for message exchanges
       *zmq* arg = socket-ID
         socket-ID for client = localhost:5555, see description below
         socket-ID for server = \*:5555, see description below
       *mpi/one* arg = none
       *mpi/two* arg = filename
         filename = file used to establish communication between 2 MPI jobs



Examples
""""""""


.. parsed-literal::

   message client md file tmp.couple
   message server md file tmp.couple

   message client md zmq localhost:5555
   message server md zmq \*:5555

   message client md mpi/one
   message server md mpi/one

   message client md mpi/two tmp.couple
   message server md mpi/two tmp.couple

Description
"""""""""""

Establish a messaging protocol between LAMMPS and another code for the
purpose of client/server coupling.

The :doc:`Howto client/server <Howto_client_server>` doc page gives an
overview of client/server coupling of LAMMPS with another code where
one code is the "client" and sends request messages to a "server"
code.  The server responds to each request with a reply message.  This
enables the two codes to work in tandem to perform a simulation.


----------


The *which* argument defines LAMMPS to be the client or the server.


----------


The *protocol* argument defines the format and content of messages
that will be exchanged between the two codes.  The current options
are:

* md = run dynamics with another code
* mc = perform Monte Carlo moves with another code

For protocol *md*\ , LAMMPS can be either a client or server.  See the
:doc:`server md <server_md>` doc page for details on the protocol.

For protocol *mc*\ , LAMMPS can be the server.  See the :doc:`server mc <server_mc>` doc page for details on the protocol.


----------


The *mode* argument specifies how messages are exchanged between the
client and server codes.  Both codes must use the same mode and use
consistent parameters.

For mode *file*\ , the 2 codes communicate via binary files.  They must
use the same filename, which is actually a file prefix.  Several files
with that prefix will be created and deleted as a simulation runs.
The filename can include a path.  Both codes must be able to access
the path/file in a common filesystem.

For mode *zmq*\ , the 2 codes communicate via a socket on the server
code's machine.  Support for socket messaging is provided by the
open-source `ZeroMQ library <http://zeromq.org>`_, which must be
installed on your system.  The client specifies an IP address (IPv4
format) or the DNS name of the machine the server code is running on,
followed by a 4-digit port ID for the socket, separated by a colon.
E.g.


.. parsed-literal::

   localhost:5555        # client and server running on same machine
   192.168.1.1:5555      # server is 192.168.1.1
   deptbox.uni.edu:5555  # server is deptbox.uni.edu

The server specifies "\*:5555" where "\*" represents all available
interfaces on the server's machine, and the port ID must match
what the client specifies.

.. note::

   What are allowed port IDs?

.. note::

   Additional explanation is needed here about how to use the *zmq*
   mode on a parallel machine, e.g. a cluster with many nodes.

For mode *mpi/one*\ , the 2 codes communicate via MPI and are launched
by the same mpirun command, e.g. with this syntax for OpenMPI:


.. parsed-literal::

   mpirun -np 2 lmp_mpi -mpicolor 0 -in in.client -log log.client : -np 4 othercode args  # LAMMPS is client
   mpirun -np 2 othercode args : -np 4 lmp_mpi -mpicolor 1 -in in.server  # LAMMPS is server

Note the use of the "-mpicolor color" command-line argument with
LAMMPS.  See the :doc:`command-line args <Run_options>` doc page for
further explanation.

For mode *mpi/two*\ , the 2 codes communicate via MPI, but are launched
be 2 separate mpirun commands.  The specified *filename* argument is a
file the 2 MPI processes will use to exchange info so that an MPI
inter-communicator can be established to enable the 2 codes to send
MPI messages to each other.  Both codes must be able to access the
path/file in a common filesystem.


----------


Normally, the message command should be used at the top of a LAMMPS
input script.  It performs an initial handshake with the other code to
setup messaging and to verify that both codes are using the same
message protocol and mode.  Assuming both codes are launched at
(nearly) the same time, the other code should perform the same kind of
initialization.

If LAMMPS is the client code, it will begin sending messages when a
LAMMPS client command begins its operation.  E.g. for the :doc:`fix client/md <fix_client_md>` command, it is when a :doc:`run <run>`
command is executed.

If LAMMPS is the server code, it will begin receiving messages when
the :doc:`server <server>` command is invoked.

A fix client command will terminate its messaging with the server when
LAMMPS ends, or the fix is deleted via the :doc:`unfix <unfix>` command.
The server command will terminate its messaging with the client when the
client signals it.  Then the remainder of the LAMMPS input script will
be processed.

If both codes do something similar, this means a new round of
client/server messaging can be initiated after termination by re-using
a 2nd message command in your LAMMPS input script, followed by a new
fix client or server command.


----------


Restrictions
""""""""""""


This command is part of the MESSAGE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`server <server>`, :doc:`fix client/md <fix_client_md>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
