.. index:: server mc

server mc command
=================

Syntax
""""""


.. parsed-literal::

   server mc

mc = the protocol argument to the :doc:`server <server>` command

Examples
""""""""


.. parsed-literal::

   server mc

Description
"""""""""""

This command starts LAMMPS running in "server" mode, where it will
expect messages from a separate "client" code that match the *mc*
protocol for format and content explained below.  For each message
LAMMPS receives it will send a message back to the client.

The :doc:`Howto client/server <Howto_client_server>` doc page gives an
overview of client/server coupling of LAMMPS with another code where
one code is the "client" and sends request messages to a "server"
code.  The server responds to each request with a reply message.  This
enables the two codes to work in tandem to perform a simulation.

When this command is invoked, LAMMPS will run in server mode in an
endless loop, waiting for messages from the client code.  The client
signals when it is done sending messages to LAMMPS, at which point the
loop will exit, and the remainder of the LAMMPS script will be
processed.

The :doc:`server <server>` doc page gives other options for using LAMMPS
See an example of how this command is used in
examples/COUPLE/lammps\_mc/in.server.


----------


When using this command, LAMMPS (as the server code) receives
instructions from a Monte Carlo (MC) driver to displace random atoms,
compute the energy before and after displacement, and run dynamics to
equilibrate the system.

The MC driver performs the random displacements on random atoms,
accepts or rejects the move in an MC sense, and orchestrates the MD
runs.

The format and content of the exchanged messages are explained here in
a conceptual sense.  Python-style pseudo code for the library calls to
the CSlib is shown, which performs the actual message exchange between
the two codes.  See the `CSlib website <http://cslib.sandia.gov>`_ doc
pages for more details on the actual library syntax.  The "cs" object
in this pseudo code is a pointer to an instance of the CSlib.

See the src/MESSAGE/server\_mc.cpp file for details on how LAMMPS uses
these messages.  See the examples/COUPLE/lammps\_mc/mc.cpp file for an
example of how an MC driver code can use these messages.

Define NATOMS=1, EINIT=2, DISPLACE=3, ACCEPT=4, RUN=5.

**Client sends one of these kinds of message**\ :


.. parsed-literal::

   cs->send(NATOMS,0)      # msgID = 1 with no fields

   cs->send(EINIT,0)       # msgID = 2 with no fields

   cs->send(DISPLACE,2)    # msgID = 3 with 2 fields
   cs->pack_int(1,ID)        # 1st field = ID of atom to displace
   cs->pack(2,3,xnew)      # 2nd field = new xyz coords of displaced atom

   cs->send(ACCEPT,1)      # msgID = 4 with 1 field
   cs->pack_int(1,flag)    # 1st field = accept/reject flag

   cs->send(RUN,1)         # msgID = 5 with 1 field
   cs->pack_int(1,nsteps)  # 1st field = # of timesteps to run MD

**Server replies**\ :


.. parsed-literal::

   cs->send(NATOMS,1)      # msgID = 1 with 1 field
   cs->pack_int(1,natoms)  # 1st field = number of atoms

   cs->send(EINIT,2)         # msgID = 2 with 2 fields
   cs->pack_double(1,poteng) # 1st field = potential energy of system
   cs->pack(2,3\*natoms,x)    # 2nd field = 3N coords of Natoms

   cs->send(DISPLACE,1)      # msgID = 3 with 1 field
   cs->pack_double(1,poteng) # 1st field = new potential energy of system

   cs->send(ACCEPT,0)      # msgID = 4 with no fields

   cs->send(RUN,0)         # msgID = 5 with no fields


----------


Restrictions
""""""""""""


This command is part of the MESSAGE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

A script that uses this command must also use the
:doc:`message <message>` command to setup the messaging protocol with
the other client code.

Related commands
""""""""""""""""

:doc:`message <message>`

**Default:** none


