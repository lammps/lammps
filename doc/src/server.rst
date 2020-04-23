.. index:: server

server command
==============

Syntax
""""""

.. parsed-literal::

   server protocol

* protocol = *md* or *mc*

Examples
""""""""

.. code-block:: LAMMPS

   server md

Description
"""""""""""

This command starts LAMMPS running in "server" mode, where it receives
messages from a separate "client" code and responds by sending a reply
message back to the client.  The specified *protocol* determines the
format and content of messages LAMMPS expects to receive and how it
responds.

The :doc:`Howto client/server <Howto_client_server>` doc page gives an
overview of client/server coupling of LAMMPS with another code where
one code is the "client" and sends request messages to a "server"
code.  The server responds to each request with a reply message.  This
enables the two codes to work in tandem to perform a simulation.

When this command is invoked, LAMMPS will run in server mode in an
endless loop, waiting for messages from the client code.  The client
signals when it is done sending messages to LAMMPS, at which point the
loop will exit, and the remainder of the LAMMPS input script will be
processed.

The *protocol* argument defines the format and content of messages
that will be exchanged between the two codes.  The current options
are:

* :doc:`md <server_md>` = run dynamics with another code
* :doc:`mc <server_mc>` = perform Monte Carlo moves with another code

For protocol *md*\ , LAMMPS can be either a client (via the :doc:`fix client/md <fix_client_md>` command) or server.  See the :doc:`server md <server_md>` doc page for details on the protocol.

For protocol *mc*\ , LAMMPS can be the server.  See the :doc:`server mc <server_mc>` doc page for details on the protocol.

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

:doc:`message <message>`, :doc:`fix client/md <fix_client_md>`

**Default:** none
