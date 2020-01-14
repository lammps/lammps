.. index:: fix client/md

fix client/md command
=====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID client/md

* ID, group-ID are documented in :doc:`fix <fix>` command
* client/md = style name of this fix command

Examples
""""""""


.. parsed-literal::

   fix 1 all client/md

Description
"""""""""""

This fix style enables LAMMPS to run as a "client" code and
communicate each timestep with a separate "server" code to perform an
MD simulation together.

The :doc:`Howto client/server <Howto_client_server>` doc page gives an
overview of client/server coupling of LAMMPS with another code where
one code is the "client" and sends request messages to a "server"
code.  The server responds to each request with a reply message.  This
enables the two codes to work in tandem to perform a simulation.

When using this fix, LAMMPS (as the client code) passes the current
coordinates of all particles to the server code each timestep, which
computes their interaction, and returns the energy, forces, and virial
for the interacting particles to LAMMPS, so it can complete the
timestep.

Note that the server code can be a quantum code, or another classical
MD code which encodes a force field (pair\_style in LAMMPS lingo) which
LAMMPS does not have.  In the quantum case, this fix is a mechanism
for running *ab initio* MD with quantum forces.

The group associated with this fix is ignored.

The protocol and :doc:`units <units>` for message format and content
that LAMMPS exchanges with the server code is defined on the :doc:`server md <server_md>` doc page.

Note that when using LAMMPS as an MD client, your LAMMPS input script
should not normally contain force field commands, like a
:doc:`pair\_style <pair_style>`, :doc:`bond\_style <bond_style>`, or
:doc:`kspace\_style <kspace_style>` command.  However it is possible for
a server code to only compute a portion of the full force-field, while
LAMMPS computes the remaining part.  Your LAMMPS script can also
specify boundary conditions or force constraints in the usual way,
which will be added to the per-atom forces returned by the server
code.

See the examples/message dir for example scripts where LAMMPS is both
the "client" and/or "server" code for this kind of client/server MD
simulation.  The examples/message/README file explains how to launch
LAMMPS and another code in tandem to perform a coupled simulation.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the potential energy computed by the server application to
the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

The :doc:`fix\_modify <fix_modify>` *virial* option is supported by this
fix to add the server application's contribution to the system's
virial as part of :doc:`thermodynamic output <thermo_style>`.  The
default is *virial yes*

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the potential
energy discussed above.  The scalar value calculated by this fix is
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the MESSAGE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

A script that uses this command must also use the
:doc:`message <message>` command to setup and shut down the messaging
protocol with the server code.

Related commands
""""""""""""""""

:doc:`message <message>`, :doc:`server <server>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
