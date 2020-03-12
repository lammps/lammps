.. index:: server md

server md command
=================

Syntax
""""""

.. parsed-literal::

   server md

md = the protocol argument to the :doc:`server <server>` command

Examples
""""""""

.. code-block:: LAMMPS

   server md

Description
"""""""""""

This command starts LAMMPS running in "server" mode, where it will
expect messages from a separate "client" code that match the *md*
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
in server mode.  See an example of how this command is used in
examples/message/in.message.server.

----------

When using this command, LAMMPS (as the server code) receives the
current coordinates of all particles from the client code each
timestep, computes their interaction, and returns the energy, forces,
and pressure for the interacting particles to the client code, so it
can complete the timestep.  This command could also be used with a
client code that performs energy minimization, using the server to
compute forces and energy each iteration of its minimizer.

When using the :doc:`fix client/md <fix_client_md>` command, LAMMPS (as
the client code) does the timestepping and receives needed energy,
forces, and pressure values from the server code.

The format and content of the exchanged messages are explained here in
a conceptual sense.  Python-style pseudo code for the library calls to
the CSlib is shown, which performs the actual message exchange between
the two codes.  See the `CSlib website <https://cslib.sandia.gov>`_ doc
pages for more details on the actual library syntax.  The "cs" object
in this pseudo code is a pointer to an instance of the CSlib.

See the src/MESSAGE/server\_md.cpp and src/MESSAGE/fix\_client\_md.cpp
files for details on how LAMMPS uses these messages.  See the
examples/COUPLE/lammps\_vasp/vasp\_wrap.py or
examples/COUPLE/lammps\_nwchem/nwchem\_wrap.py files for examples of how
a quantum code (VASP or NWChem) can use these messages.

The following pseudo-code uses these values, defined as enums.

Define:

.. parsed-literal::

   SETUP=1, STEP=2
   DIM=1, PERIODICITY=2, ORIGIN=3, BOX=4, NATOMS=5, NTYPES=6, TYPES=7, COORDS=8, UNITS-9, CHARGE=10
   FORCES=1, ENERGY=2, PRESSURE=3, ERROR=4

**Client sends 2 kinds of messages**\ :

.. parsed-literal::

   # required fields: DIM, PERIODICTY, ORIGIN, BOX, NATOMS, NTYPES, TYPES, COORDS
   # optional fields: UNITS, CHARGE

   cs->send(SETUP,nfields)        # msgID with nfields

   cs->pack_int(DIM,dim)          # dimension (2,3) of simulation
   cs->pack(PERIODICITY,3,xyz)    # periodicity flags in 3 dims
   cs->pack(ORIGIN,3,origin)      # lower-left corner of simulation box
   cs->pack(BOX,9,box)            # 3 edge vectors of simulation box
   cs->pack_int(NATOMS,natoms)    # total number of atoms
   cs->pack_int(NTYPES,ntypes)    # number of atom types
   cs->pack(TYPES,natoms,type)    # vector of per-atom types
   cs->pack(COORDS,3\*natoms,x)    # vector of 3N atom coords
   cs->pack_string(UNITS,units)   # units = "lj", "real", "metal", etc
   cs->pack(CHARGE,natoms,q)      # vector of per-atom charge

   # required fields: COORDS
   # optional fields: ORIGIN, BOX

   cs->send(STEP,nfields)         # msgID with nfields

   cs->pack(COORDS,3\*natoms,x)    # vector of 3N atom coords
   cs->pack(ORIGIN,3,origin)      # lower-left corner of simulation box
   cs->pack(BOX,9,box)            # 3 edge vectors of simulation box

**Server replies to either kind of message**\ :

.. parsed-literal::

   # required fields: FORCES, ENERGY, PRESSURE
   # optional fields: ERROR

   cs->send(msgID,nfields)      # msgID with nfields
   cs->pack(FORCES,3\*Natoms,f)  # vector of 3N forces on atoms
   cs->pack(ENERGY,1,poteng)    # total potential energy of system
   cs->pack(PRESSURE,6,press)   # global pressure tensor (6-vector)
   cs->pack_int(ERROR,flag)     # server had an error (e.g. DFT non-convergence)

----------

The units for various quantities that are sent and received iva
messages are defined for atomic-scale simulations in the table below.
The client and server codes (including LAMMPS) can use internal units
different than these (e.g. :doc:`real units <units>` in LAMMPS), so long
as they convert to these units for messaging.

* COORDS, ORIGIN, BOX = Angstroms
* CHARGE = multiple of electron charge (1.0 is a proton)
* ENERGY = eV
* FORCES = eV/Angstrom
* PRESSURE = bars

Note that these are :doc:`metal units <units>` in LAMMPS.

If you wish to run LAMMPS in another its non-atomic units, e.g. :doc:`lj units <units>`, then the client and server should exchange a UNITS
message as indicated above, and both the client and server should
agree on the units for the data they exchange.

----------

Restrictions
""""""""""""

This command is part of the MESSAGE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`message <message>`, :doc:`fix client/md <fix_client_md>`

**Default:** none
