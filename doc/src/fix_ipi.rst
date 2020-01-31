.. index:: fix ipi

fix ipi command
===============

Syntax
""""""


.. parsed-literal::

   fix ID group-ID ipi address port [unix] [reset]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ipi = style name of this fix command
* address = internet address (FQDN or IP), or UNIX socket name
* port = port number (ignored for UNIX sockets)
* optional keyword = *unix*\ , if present uses a unix socket
* optional keyword = *reset*\ , if present reset electrostatics at each call

Examples
""""""""

fix 1 all ipi my.server.com 12345
fix 1 all ipi mysocket 666 unix reset

Description
"""""""""""

This fix enables LAMMPS to be run as a client for the i-PI Python
wrapper :ref:`(IPI) <ipihome>` for performing a path integral molecular dynamics
(PIMD) simulation.  The philosophy behind i-PI is described in the
following publication :ref:`(IPI-CPC) <IPICPC>`.

A version of the i-PI package, containing only files needed for use
with LAMMPS, is provided in the tools/i-pi directory.  See the
tools/i-pi/manual.pdf for an introduction to i-PI.  The
examples/USER/i-pi directory contains example scripts for using i-PI
with LAMMPS.

In brief, the path integral molecular dynamics is performed by the
Python wrapper, while the client (LAMMPS in this case) simply computes
forces and energy for each configuration. The communication between
the two components takes place using sockets, and is reduced to the
bare minimum. All the parameters of the dynamics are specified in the
input of i-PI, and all the parameters of the force field must be
specified as LAMMPS inputs, preceding the *fix ipi* command.

The server address must be specified by the *address* argument, and
can be either the IP address, the fully-qualified name of the server,
or the name of a UNIX socket for local, faster communication. In the
case of internet sockets, the *port* argument specifies the port
number on which i-PI is listening, while the *unix* optional switch
specifies that the socket is a UNIX socket.

Note that there is no check of data integrity, or that the atomic
configurations make sense. It is assumed that the species in the i-PI
input are listed in the same order as in the data file of LAMMPS. The
initial configuration is ignored, as it will be substituted with the
coordinates received from i-PI before forces are ever evaluated.

A note of caution when using potentials that contain long-range
electrostatics, or that contain parameters that depend on box size:
all of these options will be initialized based on the cell size in the
LAMMPS-side initial configuration and kept constant during the run.
This is required to e.g. obtain reproducible and conserved forces.
If the cell varies too wildly, it may be advisable to re-initialize
these interactions at each call. This behavior can be requested by
setting the *reset* switch.

**Restart, fix\_modify, output, run start/stop, minimize info:**

There is no restart information associated with this fix, since all
the dynamical parameters are dealt with by i-PI.

Restrictions
""""""""""""


Using this fix on anything other than all atoms requires particular
care, since i-PI will know nothing on atoms that are not those whose
coordinates are transferred. However, one could use this strategy to
define an external potential acting on the atoms that are moved by
i-PI.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.  Because of the
use of UNIX domain sockets, this fix will only work in a UNIX
environment.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`


----------


.. _IPICPC:



**(IPI-CPC)** Ceriotti, More and Manolopoulos, Comp Phys Comm, 185,
1019-1026 (2014).

.. _ipihome:



**(IPI)**
`http://epfl-cosmo.github.io/gle4md/index.html?page=ipi <http://epfl-cosmo.github.io/gle4md/index.html?page=ipi>`_


