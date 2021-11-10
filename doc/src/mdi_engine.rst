.. index:: mdi/engine

mdi_engine command
==================

Syntax
""""""

.. parsed-literal::

   mdi_engine

Description
"""""""""""

This command is used to have LAMMPS act as a server with another
client code to effectively couple the two codes together in
client/server mode.

More specifically, this command causes LAMMPS to begin using the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_
to run as an MDI engine (server), responding to commands made by an
external MDI driver code (client).  See the :doc:`Howto mdi
<Howto_mdi>` page for more information about how LAMMPS can work
as both an MDI driver or engine.

General information about launching codes that communicate using the
MDI Library can be found in the `corresponding page
<https://molssi-mdi.github.io/MDI_Library/html/library_page.html#library_launching_sec>`_
of the MDI Library's documentation.

----------

This command should typically be used in an input script after LAMMPS
has setup the system it is going to model in collaboration with the
driver code.  Depending on how the driver code tells the LAMMPS engine
to exit, other commands can be executed after this command, but
typically it should be used at the end of the LAMMPS input script.

To act as a MD-based MDI engine, this is the list of MDI commands from
a driver code which LAMMPS currently recognizes.  See more details
about these commands in the `MDI library documentation
<https://molssi-mdi.github.io/MDI_Library/html/mdi_standard.html>`_
.. NOTE: Taylor - is this the best link for this info?  Can we flesh this
.. out with the full list of supported commands?  Maybe the distinction
.. of what "node" the commands refer to is not needed in this table?

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Command name
     - Action
   * - >NATOMS
     - Driver sends the number of atoms in the system
   * - <NATOMS
     - Driver requests the number of atoms in the system
   * - <COORDS
     - Driver requests 3*N double-precision atom coordinates
   * - >FORCES
     - Driver sends 3*N double-precision atom forces
   * - <COORDS
     - Driver requests 3*N double-precision atom forces
   * - EXIT
     - Driver tells the engine (LAMMPS) to exit engine mode

If these commands are not sufficient to support what a driver which
you write needs, additional commands can be defined by simply using a
new command name not in this list.  Code to support the new command
needs to be added to the MDI package within LAMMPS; see its
src/MDI/mdi_engine.cpp and fix_mdi_engine.cpp files.

Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix mdi/engine <fix_mdi_engine>`

Default
"""""""

None
