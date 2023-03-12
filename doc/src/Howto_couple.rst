Coupling LAMMPS to other codes
==============================

LAMMPS is designed to support being coupled to other codes.  For
example, a quantum mechanics code might compute forces on a subset of
atoms and pass those forces to LAMMPS.  Or a continuum finite element
(FE) simulation might use atom positions as boundary conditions on FE
nodal points, compute a FE solution, and return interpolated forces on
MD atoms.

LAMMPS can be coupled to other codes in at least 4 different ways.  Each
has advantages and disadvantages, which you will have to think about in
the context of your application.

1. Define a new :doc:`fix <fix>` command that calls the other code.  In
   this scenario, LAMMPS is the driver code.  During timestepping, the
   fix is invoked, and can make library calls to the other code, which
   has been linked to LAMMPS as a library.  This is the way the
   :ref:`LATTE <PKG-LATTE>` package, which performs density-functional
   tight-binding calculations using the `LATTE software
   <https://github.com/lanl/LATTE>`_ to compute forces, is interfaced to
   LAMMPS.  See the :doc:`fix latte <fix_latte>` command for more
   details.  Also see the :doc:`Modify <Modify>` pages for information
   on how to add a new fix to LAMMPS.

.. spacer

2. Define a new LAMMPS command that calls the other code.  This is
   conceptually similar to method (1), but in this case LAMMPS and the
   other code are on a more equal footing.  Note that now the other code
   is not called during the timestepping of a LAMMPS run, but between
   runs.  The LAMMPS input script can be used to alternate LAMMPS runs
   with calls to the other code, invoked via the new command.  The
   :doc:`run <run>` command facilitates this with its *every* option,
   which makes it easy to run a few steps, invoke the command, run a few
   steps, invoke the command, etc.

   In this scenario, the other code can be called as a library, as in
   1., or it could be a stand-alone code, invoked by a system() call
   made by the command (assuming your parallel machine allows one or
   more processors to start up another program).  In the latter case the
   stand-alone code could communicate with LAMMPS through files that the
   command writes and reads.

   See the :doc:`Modify command <Modify_command>` page for information
   on how to add a new command to LAMMPS.

.. spacer

3. Use LAMMPS as a library called by another code.  In this case, the
   other code is the driver and calls LAMMPS as needed.  Alternately, a
   wrapper code could link and call both LAMMPS and another code as
   libraries.  Again, the :doc:`run <run>` command has options that
   allow it to be invoked with minimal overhead (no setup or clean-up)
   if you wish to do multiple short runs, driven by another program.
   Details about using the library interface are given in the
   :doc:`library API <Library>` documentation.

.. spacer

4. Couple LAMMPS with another code in a client/server fashion, using the
   `MDI Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_
   developed by the `Molecular Sciences Software Institute (MolSSI)
   <https://molssi.org>`_ to run LAMMPS as either an MDI driver (client)
   or an MDI engine (server).  The MDI driver issues commands to the MDI
   server to exchange data between them.  See the :doc:`Howto_mdi` page for
   more information about how LAMMPS can operate in either of these modes.
