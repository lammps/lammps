Coupling LAMMPS to other codes
==============================

LAMMPS is designed to allow it to be coupled to other codes.  For
example, a quantum mechanics code might compute forces on a subset of
atoms and pass those forces to LAMMPS.  Or a continuum finite element
(FE) simulation might use atom positions as boundary conditions on FE
nodal points, compute a FE solution, and return interpolated forces on
MD atoms.

LAMMPS can be coupled to other codes in at least 4 ways.  Each has
advantages and disadvantages, which you will have to think about in the
context of your application.

1. Define a new :doc:`fix <fix>` command that calls the other code.  In
   this scenario, LAMMPS is the driver code.  During its timestepping,
   the fix is invoked, and can make library calls to the other code,
   which has been linked to LAMMPS as a library.  This is the way the
   :ref:`POEMS <PKG-POEMS>` package that performs constrained rigid-body motion
   on groups of atoms is hooked to LAMMPS.  See the :doc:`fix poems
   <fix_poems>` command for more details.  See the :doc:`Modify
   <Modify>` doc pages for info on how to add a new fix to LAMMPS.

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

   See the :doc:`Modify command <Modify_command>` doc page for info on how
   to add a new command to LAMMPS.

.. spacer

3. Use LAMMPS as a library called by another code.  In this case the
   other code is the driver and calls LAMMPS as needed.  Or a wrapper
   code could link and call both LAMMPS and another code as libraries.
   Again, the :doc:`run <run>` command has options that allow it to be
   invoked with minimal overhead (no setup or clean-up) if you wish to
   do multiple short runs, driven by another program.  Details about
   using the library interface are given in the :doc:`library API
   <pg_library>` documentation.

.. spacer

4. Couple LAMMPS with another code in a client/server mode.  This is
   described on the :doc:`Howto client/server <Howto_client_server>` doc
   page.
