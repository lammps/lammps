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

----------

(1) Define a new :doc:`fix <fix>` command that calls the other code.  In
this scenario, LAMMPS is the driver code.  During its timestepping,
the fix is invoked, and can make library calls to the other code,
which has been linked to LAMMPS as a library.  This is the way the
`POEMS <poems_>`_ package that performs constrained rigid-body motion on
groups of atoms is hooked to LAMMPS.  See the :doc:`fix poems <fix_poems>` command for more details.  See the
:doc:`Modify <Modify>` doc pages for info on how to add a new fix to
LAMMPS.

.. _poems: http://www.rpi.edu/~anderk5/lab

----------

(2) Define a new LAMMPS command that calls the other code.  This is
conceptually similar to method (1), but in this case LAMMPS and the
other code are on a more equal footing.  Note that now the other code
is not called during the timestepping of a LAMMPS run, but between
runs.  The LAMMPS input script can be used to alternate LAMMPS runs
with calls to the other code, invoked via the new command.  The
:doc:`run <run>` command facilitates this with its *every* option, which
makes it easy to run a few steps, invoke the command, run a few steps,
invoke the command, etc.

In this scenario, the other code can be called as a library, as in
(1), or it could be a stand-alone code, invoked by a system() call
made by the command (assuming your parallel machine allows one or more
processors to start up another program).  In the latter case the
stand-alone code could communicate with LAMMPS through files that the
command writes and reads.

See the :doc:`Modify command <Modify_command>` doc page for info on how
to add a new command to LAMMPS.

----------

(3) Use LAMMPS as a library called by another code.  In this case the
other code is the driver and calls LAMMPS as needed.  Or a wrapper
code could link and call both LAMMPS and another code as libraries.
Again, the :doc:`run <run>` command has options that allow it to be
invoked with minimal overhead (no setup or clean-up) if you wish to do
multiple short runs, driven by another program.

Examples of driver codes that call LAMMPS as a library are included in
the examples/COUPLE directory of the LAMMPS distribution; see
examples/COUPLE/README for more details:

* simple: simple driver programs in C++ and C which invoke LAMMPS as a
  library
* plugin: simple driver program in C which invokes LAMMPS as a plugin
  from a shared library.
* lammps_quest: coupling of LAMMPS and `Quest <quest_>`_, to run classical
  MD with quantum forces calculated by a density functional code
* lammps_spparks: coupling of LAMMPS and `SPPARKS <spparks_>`_, to couple
  a kinetic Monte Carlo model for grain growth using MD to calculate
  strain induced across grain boundaries

.. _quest: http://dft.sandia.gov/Quest

.. _spparks: http://www.sandia.gov/~sjplimp/spparks.html

The :doc:`Build basics <Build_basics>` doc page describes how to build
LAMMPS as a library.  Once this is done, you can interface with LAMMPS
either via C++, C, Fortran, or Python (or any other language that
supports a vanilla C-like interface).  For example, from C++ you could
create one (or more) "instances" of LAMMPS, pass it an input script to
process, or execute individual commands, all by invoking the correct
class methods in LAMMPS.  From C or Fortran you can make function
calls to do the same things.  See the :doc:`Python <Python_head>` doc
pages for a description of the Python wrapper provided with LAMMPS
that operates through the LAMMPS library interface.

The files src/library.cpp and library.h contain the C-style interface
to LAMMPS.  See the :doc:`Howto library <Howto_library>` doc page for a
description of the interface and how to extend it for your needs.

Note that the lammps_open() function that creates an instance of
LAMMPS takes an MPI communicator as an argument.  This means that
instance of LAMMPS will run on the set of processors in the
communicator.  Thus the calling code can run LAMMPS on all or a subset
of processors.  For example, a wrapper script might decide to
alternate between LAMMPS and another code, allowing them both to run
on all the processors.  Or it might allocate half the processors to
LAMMPS and half to the other code and run both codes simultaneously
before syncing them up periodically.  Or it might instantiate multiple
instances of LAMMPS to perform different calculations.

----------

(4) Couple LAMMPS with another code in a client/server mode.  This is
described on the :doc:`Howto client/server <Howto_client_server>` doc
page.
