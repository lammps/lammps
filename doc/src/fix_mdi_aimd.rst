.. index:: fix mdi/aimd

fix mdi/aimd command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mdi/aimd keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* mdi/aimd = style name of this fix command
* optional keyword = *plugin*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mdi/aimd
   fix 1 all mdi/aimd plugin

Description
"""""""""""

This command enables LAMMPS to act as a client with another server
code to couple the two codes together to perform ab initio MD (AIMD)
simulations.

More specifically, this command causes LAMMPS to begin using the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_
to run as an MDI driver (client), which sends MDI commands to an
external MDI engine code (server) which in the case of AIMD is a
quantum mechanics (QM) code, or could be LAMMPS itself, acting as a
surrogate for a QM code.  See the :doc:`Howto mdi <Howto_mdi>` page
for more information about how LAMMPS can operate as either an MDI
driver or engine.

The examples/mdi directory contains input scripts performing AIMD in
this manner with LAMMPS acting as both a driver and an engine
(surrogate for a QM code).  The examples/mdi/README file explains how
to launch both driver and engine codes so that they communicate using
the MDI library via either MPI or sockets.  Any QM code that supports
MDI could be used in place of LAMMPS acting as a QM surrogate.  See
the :doc:`Howto mdi <Howto_mdi>` page for a current list (March 2022)
of such QM codes.

The engine code can run either as a stand-alone code, launched at the
same time as LAMMPS, or as a plugin library.  See the :doc:`mdi plugin
<mdi>` command for how to trigger LAMMPS to load the plugin library.
Again, the examples/mdi/README file explains how to launch both driver
and engine codes so that engine is used in plugin mode.

To use this fix with a plugin engine, you must specify the
*plugin* keyword as the last argument, as illustrated above.

.. note::

   As of April 2022, the *plugin* keyword is needed.  In a future
   version of the MDI library it will no longer be necessary.

----------

This fix performs the timestepping portion of an AIMD simulation.
Both LAMMPS and the engine code (QM or LAMMPS) should define the same
system (simulation box, atoms and their types) in their respective
input scripts.  LAMMPS then begins its timestepping.

At the point in each timestep when LAMMPS needs the force on each
atom, it communicates with the engine code.  It sends the current
simulation box size and shape (if they change dynamically, e.g. during
an NPT simulation), and the current atom coordinates.  The engine code
computes quantum forces on each atom and returns them to LAMMPS.  If
LAMMPS also needs the system energy and/or virial, it requests those
values from the engine code as well.

Restrictions
""""""""""""

This command is part of the MDI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To use LAMMPS as an MDI driver in conjunction with other MDI-enabled
atomistic codes, the :doc:`units <units>` command should be used to
specify *real* or *metal* units.  This will ensure the correct unit
conversions between LAMMPS and MDI units, which the other codes will
also perform in their preferred units.

LAMMPS can also be used as an MDI driver in other unit choices it
supports, e.g. *lj*, but then no unit conversion is performed.

Related commands
""""""""""""""""

:doc:`mdi engine <mdi>`

Default
"""""""

none
