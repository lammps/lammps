.. index:: fix mdi/qm

fix mdi/qm command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mdi/qm keyword

* ID, group-ID are documented in :doc:`fix <fix>` command
* mdi/qm = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *plugin* or *add* of *every*

  .. parsed-literal::

       *plugin* args = none
       *add* args = *yes* or *no*
         yes = add returned value from server code to LAMMPS quantities
         no = do not add returned values to LAMMPS quantities
       *every* args = Nevery
         Nevery = request values from server code once every Nevery steps

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mdi/aimd
   fix 1 all mdi/aimd plugin

Description
"""""""""""

This command enables LAMMPS to act as a client with another server
code that will compute the total energy, per-atom forces, and total
virial for conformations and box size/shape that LAMMPS sends it.

Typically the server code will be a quantum code, hence the name of
the fix.  However this is not required, the server code could be
another classical molecular dynamics code or LAMMPS itself.  The
server code must support use of the `MDI Library
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_ as
explained below.

Example use cases for this fix are the following (discussed further
below):

* perform an ab intitio MD (AIMD) simulation with quantum forces
* perform an energy minimziation with quantum forces
* perform a nudged elatic band (NEB) calculation with quantum forces
* LAMMPS reads or creates a series of configurations and requests a QM calculation for each one

The code coupling performed by this command is done via the `MDI
Library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_.
LAMMPS runs as an MDI driver (client), and sends MDI commands to an
external MDI engine code (server), e.g. a quantum code which has
support for MDI.  See the :doc:`Howto mdi <Howto_mdi>` page for more
information about how LAMMPS can operate as either an MDI driver or
engine.

The examples/mdi directory contains input scripts performing AIMD in
this manner with LAMMPS acting as both a driver and an engine
(surrogate for a QM code).  The examples/mdi/README file explains how
to launch both driver and engine codes so that they communicate using
the MDI library via either MPI or sockets.  Any QM code that supports
MDI could be used in place of LAMMPS acting as a QM surrogate.  See
the :doc:`Howto mdi <Howto_mdi>` page for a current list (March 2022)
of such QM codes.

Engine codes can support MDI in either or both of two ways.  The
engine code can support being used as a stand-alone code, launched at
the same time as LAMMPS.  Or the engine code can support being used as
a plugin library, which LAMMPS loads.  See the :doc:`mdi plugin <mdi>`
command for how to trigger LAMMPS to load the plugin library.  The
examples/mdi/README file explains both use cases: launching both the
driver and engine as stand-alone codes or having the driver code load
an engine as a plugin library.

----------



To use this fix with a plugin engine, you must specify the
*plugin* keyword as the last argument, as illustrated above.

.. note::

   As of April 2022, the *plugin* keyword is needed.  In a future
   version of the MDI library it will no longer be necessary.

----------

----------

This fix performs the timestepping portion of anAIMD simulation.
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

Stored values by fix: energy, per-atoms forces, virial.



To use LAMMPS as an MDI driver in conjunction with other MDI-enabled
atomistic codes, the :doc:`units <units>` command should be used to
specify *real* or *metal* units.  This will ensure the correct unit
conversions between LAMMPS and MDI units.  The other code will also
perform similar unit conversions into its preferred units.

LAMMPS can also be used as an MDI driver in other unit choices it
supports, e.g. *lj*, but then no unit conversion is performed.

Related commands
""""""""""""""""

:doc:`mdi engine <mdi>`

Default
"""""""

The default for the optional keywords is add = yes, every = 1.

