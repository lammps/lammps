Using LAMMPS with the MDI library for code coupling
===================================================

Client/server coupling of two (or more) codes is where one code is the
"client" and sends request messages (data) to one (or more) "server"
code(s).  A server responds to each request with a reply message
(data).  This enables two (or more) codes to work in tandem to perform
a simulation.  In this context, LAMMPS can act as either a client or
server code.  It does this by using the `MolSSI Driver Interface (MDI)
library <https://molssi-mdi.github.io/MDI_Library/html/index.html>`_,
developed by the `Molecular Sciences Software Institute (MolSSI)
<https://molssi.org>`_, which is supported by the :ref:`MDI <PKG-MDI>`
package.

Alternate methods for coupling codes with LAMMPS are described on the
:doc:`Howto_couple` page.

Some advantages of client/server coupling are that the codes can run
as stand-alone executables; they need not be linked together.  Thus,
neither code needs to have a library interface.  This also makes it
easy to run the two codes on different numbers of processors.  If a
message protocol (format and content) is defined for a particular kind
of simulation, then in principle any code which implements the
client-side protocol can be used in tandem with any code which
implements the server-side protocol.  Neither code needs to know what
specific other code it is working with.

In MDI nomenclature, a client code is the "driver", and a server code is
an "engine".  One driver code can communicate with one or more instances
of one or more engine codes.  Driver and engine codes can be written in
any language: C, C++, Fortran, Python, etc.

In addition to allowing driver and engine(s) to run as stand-alone
executables, MDI also enables an engine to be a *plugin* to the client
code.  In this scenario, server code(s) are compiled as shared
libraries, and one (or more) instances of the server are instantiated
by the driver code.  If the driver code runs in parallel, it can split
its MPI communicator into multiple sub-communicators, and launch each
plugin engine instance on a sub-communicator.  Driver processors
within that sub-communicator exchange messages with the corresponding
engine instance, and can also send MPI messages to other processors in
the driver.  The driver code can also destroy engine instances and
re-instantiate them.  LAMMPS can operate as either a stand-alone or
plugin MDI engine.  When it operates as a driver, it can use either
stand-alone or plugin MDI engines.

The way in which an MDI driver communicates with an MDI engine is by
making MDI_Send() and MDI_Recv() calls, which are conceptually similar
to MPI_Send() and MPI_Recv() calls.  Each send or receive operation
uses a string to identify the command name, and optionally some data,
which can be a single value or vector of values of any data type.
Inside the MDI library, data is exchanged between the driver and
engine via MPI calls or sockets.  This is a run-time choice by the user.

----------

The :ref:`MDI <PKG-MDI>` package provides a :doc:`mdi engine <mdi>`
command, which enables LAMMPS to operate as an MDI engine.  Its doc
page explains the variety of standard and custom MDI commands which
the LAMMPS engine recognizes and can respond to.

The package also provides a :doc:`mdi plugin <mdi>` command, which
enables LAMMPS to operate as an MDI driver and load an MDI engine as a
plugin library.

The package furthermore includes a :doc:`fix mdi/qm <fix_mdi_qm>`
command, in which LAMMPS operates as an MDI driver in conjunction with a
quantum mechanics code as an MDI engine.  The post_force() method of the
``fix_mdi_qm.cpp`` file shows how a driver issues MDI commands to
another code.  This command can be used to couple to an MDI engine,
which is either a stand-alone code or a plugin library.

As explained in the :doc:`fix mdi/qm <fix_mdi_qm>` command
documentation, it can be used to perform *ab initio* MD simulations or
energy minimizations, or to evaluate the quantum energy and forces for a
series of independent systems.  The ``examples/mdi`` directory has
example input scripts for all of these use cases.

The package also has a :doc:`fix mdi/qmmm <fix_mdi_qmmm>` command in
which LAMMPS operates as an MDI driver in conjunction with a quantum
mechanics code as an MDI engine to perform QM/MM simulations.  The
LAMMPS input script partitions the system into QM and MM (molecular
mechanics) atoms.  As described below the ``examples/QUANTUM`` directory
has examples for coupling to 3 different quantum codes in this manner.

----------

The ``examples/mdi`` directory contains Python scripts and LAMMPS input
script which use LAMMPS as either an MDI driver or engine, or both.
Currently, 5 example use cases are provided:

* Run ab initio MD (AIMD) using 2 instances of LAMMPS.  As a driver,
  LAMMPS performs the timestepping in either NVE or NPT mode.  As an
  engine, LAMMPS computes forces and is a surrogate for a quantum
  code.

* LAMMPS runs an MD simulation as a driver.  Every N steps it passes the
  current snapshot to an MDI engine to evaluate the energy, virial, and
  peratom forces.  As the engine, LAMMPS is a surrogate for a quantum
  code.

* LAMMPS loops over a series of data files and passes the configuration
  to an MDI engine to evaluate the energy, virial, and peratom forces
  and thus acts as a simulation driver.  As the engine, LAMMPS is used
  as a surrogate for a quantum code.

* A Python script driver invokes a sequence of unrelated LAMMPS
  calculations.  Calculations can be single-point energy/force
  evaluations, MD runs, or energy minimizations.

* Run AIMD with a Python driver code and 2 LAMMPS instances as engines.
  The first LAMMPS instance performs MD timestepping.  The second LAMMPS
  instance acts as a surrogate QM code to compute forces.

.. note::

   In any of these examples where LAMMPS is used as an engine, an actual
   QM code (provided it has support for MDI) could be used in its place,
   without modifying the input scripts or launch commands, except to
   specify the name of the QM code.

The ``examples/mdi/Run.sh`` file illustrates how to launch both driver
and engine codes so that they communicate using the MDI library via
either MPI or sockets, or using the engine as a stand-alone code, or
as a plugin library.

-------------

As of March 2023, these are quantum codes with MDI support provided via
Python wrapper scripts included in the LAMMPS distribution.  These can
be used with the fix mdi/qm and fix mdi/qmmm commands to perform QM
calculations of an entire system (e.g. AIMD) or QM/MM simulations.  See
the ``examples/QUANTUM`` sub-directories for more details:

* LATTE - AIMD only
* PySCF - QM/MM only
* NWChem - AIMD or QM/MM

There are also at least two quantum codes which have direct MDI
support, `Quantum ESPRESSO (QE) <https://www.quantum-espresso.org/>`_
and `INQ <https://qsg.llnl.gov/node/101.html>`_.  There are also
several QM codes which have indirect support through QCEngine or i-PI.
The former means they require a wrapper program (QCEngine) with MDI
support which writes/read files to pass data to the quantum code
itself.  The list of QCEngine-supported and i-PI-supported quantum
codes is on the `MDI webpage
<https://molssi-mdi.github.io/MDI_Library/html/index.html>`_.

These direct- and indirect-support codes should be usable for full
system calculations (e.g. AIMD).  Whether they support QM/MM models
depends on the individual QM code.
