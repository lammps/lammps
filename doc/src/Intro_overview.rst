Overview of LAMMPS
------------------

LAMMPS is a classical molecular dynamics (MD) code that models
ensembles of particles in a liquid, solid, or gaseous state.  It can
model atomic, polymeric, biological, solid-state (metals, ceramics,
oxides), granular, coarse-grained, or macroscopic systems using a
variety of interatomic potentials (force fields) and boundary
conditions.  It can model 2d or 3d systems with only a few particles
up to millions or billions.

LAMMPS can be built and run on a laptop or desktop machine, but is
designed for parallel computers.  It will run in serial and on any
parallel machine that supports the `MPI <mpi_>`_ message-passing
library.  This includes shared-memory boxes and distributed-memory
clusters and supercomputers. Parts of LAMMPS also support
`OpenMP multi-threading <omp_>`_, vectorization and GPU acceleration.

.. _mpi: https://en.wikipedia.org/wiki/Message_Passing_Interface
.. _lws: https://www.lammps.org
.. _omp: https://www.openmp.org

LAMMPS is written in C++ and requires a compiler that is at least
compatible with the C++-11 standard.  Earlier versions were written in
F77, F90, and C++-98.  See the `History page
<https://www.lammps.org/history.html>`_ of the website for details.  All
versions can be downloaded as source code from the `LAMMPS website
<lws_>`_.

LAMMPS is designed to be easy to modify or extend with new capabilities,
such as new force fields, atom types, boundary conditions, or
diagnostics.  See the :doc:`Modify <Modify>` page for more details.

In the most general sense, LAMMPS integrates Newton's equations of
motion for a collection of interacting particles.  A single particle
can be an atom or molecule or electron, a coarse-grained cluster of
atoms, or a mesoscopic or macroscopic clump of material.  The
interaction models that LAMMPS includes are mostly short-range in
nature; some long-range models are included as well.

LAMMPS uses neighbor lists to keep track of nearby particles.  The
lists are optimized for systems with particles that are repulsive at
short distances, so that the local density of particles never becomes
too large.  This is in contrast to methods used for modeling plasma
or gravitational bodies (e.g. galaxy formation).

On parallel machines, LAMMPS uses spatial-decomposition techniques with
MPI parallelization to partition the simulation domain into small
sub-domains of equal computational cost, one of which is assigned to
each processor.  Processors communicate and store "ghost" atom
information for atoms that border their sub-domain.  Multi-threading
parallelization and GPU acceleration with with particle-decomposition
can be used in addition.
