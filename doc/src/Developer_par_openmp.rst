.. _openmp_parallelism:

OpenMP Parallelism
^^^^^^^^^^^^^^^^^^

The styles in the INTEL, KOKKOS, and OPENMP packages offer to use OpenMP
thread parallelism to predominantly distribute loops over local data and
thus follow an orthogonal parallelization strategy to the decomposition
into spatial domains used by the :doc:`MPI partitioning
<Developer_par_part>`.  For clarity, this section discusses only the
implementation in the OPENMP package, as it is the simplest. The INTEL
and KOKKOS packages offer additional options and are significantly more
complex since they support more features like vectorization, accelerator
support and a choice of using double, single, and mixed precision
floating-point math.

A practical discussion of the actual implementation with source code
examples can be found in the manual section on :doc:`writing styles for
the OPENMP package <Developer_write_openmp>`.

.. _openmp_rationale:

Rationale for OpenMP Integration
""""""""""""""""""""""""""""""""

While LAMMPS is already well-parallelized using MPI-based domain
decomposition, OpenMP provides additional benefits, particularly on
modern multi-core hardware and at high node counts:

- **Efficient Resource Utilization**: On multi-core machines, MPI tasks
  can be communication-bandwidth limited because all MPI communication
  has to share the same network link for communication between nodes.
  Running fewer MPI tasks and using OpenMP threads can improve
  performance by reducing communication overhead and memory contention.
  However, the in-node memory bandwidth demands for multi-threading are
  higher, since domain decomposition promotes CPU cache locality, while
  the multi-thread implementation uses replicated per-thread storage to
  avoid data races. This has the additional overhead or requiring a
  reduction, which is most efficient for small numbers of threads.
- **PPPM Scaling**: Long-range electrostatics solvers like PPPM have
  scaling limitations with MPI communication at high node counts because
  the 3d-FFTs require all-to-all communications and with smaller
  domains the ratio of computation versus communication become less
  favorable. Using OpenMP threads can allow running these solvers on a
  subset of MPI tasks, improving overall scaling.
- **Parallelization Granularity**: OpenMP allows parallelization over
  particles/neighbors within a node, providing a finer grain of parallelism
  than domain decomposition alone.
- **Load Balancing Considerations**: The domain decomposition based
  MPI parallelization implicitly assumes that the simulated system is
  homogeneous and thus roughly the same amount of work needs to be done
  by each MPI rank. However, for inhomogeneous systems, this is not
  always the case and the available load-balancing options have limitations, too.
  By shifting the parallelization to OpenMP, the subdomains per MPI
  rank become larger and load balancing is usually more effective then.
- **Capability Computing**: Hybrid MPI+OpenMP is often essential for
  achieving optimal performance on large HPC clusters or supercomputers.

.. _openmp_design:

Design of the OPENMP Package
""""""""""""""""""""""""""""

One of the key decisions when implementing the OPENMP package was to
keep the changes to the source code small, so that it would be easier to
maintain the code and keep it in sync with the non-threaded standard
implementation.  This is achieved by:

* **Inheritance**: Making the OPENMP version a derived class from the
  regular version (e.g. ``PairLJCutOMP`` from ``PairLJCut``) and only
  overriding methods that are multi-threaded or need to be modified
  (similar to what was done before in the OPT package).
* **Minimal Modification**: Keeping the structure in the modified code
  very similar so that side-by-side comparisons are still useful.
* **Helper Classes**: Offloading multi-thread support functions into
  three separate classes:

  - ``ThrOMP``: Provides multi-thread aware functionality like "_thr"
    variants of tally functions via multiple inheritance.
  - ``ThrData``: Manages per-thread data structures to avoid "false
    sharing" slowdowns.
  - ``FixOMP``: Manages the global multi-thread state and settings
    activated by the :doc:`package omp <package>` command.

.. _openmp_avoid_races:

Avoiding data races
"""""""""""""""""""

A key problem when implementing thread parallelism in an MD code is to
avoid data races when updating accumulated per-atom properties like
forces, energies, stresses, and torques.  When interactions are
computed, they always involve multiple atoms and thus there are race
conditions when multiple threads want to update per-atom data of the
same atoms.

Without a significant rewrite of the code, there are three main approaches
to avoid data races.

* Use locks or **atomic operations** to guarantee that only one thread
  at a time updates per-thread data.  This is most commonly used when
  the access conflicts are rare.
* Compute per-atom properties **multiple times** for all threads that
  "own" one or more atoms of a tuple and store only the data for those
  "owned" atom.  This is most effective for a large number of threads
  and thus commonly used with GPU acceleration.
* Have **per-thread copies** of the per-atom data and update them
  independently and use a reduction to combine the per-thread data
  into the regular per-atom storage.  This requires additional steps
  to set up, manage, clear, and combine the collected data.

The OPENMP package uses **replicated per-thread data structures** because
this approach:

* retains the performance for the single-thread case unlike the other options,
* keeps the code maintainable and similar to the non-threaded version,
* and is most efficient for a small number of threads (2-8), which is a
  common use case, since OpenMP is typically the secondary parallelization
  option after MPI.

.. _openmp_scheduling:

Loop scheduling
"""""""""""""""

Multi-thread parallelization is applied by distributing (outer) loops
statically across threads.  Typically, this is the loop over local atoms
*i* when processing *i,j* pairs of atoms from a neighbor list or the
list of bonds, angles, dihedrals, or impropers.  Since neighbor lists
typically result in a similar number of neighbors per atom (for
homogeneous systems), significant load imbalances across threads are
uncommon.

.. _openmp_neighbor:

Neighbor list parallelization
"""""""""""""""""""""""""""""

In addition to force computations, neighbor list generation is
parallelized. Each thread operates on a different chunk of "owned"
atoms and manages its own set of neighbor "pages" using an instance
of the :cpp:class:`MyPage <LAMMPS_NS::MyPage>` page allocator.

