Creating styles for the OPENMP package
-----------------------------------------

The OPENMP package in LAMMPS provides an accelerated execution of many
LAMMPS styles with multi-threading using `OpenMP pragmas
<https://www.openmp.org/>`_.  Adding support for a new style to the
OPENMP package involves creating a derived class that overrides the
standard serial method.

For a high-level discussion on the rationale and parallelization
strategy of the OpenMP package, see :doc:`OpenMP Parallelism
<Developer_par_openmp>`.

In the text below, we will describe important steps of the procedure of
adding multi-thread style to the OPENMP package by looking at a few
selected example sections of source files in the OPENMP package.

Package and build system considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The OPENMP package is one of the accelerator packages. New styles should
be added to the `src/OPENMP` directory. The contributed code needs to
support both the CMake build process and (if applicable) the traditional
GNU make process.

Case study: a pairwise additive model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many pair styles, adding OpenMP support involves creating a derived
class that inherits from the serial pair style class. Let's look at the
`pair_lj_cut` style and its OpenMP-accelerated variant.

Header file
"""""""""""

The header file for an OPENMP-enabled pair style (e.g., `pair_lj_cut_omp.h`)
typically derives from the original style class (e.g., `PairLJCut`).

.. code-block:: c++

   #ifndef LMP_PAIR_LJ_CUT_OMP_H
   #define LMP_PAIR_LJ_CUT_OMP_H

   #include "pair_lj_cut.h"
   #include "thr_omp.h"

   namespace LAMMPS_NS {

   class PairLJCutOMP : public PairLJCut, public ThrOMP {
    public:
     PairLJCutOMP(class LAMMPS *);
     virtual ~PairLJCutOMP() {}

     void compute(int, int) override;
     double memory_usage() override;

    private:
     template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
     void eval(int ifrom, int ito, ThrData *const thr);
   };
   }

The class inherits from `PairLJCut` (the base serial style) and `ThrOMP`,
a utility class that manages the OpenMP-specific data, such as thread-local
force and energy accumulation arrays.

Constructor and methods
^^^^^^^^^^^^^^^^^^^^^^^

The constructor for an OMP style typically just calls the base class
constructor and initializes the `ThrOMP` base:

.. code-block:: c++

   PairLJCutOMP::PairLJCutOMP(LAMMPS *lmp) :
     PairLJCut(lmp), ThrOMP(lmp, THR_PAIR)
   {
     respa_enable = 0;
   }

- The `THR_PAIR` argument tells the `ThrOMP` base class that this is a
  pair-style acceleration, which influences how data is allocated and
  managed.
- Setting `respa_enable = 0` is common for OMP accelerated styles as
  RESPA acceleration often requires separate implementation logic.

Why only few methods?
"""""""""""""""""""""

You may notice that only a few methods are overridden. Most methods
(e.g., `coeff()`, `init_one()`, `write_restart()`) are inherited directly
from `PairLJCut`. This is possible because:

* The serial class already handles the data structures and input parsing.
* The `compute()` method is the part that consumes most of the CPU time and is called in every step during a run.
* `single()` is not overridden because it cannot be efficiently
  parallelized with OpenMP.  It represents just one iteration of the
  inner loop of the main loop constructs of the `compute()` function.

Utility Functions for Threading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `ThrOMP` class and related LAMMPS infrastructure provide several key methods to manage parallel execution:

* `loop_setup_thr(int nlocal)`: This method calculates the range of atoms (or neighbors)
  each thread should process based on the number of threads and the total count.
  It ensures an even workload distribution across threads.
* `FixOMP::get_thr(int tid)`: This method retrieves the `ThrData` object for a specific
  thread ID (`tid`). `ThrData` maintains thread-private buffers for forces and energies,
  avoiding race conditions during the force computation.
* `reduce_thr(ThrData* thr, int eflag, int vflag)`: After the parallel region, this
  method is used to perform the reduction of per-thread force, energy, and virial
  contributions into the global arrays (e.g., `atom->f`).

`ev_setup_thr()` and `ev_tally_thr()`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These methods replace their non-OpenMP counterparts (`ev_setup()` and `ev_tally()`)
to support thread-local data:

* `ev_setup_thr(eflag, vflag, thr)`: Unlike `ev_setup()`, which initializes global
  energy/virial flags, `ev_setup_thr()` initializes the thread-specific `ThrData`
  structures to track energy and virial contributions independently per thread.
* `ev_tally_thr(i, j, nlocal, newton, evdwl, ecoul, fpair, delx, dely, delz, thr)`:
  Instead of updating global force/energy arrays directly, `ev_tally_thr()` updates
  the per-thread tally arrays stored within the provided `ThrData` object (`thr`).
  This allows multiple threads to record force/energy contributions simultaneously
  without locking or memory contention.

Performance Optimization with Templates and `_noalias`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The template-based approach in `eval<EVFLAG, EFLAG, NEWTON_PAIR>()` is a powerful
optimization technique:

* **Branch Elimination**: By passing state as template parameters, the compiler
  can generate specialized machine code for each combination of flags (e.g.,
  `NEWTON_PAIR=1`, `EFLAG=0`). This removes conditional branches (`if` statements)
  from the inner loop, which is critical for performance because CPU branch
  predictors often struggle with these toggles in tight loops.
* **`_noalias`**: The `_noalias` keyword (or compiler hint) tells the compiler
  that pointers (e.g., force arrays) do not point to the same memory addresses
  as other pointers used in the same context. This enables the compiler to
  aggressively optimize load/store operations and loop unrolling, as it
  no longer needs to ensure memory consistency across potential pointer aliasing.

Porting complex many-body styles (e.g., PairSWOMP, PairTersoffOMP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many-body pair styles like Stillinger-Weber (`PairSW`) or Tersoff (`PairTersoff`)
differ significantly from pairwise additive potentials. In pairwise additive
potentials (like LJ), the force on atom `i` depends only on its neighbors `j`.
Many-body styles, however, compute forces based on three-body (or higher-order)
terms, meaning the force on atom `i` depends on the environment of atom `j` as well.

Differences from Pairwise Additive Potentials
"""""""""""""""""""""""""""""""""""""""""""""

1. **Neighbor List Requirements**: Pairwise additive potentials typically
   require a "half" neighbor list. Many-body potentials, however, usually
   require a "full" neighbor list to iterate over all neighbors of both atoms `i`
   and `j` when calculating the three-body interaction `i-j-k`.
2. **Computational Structure**: Instead of a simple nested loop, many-body
   potentials involve complex kernels that iterate over neighbors of neighbors.
   This significantly increases the complexity of the parallel execution
   and data locality requirements.
3. **Communication**: Due to the three-body nature, the potential needs access
   to the positions of all atoms within the cutoff of `j`, even if `j` is a
   ghost atom. This increases the amount of ghost-atom data that must be kept
   synchronized.

Case Study: PairEAMOMP
""""""""""""""""""""""

Many-body pair styles like EAM are more complex because they require multiple
passes (e.g., one pass to accumulate density, one to compute forces).
In `PairEAMOMP`, this is handled by having the `compute()` function
orchestrate multiple parallel loops, each calling its own `eval()` variant:

.. code-block:: c++

   // First loop: accumulate density
   #pragma omp parallel
   {
     eval_rho<...>(...);
   }

   // Communication step:
   comm->forward_comm(this);

   // Second loop: compute forces
   #pragma omp parallel
   {
     eval_force<...>(...);
   }

Case Study: PairSWOMP and PairTersoffOMP
""""""""""""""""""""""""""""""""""""""""

In `PairSWOMP` and `PairTersoffOMP`, the acceleration is handled by parallelizing
the nested loops while maintaining thread-safety.

- **Data Locality**: Because these styles require a full neighbor list, the
  `PairTersoffOMP::compute()` function must carefully partition the neighbor
  list among threads, ensuring that every thread processes a unique set of
  interactions without redundant calculations.
- **Tallying**: Just as in pairwise styles, `ev_tally_thr()` is used to
  accumulate energy/virial. However, since the many-body calculation is more
  involved, the tallying is often called multiple times or with different
  arguments to properly decompose the contributions of the 2-body and 3-body terms.
- **Kernel Specialization**: Similar to the pairwise case, the `eval()` kernels
  for many-body styles are heavily templated to handle various scenarios (e.g.,
  different Newton's third law configurations).

When porting these styles, the primary challenge is to ensure that the
many-body logic remains correct within a multi-threaded environment. Developers
must carefully synchronize data access, especially for shared neighbor list
information and thread-local atom data, by utilizing the `ThrData` structure.

Accelerating bonded force styles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accelerating bonded force styles follows the same philosophy as pair styles.
A `BondLJCutOMP` class would inherit from `BondLJCut` and `ThrOMP`.  The
same principle applies to angle, dihedral, and improper styles.
The `compute()` function again loops over the bonds using OpenMP pragmas,
using the same thread-local tallying mechanisms as pair styles.

.. code-block:: c++

   void BondLJCutOMP::compute(int eflag, int vflag)
   {
     // ...
     #pragma omp parallel
     {
        // loop over bond list in chunks
        // use ev_tally_thr for tallies
     }
     // ...
   }

The reduction logic and templated kernels remain nearly identical to
those used in pair styles, making it relatively straightforward to
port any bonded interaction style to OpenMP.

Implementation file
"""""""""""""""""""

The core of the OPENMP implementation resides in the `.cpp` file (e.g.,
`pair_lj_cut_omp.cpp`).

The `compute()` function is overridden to trigger the parallel
execution and handle thread synchronization:

.. code-block:: c++

   void PairLJCutOMP::compute(int eflag, int vflag)
   {
     // ... (initialization of evflag, vflag, etc.)

     const int nthreads = comm->nthreads;
     const int inum = list->inum;

     #if defined(_OPENMP)
     #pragma omp parallel default(none) shared(eflag,vflag)
     #endif
     {
       // ... (setup thread-local data)

       if (newton_pair) eval<1,1,1>(...);
       else eval<1,1,0>(...);
     }
     // ... (reduce per-thread forces/energies/virials)
   }

1. The `compute()` function initializes the OpenMP parallel region
   using `#pragma omp parallel`.
2. The `eval()` method (a template function) encapsulates the actual
   force-calculation loop. It uses template parameters (`EVFLAG`, `EFLAG`,
   `NEWTON_PAIR`) to optimize the loop based on whether energy/virial
   calculations are enabled and if Newton's third law is used,
   avoiding runtime conditionals inside the hot loop.
3. The `eval()` method uses `ThrData` to accumulate force and energy
   contributions safely within each thread.
4. After the parallel loop, the results from each thread's local buffers
   are reduced into the main force arrays, handling thread-safety during
   the reduction process.

By templating the `eval()` method, we can compile multiple variants of the
kernel at build time, significantly improving performance by removing branches
inside the main force loop.

Refer to existing files in `src/OPENMP/` for comprehensive implementation
details and complex scenarios involving many-body interactions or
communication.
