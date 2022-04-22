OpenMP Parallelism
^^^^^^^^^^^^^^^^^^

The styles in the INTEL, KOKKOS, and OPENMP package offer to use OpenMP
thread parallelism to predominantly distribute loops over local data
and thus follow an orthogonal parallelization strategy to the
decomposition into spatial domains used by the :doc:`MPI partitioning
<Developer_par_part>`.  For clarity, this section discusses only the
implementation in the OPENMP package as it is the simplest. The INTEL
and KOKKOS package offer additional options and are more complex since
they support more features and different hardware like co-processors
or GPUs.

One of the key decisions when implementing the OPENMP package was to
keep the changes to the source code small, so that it would be easier to
maintain the code and keep it in sync with the non-threaded standard
implementation.  this is achieved by a) making the OPENMP version a
derived class from the regular version (e.g. ``PairLJCutOMP`` from
``PairLJCut``) and overriding only methods that are multi-threaded or
need to be modified to support multi-threading (similar to what was done
in the OPT package), b) keeping the structure in the modified code very
similar so that side-by-side comparisons are still useful, and c)
offloading additional functionality and multi-thread support functions
into three separate classes ``ThrOMP``, ``ThrData``, and ``FixOMP``.
``ThrOMP`` provides additional, multi-thread aware functionality not
available in the corresponding base class (e.g. ``Pair`` for
``PairLJCutOMP``) like multi-thread aware variants of the "tally"
functions. Those functions are made available through multiple
inheritance so those new functions have to have unique names to avoid
ambiguities; typically ``_thr`` is appended to the name of the function.
``ThrData`` is a classes that manages per-thread data structures.
It is used instead of extending the corresponding storage to per-thread
arrays to avoid slowdowns due to "false sharing" when multiple threads
update adjacent elements in an array and thus force the CPU cache lines
to be reset and re-fetched.  ``FixOMP`` finally manages the "multi-thread
state" like settings and access to per-thread storage, it is activated
by the :doc:`package omp <package>` command.

Avoiding data races
"""""""""""""""""""

A key problem when implementing thread parallelism in an MD code is
to avoid data races when updating accumulated properties like forces,
energies, and stresses.  When interactions are computed, they always
involve multiple atoms and thus there are race conditions when multiple
threads want to update per-atom data of the same atoms.  Five possible
strategies have been considered to avoid this:

1) restructure the code so that there is no overlapping access possible
   when computing in parallel, e.g. by breaking lists into multiple
   parts and synchronizing threads in between.
2) have each thread be "responsible" for a specific group of atoms and
   compute these interactions multiple times, once on each thread that
   is responsible for a given atom and then have each thread only update
   the properties of this atom.
3) use mutexes around functions and regions of code where the data race
   could happen
4) use atomic operations when updating per-atom properties
5) use replicated per-thread data structures to accumulate data without
   conflicts and then use a reduction to combine those results into the
   data structures used by the regular style.

Option 5 was chosen for the OPENMP package because it would retain the
performance for the case of 1 thread and the code would be more
maintainable.  Option 1 would require extensive code changes,
particularly to the neighbor list code; options 2 would have incurred a
2x or more performance penalty for the serial case; option 3 causes
significant overhead and would enforce serialization of operations in
inner loops and thus defeat the purpose of multi-threading; option 4
slows down the serial case although not quite as bad as option 2.  The
downside of option 5 is that the overhead of the reduction operations
grows with the number of threads used, so there would be a crossover
point where options 2 or 4 would result in faster executing.  That is
why option 2 for example is used in the GPU package because a GPU is a
processor with a massive number of threads.  However, since the MPI
parallelization is generally more effective for typical MD systems, the
expectation is that thread parallelism is only used for a smaller number
of threads (2-8).  At the time of its implementation, that number was
equivalent to the number of CPU cores per CPU socket on high-end
supercomputers.

Thus arrays like the force array are dimensioned to the number of atoms
times the number of threads when enabling OpenMP support and inside the
compute functions a pointer to a different chunk is obtained by each thread.
Similarly, accumulators like potential energy or virial are kept in
per-thread instances of the ``ThrData`` class and then only reduced and
stored in their global counterparts at the end of the force computation.


Loop scheduling
"""""""""""""""

Multi-thread parallelization is applied by distributing (outer) loops
statically across threads.  Typically this would be the loop over local
atoms *i* when processing *i,j* pairs of atoms from a neighbor list.
The design of the neighbor list code results in atoms having a similar
number of neighbors for homogeneous systems and thus load imbalances
across threads are not common and typically happen for systems where
also the MPI parallelization would be unbalanced, which would typically
have a more pronounced impact on the performance.  This same loop
scheduling scheme can also be applied to the reduction operations on
per-atom data to try and reduce the overhead of the reduction operation.

Neighbor list parallelization
"""""""""""""""""""""""""""""

In addition to the parallelization of force computations, also the
generation of the neighbor lists is parallelized.  As explained
previously, neighbor lists are built by looping over "owned" atoms and
storing the neighbors in "pages".  In the OPENMP variants of the
neighbor list code, each thread operates on a different chunk of "owned"
atoms and allocates and fills its own set of pages with neighbor list
data.  This is achieved by each thread keeping its own instance of the
:cpp:class:`MyPage <LAMMPS_NS::MyPage>` page allocator class.
