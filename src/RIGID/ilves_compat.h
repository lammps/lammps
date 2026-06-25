/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   ILVES constraint solver: compatibility shim.

   The ILVES solver core (ilves_schur_solver, ilves_molecule, ilves) is ported
   near-verbatim from the GROMACS 2021 reference implementation (LGPL-2.1).  To
   keep that code changed as little as possible, this header supplies the few
   GROMACS spellings the solver relies on:

     * "real"            -> double (ILVES is solved in double precision)
     * "AlignedAllocator" -> std::allocator (the scalar port drops SIMD, so the
                             special alignment is unnecessary)
     * the DIM / XX / YY / ZZ dimension constants and square()
     * OpenMP entry points: the real ones when built with OpenMP, otherwise
       no-op stubs so the per-partition OpenMP code the solver carries compiles
       and runs correctly on a single thread.

   See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_COMPAT_H
#define LMP_ILVES_COMPAT_H

#include <memory>

#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_lock_t;
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_max_threads() { return 1; }
inline void omp_init_lock(omp_lock_t *) {}
inline void omp_destroy_lock(omp_lock_t *) {}
inline void omp_set_lock(omp_lock_t *) {}
inline void omp_unset_lock(omp_lock_t *) {}
#endif

namespace LAMMPS_NS {
namespace ILVES {

// ILVES is solved in double precision in LAMMPS.
using real = double;

// GROMACS used an aligned allocator for its SIMD kernels.  The scalar LAMMPS
// port does not use SIMD, so the default allocator is sufficient.  Keeping the
// name lets the ported containers stay textually identical to upstream.
template <class T> using AlignedAllocator = std::allocator<T>;

// dimension constants, matching GROMACS DIM / XX / YY / ZZ
enum { XX = 0, YY = 1, ZZ = 2, DIM = 3 };

template <class T> inline T square(T x) { return x * x; }

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif
