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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(orientorder/atom/omp,ComputeOrientOrderAtomOMP);
// clang-format on
#else

#ifndef LMP_COMPUTE_ORIENTORDER_ATOM_OMP_H
#define LMP_COMPUTE_ORIENTORDER_ATOM_OMP_H

#include "compute_orientorder_atom.h"

namespace LAMMPS_NS {

class ComputeOrientOrderAtomOMP : public ComputeOrientOrderAtom {
 public:
  ComputeOrientOrderAtomOMP(class LAMMPS *, int, char **);
  ~ComputeOrientOrderAtomOMP() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nthreads_alloc;    // number of per-thread scratch slots allocated
  int maxneigh_thr;      // per-thread neighbor scratch capacity
  double **distsq_thr;   // [nthreads][maxneigh_thr]
  int **nearest_thr;     // [nthreads][maxneigh_thr]
  double ***rlist_thr;   // [nthreads][maxneigh_thr][3]
  double ***qnm_r_thr;   // [nthreads][nqlist][qmax+1]
  double ***qnm_i_thr;   // [nthreads][nqlist][qmax+1]

  void free_thr_arrays();
};

}    // namespace LAMMPS_NS

#endif
#endif
