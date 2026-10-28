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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/omp, PairPACEOMP);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_OMP_H
#define LMP_PAIR_PACE_OMP_H

#include "pair_pace.h"
#include "thr_omp.h"

#include <vector>

namespace LAMMPS_NS {
struct ACEImplThr;

// inherit from PairPACE and ThrOMP
class PairPACEOMP : public PairPACE, public ThrOMP {
 public:
  PairPACEOMP(class LAMMPS *);
  ~PairPACEOMP() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  double memory_usage() override;

 private:
  ACEImplThr *ace_thr;

  // ACEImpl contains two separate objects that are not thread-safe by default, allocate copies on each thread
  // Helper to create one ACERecursiveEvaluator per thread
  void build_thread_evaluators();
};

}    // namespace LAMMPS_NS

#endif
#endif
