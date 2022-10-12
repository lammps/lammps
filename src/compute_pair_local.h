/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pair/local,ComputePairLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_PAIR_LOCAL_H
#define LMP_COMPUTE_PAIR_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePairLocal : public Compute {
 public:
  ComputePairLocal(class LAMMPS *, int, char **);
  ~ComputePairLocal() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_local() override;
  double memory_usage() override;

 private:
  int nvalues, ncount, cutstyle;

  int *pstyle;    // style of each requested output
  int *pindex;    // for pI, index of the output (0 to M-1)
  int singleflag;

  int nmax;
  double *vlocal;
  double **alocal;

  class NeighList *list;

  int compute_pairs(int);
  void reallocate(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
