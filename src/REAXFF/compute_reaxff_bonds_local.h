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
   Contributing author: Richard Berger (LANL)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(reaxff/bonds/local,ComputeReaxFFBondsLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_REAXFF_BONDS_LOCAL_H
#define LMP_COMPUTE_REAXFF_BONDS_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReaxFFBondsLocal : public Compute {
 public:
  ComputeReaxFFBondsLocal(class LAMMPS *, int, char **);
  ~ComputeReaxFFBondsLocal() override;
  void init() override;
  void compute_local() override;
  double memory_usage() override;

 private:
  int nlocal;
  int nvalues;
  int prev_nvalues;

  double **alocal;
  tagint **neighid;
  double **abo;
  int *numneigh;
  class PairReaxFF *reaxff;

  int FindBond();

  void allocate(int);
  void destroy();
  void reallocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
