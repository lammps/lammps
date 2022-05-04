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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/charmm/coul/msm,PairLJCharmmCoulMSM);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_MSM_H
#define LMP_PAIR_LJ_CHARMM_COUL_MSM_H

#include "pair_lj_charmm_coul_long.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulMSM : public PairLJCharmmCoulLong {
 public:
  PairLJCharmmCoulMSM(class LAMMPS *);
  ~PairLJCharmmCoulMSM() override;
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void compute_outer(int, int) override;
  void *extract(const char *, int &) override;

 protected:
  int nmax;
  double **ftmp;
};

}    // namespace LAMMPS_NS

#endif
#endif
