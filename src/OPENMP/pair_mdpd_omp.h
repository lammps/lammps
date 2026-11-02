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
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mdpd/omp,PairMDPDOMP);
// clang-format on
#else

#ifndef LMP_PAIR_MDPD_OMP_H
#define LMP_PAIR_MDPD_OMP_H

#include "pair_mdpd.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairMDPDOMP : public PairMDPD, public ThrOMP {

 public:
  PairMDPDOMP(class LAMMPS *);
  ~PairMDPDOMP() override;

  void settings(int, char **) override;
  void compute(int, int) override;
  double memory_usage() override;

 protected:
  class RanMars **random_thr;
  int nthreads;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
