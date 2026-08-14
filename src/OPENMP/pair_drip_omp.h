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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(drip/omp,PairDRIPOMP);
// clang-format on
#else

#ifndef LMP_PAIR_DRIP_OMP_H
#define LMP_PAIR_DRIP_OMP_H

#include "pair_drip.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairDRIPOMP : public PairDRIP, public ThrOMP {
 public:
  PairDRIPOMP(class LAMMPS *);
  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int VFLAG_EITHER>
  void eval(int ifrom, int ito, ThrData *const thr);

  double calc_repulsive_thr(int const i, int const j, Param &p, double const rsq,
                            double const *rvec, double const *ni, V3 const *dni_dri,
                            V3 const *dni_drnb1, V3 const *dni_drnb2, V3 const *dni_drnb3,
                            double *const fi, double *const fj, double **f_thr,
                            ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
