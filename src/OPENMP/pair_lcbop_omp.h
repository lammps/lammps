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
PairStyle(lcbop/omp,PairLCBOPOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LCBOP_OMP_H
#define LMP_PAIR_LCBOP_OMP_H

#include "pair_lcbop.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLCBOPOMP : public PairLCBOP, public ThrOMP {

 public:
  PairLCBOPOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 protected:
  void SR_neigh_thr();
  void FSR_thr(int ifrom, int ito, int eflag, ThrData *const thr);
  void FLR_thr(int ifrom, int ito, int eflag, ThrData *const thr);

  void FNij_thr(int i, int j, double factor, double **f, ThrData *const thr);
  void FMij_thr(int i, int j, double factor, double **f, ThrData *const thr);
  double bondorder_thr(int i, int j, double *rij, double rijmag, double VA,
                       double **f, ThrData *const thr);
  double b_thr(int i, int j, double *rij, double rijmag, double VA,
               double **f, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
