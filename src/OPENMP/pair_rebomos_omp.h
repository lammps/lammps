/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(rebomos/omp,PairREBOMoSOMP);
// clang-format on
#else

#ifndef LMP_PAIR_REBOMOS_OMP_H
#define LMP_PAIR_REBOMOS_OMP_H

#include "pair_rebomos.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairREBOMoSOMP : public PairREBOMoS, public ThrOMP {
 public:
  PairREBOMoSOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 protected:
  void FREBO_thr(int ifrom, int ito, int eflag, ThrData *const thr);
  void FLJ_thr(int ifrom, int ito, int eflag, ThrData *const thr);

  void REBO_neigh_thr();

  double bondorder_thr(int, int, double *, double, double, ThrData *const thr);
};
}    // namespace LAMMPS_NS

#endif
#endif
