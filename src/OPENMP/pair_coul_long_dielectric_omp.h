/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
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
PairStyle(coul/long/dielectric/omp,PairCoulLongDielectricOMP);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_LONG_DIELECTRIC_OMP_H
#define LMP_PAIR_COUL_LONG_DIELECTRIC_OMP_H

#include "pair_coul_long_dielectric.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairCoulLongDielectricOMP : public PairCoulLongDielectric, public ThrOMP {

 public:
  PairCoulLongDielectricOMP(class LAMMPS *);
  ~PairCoulLongDielectricOMP() override = default;
  void compute(int, int) override;

 protected:
  template <int EVFLAG, int EFLAG> void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
