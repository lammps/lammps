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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/cut/dielectric,PairLJCutCoulCutDielectric);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_DIELECTRIC_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_DIELECTRIC_H

#include "pair_lj_cut_coul_cut.h"

namespace LAMMPS_NS {

class PairLJCutCoulCutDielectric : public PairLJCutCoulCut {
 public:
  PairLJCutCoulCutDielectric(class LAMMPS *);
  ~PairLJCutCoulCutDielectric() override;
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void init_style() override;

  double **efield;
  double *epot;

 protected:
  class AtomVecDielectric *avec;
  int nmax;
};

}    // namespace LAMMPS_NS

#endif
#endif
