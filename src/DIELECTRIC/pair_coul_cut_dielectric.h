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
PairStyle(coul/cut/dielectric,PairCoulCutDielectric);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_CUT_DIELECTRIC_H
#define LMP_PAIR_COUL_CUT_DIELECTRIC_H

#include "pair_coul_cut.h"

namespace LAMMPS_NS {

class PairCoulCutDielectric : public PairCoulCut {
 public:
  PairCoulCutDielectric(class LAMMPS *);
  ~PairCoulCutDielectric() override;
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void init_style() override;

  double **efield;

 protected:
  class AtomVecDielectric *avec;
  int nmax;
};

}    // namespace LAMMPS_NS

#endif
#endif
