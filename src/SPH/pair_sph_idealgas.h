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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sph/idealgas,PairSPHIdealGas);
// clang-format on
#else

#ifndef LMP_PAIR_IDEALGAS_H
#define LMP_PAIR_IDEALGAS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHIdealGas : public Pair {
 public:
  PairSPHIdealGas(class LAMMPS *);
  ~PairSPHIdealGas() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double **cut, **viscosity;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
