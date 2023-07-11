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
PairStyle(lepton/sphere,PairLeptonSphere);
// clang-format on
#else

#ifndef LMP_PAIR_LEPTON_SPHERE_H
#define LMP_PAIR_LEPTON_SPHERE_H

#include "pair_lepton.h"

namespace LAMMPS_NS {

class PairLeptonSphere : public PairLepton {
 public:
  PairLeptonSphere(class LAMMPS *_lmp) : PairLepton(_lmp){};

  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
};
}    // namespace LAMMPS_NS
#endif
#endif
