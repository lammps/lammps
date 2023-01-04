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

#ifdef FIX_CLASS
// clang-format off
FixStyle(tfmc,FixTFMC);
// clang-format on
#else

#ifndef LMP_FIX_TFMC_H
#define LMP_FIX_TFMC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTFMC : public Fix {
 public:
  FixTFMC(class LAMMPS *, int, char **);
  ~FixTFMC() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;

 private:
  double d_max;
  double T_set;
  double mass_min;
  double **xd;
  int mass_require;
  int seed;
  int comflag, rotflag, xflag, yflag, zflag;
  int nmax;
  class RanMars *random_num;
};

}    // namespace LAMMPS_NS

#endif
#endif
