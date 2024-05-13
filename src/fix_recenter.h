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
FixStyle(recenter,FixRecenter);
// clang-format on
#else

#ifndef LMP_FIX_RECENTER_H
#define LMP_FIX_RECENTER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRecenter : public Fix {
 public:
  FixRecenter(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void initial_integrate_respa(int, int, int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  int group2bit, scaleflag;
  int xflag, yflag, zflag;
  int xinitflag, yinitflag, zinitflag;
  int nlevels_respa;
  double xcom, ycom, zcom, xinit, yinit, zinit, masstotal, distance, shift[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
