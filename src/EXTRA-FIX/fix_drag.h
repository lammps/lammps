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
FixStyle(drag,FixDrag);
// clang-format on
#else

#ifndef LMP_FIX_DRAG_H
#define LMP_FIX_DRAG_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDrag : public Fix {
 public:
  FixDrag(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_vector(int) override;

 private:
  double xc, yc, zc;
  double f_mag;
  int xflag, yflag, zflag;
  double delta;
  int ilevel_respa;
  double ftotal[3], ftotal_all[3];
  int force_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
