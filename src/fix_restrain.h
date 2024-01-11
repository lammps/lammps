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
FixStyle(restrain,FixRestrain);
// clang-format on
#else

#ifndef LMP_FIX_RESTRAIN_H
#define LMP_FIX_RESTRAIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRestrain : public Fix {
 public:
  FixRestrain(class LAMMPS *, int, char **);
  ~FixRestrain() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  int ilevel_respa;
  int nrestrain, maxrestrain;
  int *rstyle;
  int *mult;
  tagint **ids;
  double *kstart, *kstop, *deqstart, *deqstop, *target;
  double *cos_target, *sin_target;
  double energy, energy_all;
  double ebond, ebond_all;
  double elbound, elbound_all;
  double eangle, eangle_all;
  double edihed, edihed_all;

  void restrain_bond(int);
  void restrain_lbound(int);
  void restrain_angle(int);
  void restrain_dihedral(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
