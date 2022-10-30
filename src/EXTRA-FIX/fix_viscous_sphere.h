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
FixStyle(viscous/sphere,FixViscousSphere);
// clang-format on
#else

#ifndef LMP_FIX_VISCOUS_SPHERE_H
#define LMP_FIX_VISCOUS_SPHERE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixViscousSphere : public Fix {
 public:
  FixViscousSphere(class LAMMPS *, int, char **);
  virtual ~FixViscousSphere();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double memory_usage();

 protected:
  double gamma, *scalegamma, *scaleval;
  const char *scalevarid;
  int scalestyle, scalevar;
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
