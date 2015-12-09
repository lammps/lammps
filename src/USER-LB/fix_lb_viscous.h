/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(lb/viscous,FixLbViscous)

#else

#ifndef LMP_FIX_LB_VISCOUS_H
#define LMP_FIX_LB_VISCOUS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLbViscous : public Fix {
 public:
  FixLbViscous(class LAMMPS *, int, char **);
  ~FixLbViscous();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 protected:
  int nlevels_respa;

 private:
  class FixLbFluid *fix_lb_fluid;
};

}

#endif
#endif
