/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_SPRING_H
#define FIX_SPRING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpring : public Fix {
 public:
  FixSpring(class LAMMPS *, int, char **);
  ~FixSpring();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  double xc,yc,zc,r0;
  double k_spring;
  int xflag,yflag,zflag;
  int styleflag;
  char *group2;
  int igroup2,group2bit;
  double masstotal,masstotal2;
  int nlevels_respa;
  double espring,ftotal[4];
  int force_flag;

  void spring_tether();
  void spring_couple();
};

}

#endif
