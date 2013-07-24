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

FixStyle(spring/pull,FixSpringPull)

#else

#ifndef LMP_FIX_SPRING_PULL_H
#define LMP_FIX_SPRING_PULL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringPull : public Fix {
 public:
  FixSpringPull(class LAMMPS *, int, char **);
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void min_setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  virtual void min_post_force(int);
  virtual double compute_scalar();
  virtual double compute_vector(int);

 private:
  double xc,yc,zc,r0;
  double xv,yv,zv;
  double k_spring;
  int xflag,yflag,zflag;
  double masstotal;
  int nlevels_respa;
  double espring,ftotal[8];
};

}

#endif
#endif
