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

FixStyle(gle,FixGLE)

#else

#ifndef LMP_FIX_GLE_H
#define LMP_FIX_GLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGLE : public Fix {
 public:
  FixGLE(class LAMMPS *, int, char **);
  virtual ~FixGLE();
  int setmask();
  void init();
  void setup(int);
  void gle_integrate();
  void initial_integrate_respa(int vflag, int ilevel, int iloop);
  void final_integrate_respa(int ilevel, int iloop);
  void initial_integrate(int vflag);
  void final_integrate();
  double compute_scalar();
  virtual void *extract(const char *, int &);

 protected:
  int ns;
  double *A, *C, *S, *T;
  double *gle_tmp, *gle_rnd; int gle_buff;
  double temp, tau, dtv, dtf;
  int dogle;
  int dorattle;
  class RanMars *random;
  double *sqrt_m;
  double *step_respa;
  double energy;
  int nlevels_respa;

  double **vaux;
};

}

#endif
#endif
