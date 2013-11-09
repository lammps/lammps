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

#ifdef FIX_CLASS

FixStyle(lb/pc,FixLbPC)

#else

#ifndef LMP_FIX_LB_PC_H
#define LMP_FIX_LB_PC_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixLbPC : public Fix {
 public:
  FixLbPC(class LAMMPS *, int, char **);
  ~FixLbPC();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  //  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);


 private:
  double dtv,dtf;
  int me;
  double *Gamma_MD;
  double expminusdttimesgamma;
  double DMDcoeff;

  double **force_old;
  double **up;
  double **up_old;

  void compute_up(void);
  class FixLbFluid *fix_lb_fluid; 
};

}

#endif
#endif
