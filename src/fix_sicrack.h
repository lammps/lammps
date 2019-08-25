/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Mode I, II or III semi-infinite crack, with T-stress
/*------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sicrack,FixSICrack)

#else

#ifndef LMP_FIX_SICRACK_H
#define LMP_FIX_SICRACK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSICrack : public Fix {
 public:
  FixSICrack(class LAMMPS *, int, char **);
  ~FixSICrack();  
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:
  double xtip,ytip,mu,vu,K1,dK1,K2,dK2,K3,dK3,T,dT;
  double PI,Xk;

  double xtemp,ytemp,rdist,theta,theta2,cost,sint,ux,uy,uz,uTx,uTy;
  double **xoriginal;         // original coords of atoms
};

}

#endif
#endif
