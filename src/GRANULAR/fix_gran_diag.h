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

#ifndef FIX_GRAN_DIAG_H
#define FIX_GRAN_DIAG_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixGranDiag : public Fix {
 public:
  FixGranDiag(class LAMMPS *, int, char **);
  ~FixGranDiag();
  int setmask();
  void init();
  void setup();
  void end_of_step();

 private:
  int me,first,pairstyle,nlayers,maxlayers,dim;
  FILE *fpden,*fpvel,*fpstr;
  double step,stepinv,boxlo,PI;
  double dt,xkk,xkkt,xmu,gamman_dl,gammas_dl;
  int freeze_group_bit;

  int *numdens;
  double *dendens;
  double *velx,*vely,*velz;
  double *velxx,*velyy,*velzz,*velxy,*velxz,*velyz;
  double *sigxx,*sigyy,*sigzz,*sigxy,*sigxz,*sigyz;
  double *velx11,*vely11,*velz11;
  double *velxx11,*velyy11,*velzz11,*velxy11,*velxz11,*velyz11;
  double *velfxx,*velfyy,*velfzz,*velfxy,*velfxz,*velfyz;
  double *sigx2,*sigy2,*sigz2,*sigxy2,*sigxz2,*sigyz2;

  void allocate();
  void deallocate();
  void stress_no_history();
  void stress_history();
  void stress_hertzian();
};

}

#endif
