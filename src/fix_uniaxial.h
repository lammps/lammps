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

#ifndef FIX_UNIAXIAL_H
#define FIX_UNIAXIAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixUniaxial : public Fix {
 public:
  FixUniaxial(class LAMMPS *, int, char **);
  ~FixUniaxial();
  int setmask();
  void init();
  void end_of_step();

 private:
  int dir;
  double lambda_final;
  double alpha0;                   // asymmetry parameter
  
  double lambdai[3];               // initial
  double lambdaf[3];               // final
  double domainloi[3];             // initial box lo/hi
  double domainhii[3];
  double *domainlo[3];             // pointers to make loop over dims possible
  double *domainhi[3];
  double *domainprd[3];
  
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes
};

#endif

}


