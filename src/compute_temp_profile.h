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

#ifndef COMPUTE_TEMP_PROFILE_H
#define COMPUTE_TEMP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempProfile : public Compute {
 public:
  ComputeTempProfile(class LAMMPS *, int, char **);
  ~ComputeTempProfile();
  void init();
  double compute_scalar();
  void compute_vector();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 private:
  int xflag,yflag,zflag,ncount;
  int nbinx,nbiny,nbinz,nbins;
  int ivx,ivy,ivz;
  int fix_dof;
  double tfactor;

  int box_change,triclinic;
  int *periodicity;
  double *boxlo,*boxhi,*prd;
  double invdelta[3];

  int maxatom;
  int *bin;
  double **vbin,**binave;
  
  void dof_compute();
  void bin_average();
  void bin_setup();
  void bin_assign();
};

}

#endif
