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

#ifdef COMPUTE_CLASS

ComputeStyle(temp/profile,ComputeTempProfile)

#else

#ifndef LMP_COMPUTE_TEMP_PROFILE_H
#define LMP_COMPUTE_TEMP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempProfile : public Compute {
 public:
  ComputeTempProfile(class LAMMPS *, int, char **);
  ~ComputeTempProfile();
  void init();
  void setup();
  double compute_scalar();
  void compute_vector();
  void compute_array();

  void remove_bias(int, double *);
  void remove_bias_all();
  void restore_bias(int, double *);
  void restore_bias_all();
  double memory_usage();

 private:
  int xflag,yflag,zflag,ncount,outflag;
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
  double *tbin,*tbinall;

  void dof_compute();
  void bin_average();
  void bin_setup();
  void bin_assign();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute temp/profile cannot use vz for 2d systemx

Self-explanatory.

E: Compute temp/profile cannot bin z for 2d systems

Self-explanatory.

*/
