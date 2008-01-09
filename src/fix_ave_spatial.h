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

#ifndef FIX_AVE_SPATIAL_H
#define FIX_AVE_SPATIAL_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveSpatial : public Fix {
 public:
  FixAveSpatial(class LAMMPS *, int, char **);
  ~FixAveSpatial();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_vector(int);
  double memory_usage();

 private:
  int me,nvalues;
  int nrepeat,nfreq,nvalid,irepeat;
  int dim,originflag,normflag;
  double origin,delta;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;

  int nlayers,ave,nwindow;
  int maxlayer,scaleflag;
  double xscale,yscale,zscale;
  double layer_volume;
  double *coord;
  double *count_one,*count_many,*count_sum;
  double **values_one,**values_many,**values_sum;
  double offset,invdelta; 

  int maxatomvar;
  double *varatom;

  int maxatomlayer;
  int *layer;

  int norm,iwindow,window_limit;
  double *count_total;
  double **count_list;
  double **values_total;
  double ***values_list;
};

}

#endif
