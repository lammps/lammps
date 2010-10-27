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

FixStyle(ave/spatial,FixAveSpatial)

#else

#ifndef LMP_FIX_AVE_SPATIAL_H
#define LMP_FIX_AVE_SPATIAL_H

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
  double compute_array(int,int);
  double memory_usage();

 private:
  int me,nvalues;
  int nrepeat,nfreq,nvalid,irepeat;
  int ndim,normflag,regionflag,iregion;
  char *tstring,*sstring,*idregion;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  class Region *region;

  int ave,nwindow,scaleflag;
  int norm,iwindow,window_limit;
  double xscale,yscale,zscale;
  double bin_volume;

  int dim[3],originflag[3],nlayers[3];
  double origin[3],delta[3];
  double offset[3],invdelta[3]; 

  int maxvar;
  double *varatom;

  int maxatom;
  int *bin;

  int nbins,maxbin;
  double **coord;
  double *count_one,*count_many,*count_sum;
  double **values_one,**values_many,**values_sum;
  double *count_total,**count_list;
  double **values_total,***values_list;

  void setup_bins();
  void atom2bin1d();
  void atom2bin2d();
  void atom2bin3d();
};

}

#endif
#endif
