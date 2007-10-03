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
  void end_of_step();

 private:
  int me;
  int nrepeat,nfreq,irepeat,nvalid;
  int dim,originflag,which,normflag;
  double origin,delta;
  char *id_compute,*id_fix;
  FILE *fp;

  int nlayers,nvalues,maxlayer,scaleflag,size_peratom;
  double xscale,yscale,zscale;
  double layer_volume;
  double *coord;
  double *count_one,*count_many,*count_total;
  double **values_one,**values_many,**values_total;
  double offset,invdelta; 
  int ncompute;
  class Compute **compute;
  class Fix *fix;
};

}

#endif
