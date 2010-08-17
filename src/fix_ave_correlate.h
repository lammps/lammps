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

FixStyle(ave/correlate,FixAveCorrelate)

#else

#ifndef LMP_FIX_AVE_CORRELATE_H
#define LMP_FIX_AVE_CORRELATE_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveCorrelate : public Fix {
 public:
  FixAveCorrelate(class LAMMPS *, int, char **);
  ~FixAveCorrelate();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_array(int,int);

 private:
  int me,nvalues;
  int nrepeat,nfreq,nvalid;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;

  int type,ave,startstep;
  double prefactor;
  char *title1,*title2,*title3;

  int firstindex;      // index in values ring of earliest time sample
  int lastindex;       // index in values ring of latest time sample
  int nsample;         // number of time samples in values ring

  int npair;           // number of correlation pairs to calculate
  int *count;
  double **values,**corr;

  int *save_count;     // saved values at Nfreq for output via compute_array()
  double **save_corr;
    
  void accumulate();
};

}

#endif
#endif
