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

FixStyle(ave/time,FixAveTime)

#else

#ifndef LMP_FIX_AVE_TIME_H
#define LMP_FIX_AVE_TIME_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class LAMMPS *, int, char **);
  ~FixAveTime();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  int me,nvalues;
  int nrepeat,nfreq,nvalid,irepeat;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  int nrows;

  int ave,nwindow,nsum,startstep,mode;
  int *offcol;
  char *title1,*title2,*title3;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;
  double *column;
  double **array;
  double **array_total;
  double ***array_list;

  void invoke_scalar(int);
  void invoke_vector(int);
  void options(int, char **);
};

}

#endif
#endif
