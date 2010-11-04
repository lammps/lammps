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

FixStyle(adapt,FixAdapt)

#else

#ifndef LMP_FIX_ADAPT_H
#define LMP_FIX_ADAPT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdapt : public Fix {
 public:
  int diamflag;        // 1 if atom diameters will vary, for AtomVecGranular

  FixAdapt(class LAMMPS *, int, char **);
  ~FixAdapt();
  int setmask();
  void init();
  void setup_pre_force(int);
  void pre_force(int);
  void post_run();

 private:
  int nadapt,resetflag,scaleflag;
  int anypair;

  struct Adapt {
    int which,ivar;
    char *var;
    char *pstyle,*pparam;
    int ilo,ihi,jlo,jhi;
    int pdim;
    double *scalar,scalar_orig;
    double **array,**array_orig;
    int aparam;
  };

  Adapt *adapt;
  double *kspace_scale;

  void change_settings();
  void restore_settings();
};

}

#endif
#endif
