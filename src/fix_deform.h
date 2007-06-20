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

#ifndef FIX_DEFORM_H
#define FIX_DEFORM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeform : public Fix {
 public:
  int remapflag;

  FixDeform(class LAMMPS *, int, char **);
  ~FixDeform();
  int setmask();
  void init();
  void pre_exchange();
  void end_of_step();

 private:
  double xlo_start,xhi_start,ylo_start,yhi_start,zlo_start,zhi_start;
  double xlo_stop,xhi_stop,ylo_stop,yhi_stop,zlo_stop,zhi_stop;
  double xy_start,xz_start,yz_start;
  double xy_stop,xz_stop,yz_stop;
  double xscale,yscale,zscale;
  int triclinic,scaleflag,flip;
  double *h_rate,*h_ratelo;

  struct Set {
    int style,substyle;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale,vel,rate;
    double lo_start,hi_start;
    double lo_stop,hi_stop;
    double lo_target,hi_target;
    double tilt_start,tilt_stop,tilt_target;
    double vol_start,tilt_flip;
    int fixed,dynamic1,dynamic2;
  };
  Set *set;

  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes

  void options(int, char **);
};

}

#endif
