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

#ifndef FIX_VOL_RESCALE_H
#define FIX_VOL_RESCALE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixVolRescale : public Fix {
 public:
  FixVolRescale(class LAMMPS *, int, char **);
  ~FixVolRescale();
  int setmask();
  void init();
  void end_of_step();

 private:
  int xflag,yflag,zflag;
  double xlo_start,xlo_stop,xhi_start,xhi_stop;
  double ylo_start,ylo_stop,yhi_start,yhi_stop;
  double zlo_start,zlo_stop,zhi_start,zhi_stop;
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes
};

}

#endif
