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

#ifndef FIX_LANGEVIN_H
#define FIX_LANGEVIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevin : public Fix {
 public:
  FixLangevin(class LAMMPS *, int, char **);
  ~FixLangevin();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void reset_target(double);
  void reset_dt();

 private:
  double t_start,t_stop,t_period;
  int flagx,flagy,flagz,iregion;
  double *gfactor1,*gfactor2,*ratio;

  int nlevels_respa;
  class RanMars *random;
};

}

#endif
