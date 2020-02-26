/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(bd/euler,FixBDEuler)

#else

#ifndef LMP_FIX_BD_EULER_H
#define LMP_FIX_BD_EULER_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixBDEuler : public FixNVE {
 public:
  FixBDEuler(class LAMMPS *, int, char **);
  virtual ~FixBDEuler();
  void init();
  void initial_integrate(int);
  //void final_integrate();

 private:
  double dt, sqrtdt;
 protected:
  int seed;
  double t_start,t_stop,t_target,tsqrt;
  double diff;
  double gamma1,gamma2, gamma3, gamma4;
  double cosda, sinda, da, dar;
  class RanMars *random;
  void compute_target();
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Compute bd/euler requires atom style ellipsoid

Self-explanatory.

E: Fix bd/euler requires extended particles

This fix can only be used for particles with a shape setting.

*/
