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

FixStyle(nve/eff,FixNVEEff)

#else

#ifndef LMP_FIX_NVE_EFF_H
#define LMP_FIX_NVE_EFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVEEff : public Fix {
 public:
  FixNVEEff(class LAMMPS *, int, char **);
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  void reset_dt();

 protected:
  double dtv,dtf;
  double *step_respa;
  int mass_require;
};

}

#endif
#endif
