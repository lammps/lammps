/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_NVE_DIPOLE_H
#define FIX_NVE_DIPOLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVEDipole : public Fix {
 public:
  FixNVEDipole(class LAMMPS *, int, char **);
  ~FixNVEDipole();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int);
  void reset_dt();

 private:
  double dtv,dtf;
  double *step_respa;
  double *inertia;
};

}

#endif
