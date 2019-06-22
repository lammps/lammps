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

FixStyle(cac/nve,FixNVECAC)

#else

#ifndef LMP_FIX_NVE_CAC_H
#define LMP_FIX_NVE_CAC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVECAC : public Fix {
 public:
  FixNVECAC(class LAMMPS *, int, char **);
  virtual ~FixNVECAC() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();

 protected:
  double dtv,dtf;
  double *step_respa;
  int mass_require;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: CAC fix styles require a CAC atom style

Self-explanatory.

*/
