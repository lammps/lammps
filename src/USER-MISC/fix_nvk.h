/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvk,FixNVK);
// clang-format on
#else

#ifndef LMP_FIX_NVK_H
#define LMP_FIX_NVK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVK : public Fix {
 public:
  FixNVK(class LAMMPS *, int, char **);
  virtual ~FixNVK() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void reset_dt();

 protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;
  double K_target;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nvk only supports group all

Self-explanatory.

E: Fix nvk not yet enabled for RESPA

Self-explanatory.

*/
