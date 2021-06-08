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
FixStyle(sph/stationary,FixSPHStationary);
// clang-format on
#else

#ifndef LMP_FIX_SPH_STATIONARY_H
#define LMP_FIX_SPH_STATIONARY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSPHStationary : public Fix {
 public:
  FixSPHStationary(class LAMMPS *, int, char **);
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void reset_dt();

 private:
  class NeighList *list;

 protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;

  class Pair *pair;
};

}    // namespace LAMMPS_NS

#endif
#endif
