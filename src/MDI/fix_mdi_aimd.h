/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdi/aimd,FixMDIAimd);
// clang-format on
#else

#ifndef LMP_FIX_MDI_AIMD_H
#define LMP_FIX_MDI_AIMD_H

#include "fix.h"
#include <mdi.h>

namespace LAMMPS_NS {

class FixMDIAimd : public Fix {
 public:
  FixMDIAimd(class LAMMPS *, int, char **);
  ~FixMDIAimd();
  int setmask();

  void setup(int);
  void setup_pre_reverse(int, int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int eflag_caller;
  double engine_energy;

  MDI_Comm engine;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
