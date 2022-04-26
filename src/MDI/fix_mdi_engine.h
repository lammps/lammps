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
FixStyle(MDI/ENGINE, FixMDIEngine);
// clang-format on
#else

#ifndef LMP_FIX_MDI_ENGINE_H
#define LMP_FIX_MDI_ENGINE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMDIEngine : public Fix {
 public:
  class MDIEngine *mdi_engine;

  FixMDIEngine(class LAMMPS *, int, char **);
  ~FixMDIEngine() {}
  int setmask();
  void init() {}
  void setup(int);
  void post_integrate();
  void min_pre_force(int);
  void post_force(int);
  void min_post_force(int);
  void end_of_step();
};

}    // namespace LAMMPS_NS

#endif
#endif
