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
FixStyle(controller,FixController);
// clang-format on
#else

#ifndef LMP_FIX_CONTROLLER_H
#define LMP_FIX_CONTROLLER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixController : public Fix {
 public:
  FixController(class LAMMPS *, int, char **);
  ~FixController() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  void reset_dt() override;
  double compute_vector(int) override;

 private:
  double kp, ki, kd, alpha, tau;
  double setpoint;
  int pvwhich, pvindex;
  char *pvID, *cvID;
  int firsttime;

  double control, err, olderr, deltaerr, sumerr;

  class Compute *pcompute;
  class Fix *pfix;
  int pvar, cvar;
};

}    // namespace LAMMPS_NS

#endif
#endif
