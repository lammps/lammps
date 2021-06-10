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
FixStyle(mvv/edpd,FixMvvEDPD);
// clang-format on
#else

#ifndef LMP_FIX_MVV_EDPD_H
#define LMP_FIX_MVV_EDPD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMvvEDPD : public Fix {
 public:
  FixMvvEDPD(class LAMMPS *, int, char **);
  virtual ~FixMvvEDPD() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void reset_dt();

 protected:
  double dtv, dtf;
  double verlet;
};

}    // namespace LAMMPS_NS

#endif
#endif
