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
FixStyle(dpd/energy,FixDPDenergy);
// clang-format on
#else

#ifndef LMP_FIX_DPDE_H
#define LMP_FIX_DPDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDPDenergy : public Fix {
 public:
  FixDPDenergy(class LAMMPS *, int, char **);
  virtual ~FixDPDenergy() {}
  int setmask();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  class PairDPDfdtEnergy *pairDPDE;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
