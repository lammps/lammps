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
FixStyle(SRPREACT,FixSRPREACT);
// clang-format on
#else

#ifndef LMP_FIX_SRP_REACT_H
#define LMP_FIX_SRP_REACT_H

#include "fix_srp.h"

namespace LAMMPS_NS {

class FixSRPREACT : public FixSRP {
 public:
  FixSRPREACT(class LAMMPS *, int, char **);
  ~FixSRPREACT() override;
  int setmask() override;
  void init() override;
  void post_neighbor() override;
  int modify_param(int, char **) override;

 private:
  class FixBondBreak *f_bb;
  char *idbreak;
  class FixBondCreate *f_bc;
  char *idcreate;
};
}    // namespace LAMMPS_NS
#endif
#endif
