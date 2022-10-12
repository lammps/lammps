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
FixStyle(oneway,FixOneWay);
// clang-format on
#else

#ifndef LMP_FIX_ONEWAY_H
#define LMP_FIX_ONEWAY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOneWay : public Fix {
 public:
  FixOneWay(class LAMMPS *, int, char **);
  ~FixOneWay() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

 protected:
  int direction;
  class Region *region;
  char *idregion;
};

}    // namespace LAMMPS_NS

#endif
#endif
