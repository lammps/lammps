/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(set,FixSet);
// clang-format on
#else

#ifndef LMP_FIX_SET_H
#define LMP_FIX_SET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSet : public Fix {
 public:
  FixSet(class LAMMPS *, int, char **);
  ~FixSet() override;
  int setmask() override;
  void end_of_step() override;

 private:
  int reneighbor;

  class Set *set;
};

}    // namespace LAMMPS_NS

#endif
#endif
