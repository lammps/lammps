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

#ifdef BOND_CLASS
// clang-format off
BondStyle(quartic/unbreakable,BondQuarticUnbreakable);
// clang-format on
#else

#include "bond_quartic_breakable.h"

namespace LAMMPS_NS {

class BondQuarticUnbreakable : public BondQuarticBreakable {
 public:
  BondQuarticUnbreakable(class LAMMPS *);
  ~BondQuarticUnbreakable() override;

};

}    // namespace LAMMPS_NS

#endif
