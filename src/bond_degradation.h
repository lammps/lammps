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

#ifndef LMP_BOND_DEGRADATION_H
#define LMP_BOND_DEGRADATION_H

#include "lmptype.h"

namespace LAMMPS_NS {

// Internal interface returned by Fix::extract("bond_degradation", dim).
// Bond styles can query damage without depending on an optional fix package.
class BondDegradation {
 public:
  virtual ~BondDegradation() noexcept(false) = default;
  virtual double bond_scale(tagint, tagint, int) const = 0;
};

}    // namespace LAMMPS_NS

#endif
