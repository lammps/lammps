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

#ifdef GRAN_SUB_MOD_CLASS
// clang-format off
GranSubModStyle(none,
         GranSubModRollingNone,
         ROLLING);

GranSubModStyle(sds,
         GranSubModRollingSDS,
         ROLLING);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_ROLLING_H
#define GRAN_SUB_MOD_ROLLING_H

#include "gran_sub_mod.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GranSubModRolling : public GranSubMod {
 public:
  GranSubModRolling(class GranularModel *, class LAMMPS *);
  ~GranSubModRolling() {};
  virtual void calculate_forces() = 0;
};

/* ---------------------------------------------------------------------- */

class GranSubModRollingNone : public GranSubModRolling {
 public:
  GranSubModRollingNone(class GranularModel *, class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class GranSubModRollingSDS : public GranSubModRolling {
 public:
  GranSubModRollingSDS(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double k, mu, gamma;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_ROLLING_H */
#endif /*GRAN_SUB_MOD_CLASS_H */
