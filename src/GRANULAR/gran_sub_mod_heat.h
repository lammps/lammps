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
GranSubModStyle(none,GranSubModHeatNone,HEAT);
GranSubModStyle(radius,GranSubModHeatRadius,HEAT);
GranSubModStyle(area,GranSubModHeatArea,HEAT);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_HEAT_H
#define GRAN_SUB_MOD_HEAT_H

#include "gran_sub_mod.h"

namespace LAMMPS_NS {
namespace Granular_NS {

  class GranSubModHeat : public GranSubMod {
   public:
    GranSubModHeat(class GranularModel *, class LAMMPS *);
    virtual double calculate_heat() = 0;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModHeatNone : public GranSubModHeat {
   public:
    GranSubModHeatNone(class GranularModel *, class LAMMPS *);
    double calculate_heat() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModHeatRadius : public GranSubModHeat {
   public:
    GranSubModHeatRadius(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    double calculate_heat() override;

   protected:
    double conductivity;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModHeatArea : public GranSubModHeat {
   public:
    GranSubModHeatArea(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    double calculate_heat() override;

   protected:
    double heat_transfer_coeff;
  };

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_HEAT_H */
#endif /*GRAN_SUB_MOD_CLASS_H */
