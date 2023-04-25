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
GranSubModStyle(none,GranSubModDampingNone,DAMPING);
GranSubModStyle(velocity,GranSubModDampingVelocity,DAMPING);
GranSubModStyle(mass_velocity,GranSubModDampingMassVelocity,DAMPING);
GranSubModStyle(viscoelastic,GranSubModDampingViscoelastic,DAMPING);
GranSubModStyle(tsuji,GranSubModDampingTsuji,DAMPING);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_DAMPING_H
#define GRAN_SUB_MOD_DAMPING_H

#include "gran_sub_mod.h"
#include "pointers.h"

namespace LAMMPS_NS {
namespace Granular_NS {

  class GranSubModDamping : public GranSubMod {
   public:
    GranSubModDamping(class GranularModel *, class LAMMPS *);
    void init() override;
    virtual double calculate_forces() = 0;
    double get_damp_prefactor() const { return damp_prefactor; }

   protected:
    double damp_prefactor;
    double damp;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModDampingNone : public GranSubModDamping {
   public:
    GranSubModDampingNone(class GranularModel *, class LAMMPS *);
    void init() override{};
    double calculate_forces() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModDampingVelocity : public GranSubModDamping {
   public:
    GranSubModDampingVelocity(class GranularModel *, class LAMMPS *);
    double calculate_forces() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModDampingMassVelocity : public GranSubModDamping {
   public:
    GranSubModDampingMassVelocity(class GranularModel *, class LAMMPS *);
    double calculate_forces() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModDampingViscoelastic : public GranSubModDamping {
   public:
    GranSubModDampingViscoelastic(class GranularModel *, class LAMMPS *);
    double calculate_forces() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModDampingTsuji : public GranSubModDamping {
   public:
    GranSubModDampingTsuji(class GranularModel *, class LAMMPS *);
    void init() override;
    double calculate_forces() override;
  };

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_DAMPING_H */
#endif /*GRAN_SUB_MOD_CLASS_H */
