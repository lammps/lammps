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
GranSubModStyle(none,GranSubModTwistingNone,TWISTING);
GranSubModStyle(marshall,GranSubModTwistingMarshall,TWISTING);
GranSubModStyle(sds,GranSubModTwistingSDS,TWISTING);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_TWISTING_H
#define GRAN_SUB_MOD_TWISTING_H

#include "gran_sub_mod.h"

namespace LAMMPS_NS {
namespace Granular_NS {

  class GranSubModTwisting : public GranSubMod {
   public:
    GranSubModTwisting(class GranularModel *, class LAMMPS *);
    virtual void calculate_forces() = 0;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingNone : public GranSubModTwisting {
   public:
    GranSubModTwistingNone(class GranularModel *, class LAMMPS *);
    void calculate_forces() override{};
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingMarshall : public GranSubModTwisting {
   public:
    GranSubModTwistingMarshall(class GranularModel *, class LAMMPS *);
    void init() override;
    void calculate_forces() override;

   protected:
    double k_tang, mu_tang;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingSDS : public GranSubModTwisting {
   public:
    GranSubModTwistingSDS(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    void calculate_forces() override;

   protected:
    double k, mu, damp;
  };

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_TWISTING_H */
#endif /*GRAN_SUB_MOD_CLASS_H */
