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

#ifndef GRAN_SUB_MOD_ROLLING_H
#define GRAN_SUB_MOD_ROLLING_H

#include "gran_sub_mod.h"
#include "gran_sub_mod_rolling_kernel.h"


namespace LAMMPS_NS::Granular_NS {

  class GranSubModRolling : public GranSubMod {
   public:
    GranSubModRolling(class GranularModel *, class LAMMPS *);
    virtual void calculate_forces() = 0;

    void fill_kernel_params(GranKernel::GranRollingParams<double> &p) const
    {
      p.model = GRAN_ROLLING_NONE;    // host classes call the models directly
      p.k = get_k();
      p.gamma = get_gamma();
      p.mu = get_mu();
    }

   protected:
    // only the SDS model defines coefficients, so the base returns zero
    [[nodiscard]] virtual double get_k() const { return 0.0; }
    [[nodiscard]] virtual double get_gamma() const { return 0.0; }
    [[nodiscard]] virtual double get_mu() const { return 0.0; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModRollingNone : public GranSubModRolling {
   public:
    GranSubModRollingNone(class GranularModel *, class LAMMPS *);
    void calculate_forces() override {};
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModRollingSDS : public GranSubModRolling {
   public:
    GranSubModRollingSDS(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    void calculate_forces() override;

   protected:
    double k, mu, gamma;
    [[nodiscard]] double get_k() const override { return k; }
    [[nodiscard]] double get_gamma() const override { return gamma; }
    [[nodiscard]] double get_mu() const override { return mu; }
  };

} // namespace LAMMPS_NS::Granular_NS


#endif /*GRAN_SUB_MOD_ROLLING_H */
