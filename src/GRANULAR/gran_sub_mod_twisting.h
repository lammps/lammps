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

#ifndef GRAN_SUB_MOD_TWISTING_H
#define GRAN_SUB_MOD_TWISTING_H

#include "gran_sub_mod.h"
#include "gran_sub_mod_twisting_kernel.h"


namespace LAMMPS_NS::Granular_NS {

  class GranSubModTwisting : public GranSubMod {
   public:
    GranSubModTwisting(class GranularModel *, class LAMMPS *);
    virtual double calculate_forces() = 0;

    void fill_kernel_params(GranKernel::GranTwistingParams<double> &p) const
    {
      p.model = GRAN_TWISTING_NONE;    // host classes call the models directly
      p.k = get_k();
      p.damp = get_damp();
      p.mu = get_mu();
      p.k_tang = get_k_tang();
      p.mu_tang = get_mu_tang();
    }

   protected:
    // each model defines only part of these, so the base returns zero
    [[nodiscard]] virtual double get_k() const { return 0.0; }
    [[nodiscard]] virtual double get_damp() const { return 0.0; }
    [[nodiscard]] virtual double get_mu() const { return 0.0; }
    [[nodiscard]] virtual double get_k_tang() const { return 0.0; }
    [[nodiscard]] virtual double get_mu_tang() const { return 0.0; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingNone : public GranSubModTwisting {
   public:
    GranSubModTwistingNone(class GranularModel *, class LAMMPS *);
    double calculate_forces() override {return 0.0;};
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingMarshall : public GranSubModTwisting {
   public:
    GranSubModTwistingMarshall(class GranularModel *, class LAMMPS *);
    void init() override;
    double calculate_forces() override;

   protected:
    double k_tang, mu_tang;
    [[nodiscard]] double get_k_tang() const override { return k_tang; }
    [[nodiscard]] double get_mu_tang() const override { return mu_tang; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTwistingSDS : public GranSubModTwisting {
   public:
    GranSubModTwistingSDS(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    double calculate_forces() override;

   protected:
    double k, mu, damp;
    [[nodiscard]] double get_k() const override { return k; }
    [[nodiscard]] double get_damp() const override { return damp; }
    [[nodiscard]] double get_mu() const override { return mu; }
  };

} // namespace LAMMPS_NS::Granular_NS


#endif /*GRAN_SUB_MOD_TWISTING_H */
