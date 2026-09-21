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

#ifndef GRAN_SUB_MOD_TANGENTIAL_H
#define GRAN_SUB_MOD_TANGENTIAL_H

#include "gran_sub_mod.h"
#include "gran_sub_mod_tangential_kernel.h"


namespace LAMMPS_NS::Granular_NS {

  class GranSubModTangential : public GranSubMod {
   public:
    GranSubModTangential(class GranularModel *, class LAMMPS *);
    virtual void calculate_forces() = 0;

    [[nodiscard]] double get_k() const { return k; }
    [[nodiscard]] double get_damp() const { return damp; }
    [[nodiscard]] double get_mu() const { return mu; }

    // collect the run-constant coefficients for the shared kernels
    void fill_kernel_params(GranKernel::GranTangentialParams<double> &p) const
    {
      p.model = GRAN_TANGENTIAL_NONE;    // host classes call the models directly
      p.k = k;
      p.xt = xt_kernel();
      p.mu = mu;
      p.mindlin_force = mindlin_force;
      p.mindlin_rescale = mindlin_rescale;
      p.contact_radius_flag = contact_radius_flag;
    }

   protected:
    double k, damp, mu;    // Used by Marshall twisting model
    int mindlin_force, mindlin_rescale;

    // xt lives in the derived classes that define a tangential stiffness
    virtual double xt_kernel() const { return 0.0; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialNone : public GranSubModTangential {
   public:
    GranSubModTangentialNone(class GranularModel *, class LAMMPS *);
    void calculate_forces() override {};
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialLinearNoHistory : public GranSubModTangential {
   public:
    GranSubModTangentialLinearNoHistory(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    void calculate_forces() override;

   protected:
    double xt;
    double xt_kernel() const override { return xt; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialLinearHistory : public GranSubModTangential {
   public:
    GranSubModTangentialLinearHistory(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    void calculate_forces() override;

   protected:
    double xt;
    double xt_kernel() const override { return xt; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialLinearHistoryClassic : public GranSubModTangentialLinearHistory {
   public:
    GranSubModTangentialLinearHistoryClassic(class GranularModel *, class LAMMPS *);
    void calculate_forces() override;
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialMindlinClassic : public GranSubModTangentialLinearHistoryClassic {
   public:
    GranSubModTangentialMindlinClassic(class GranularModel *, class LAMMPS *);
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialMindlin : public GranSubModTangential {
   public:
    GranSubModTangentialMindlin(class GranularModel *, class LAMMPS *);
    void coeffs_to_local() override;
    void mix_coeffs(double *, double *) override;
    void calculate_forces() override;

   protected:
    double xt;
    double xt_kernel() const override { return xt; }
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialMindlinForce : public GranSubModTangentialMindlin {
   public:
    GranSubModTangentialMindlinForce(class GranularModel *, class LAMMPS *);
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialMindlinRescale : public GranSubModTangentialMindlin {
   public:
    GranSubModTangentialMindlinRescale(class GranularModel *, class LAMMPS *);
  };

  /* ---------------------------------------------------------------------- */

  class GranSubModTangentialMindlinRescaleForce : public GranSubModTangentialMindlin {
   public:
    GranSubModTangentialMindlinRescaleForce(class GranularModel *, class LAMMPS *);
  };

} // namespace LAMMPS_NS::Granular_NS


#endif /*GRAN_SUB_MOD_TANGENTIAL_H */
