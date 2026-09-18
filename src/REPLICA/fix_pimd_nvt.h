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
FixStyle(pimd,FixPIMDNVT);
FixStyle(pimd/nvt,FixPIMDNVT);
// clang-format on
#else

#ifndef FIX_PIMD_NVT_H
#define FIX_PIMD_NVT_H

#include "fix_pimd_nve.h"

namespace LAMMPS_NS {

class PIMDNoseHoover;

class FixPIMDNVT : public FixPIMDNVE {
  friend class PIMDNoseHoover;

 public:
  FixPIMDNVT(class LAMMPS *, int, char **, bool defer_setup = false);
  ~FixPIMDNVT() override;

  double compute_scalar() override;
  std::string get_thermo_colname(int) override;

 protected:
  bool parse_keyword(int, char **, int &) override;
  void finish_nuclear_constructor_setup();

  // The unregistered FixNH adapter owns these arrays; aliases preserve the
  // existing restart/output layout.
  PIMDNoseHoover *nhc;
  double *eta;
  double *eta_dot;
  double *eta_dotdot;
  double *eta_mass;

  int mtchain;
  int nc_tchain;
  double drag, tdrag_factor;
  double t_freq;
  double t_period;
  double tdof;
  int tdof_override_flag;
  double tdof_override;
  double ke_target;
  double ecouple_work;

  // NHC substep durations; dthalf is a full timestep for BAOAB.
  double dthalf, dt4, dt8;
  double *tau_k;
  double pilescale;
  int tstat_flag;

  void nhc_init();
  double compute_nuclear_kinetic_energy() const;
  double chain_target_energy() const;
  virtual bool thermostat_chain_active() const;
  // Twice the additional kinetic energy coupled to this chain (zero for NVT).
  virtual double thermostat_extra_kinetic_energy() const { return 0.0; }

  void o_step() override;
  virtual void thermostat_step();
  virtual void thermostat_extra_velocity_step() {}
  double thermostat_work_delta(double) const;
  virtual double chain0_target_energy() const;

  void setup_subclass_state() override;
  int base_restart_size() const override;
  int pack_base_restart(double *) const override;
  int unpack_base_restart(const double *) override;
  double compute_subclass_vector(int) const override;
};

}    // namespace LAMMPS_NS

#endif
#endif
