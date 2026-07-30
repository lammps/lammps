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

class FixPIMDNVT : public FixPIMDNVE {
 public:
  FixPIMDNVT(class LAMMPS *, int, char **);
  ~FixPIMDNVT() override;

  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override;
  std::string get_thermo_colname(int) override;

 protected:
  FixPIMDNVT(class LAMMPS *, int, char **, bool);
  void init_nvt_defaults();
  void parse_nvt_arguments(int, char **, const KeywordParser &);
  bool parse_nvt_keyword(int, char **, int &);
  void finish_nuclear_constructor_setup();

  double fixedpoint[3];

  double *eta;
  double *eta_dot;
  double *eta_dotdot;
  double *eta_mass;

  int mtchain;
  int nc_tchain;
  double factor_eta;
  double drag, tdrag_factor;
  double t_freq;
  double t_period;
  double tdof;
  int tdof_override_flag;
  double tdof_override;
  double ke_target;
  double ecouple_work;

  double dthalf, dt4, dt8;
  double *tau_k;
  double pilescale;
  int tstat_flag;

  void nhc_init();
  void o_step();
  void nhc_temp_integrate();
  double compute_nuclear_kinetic_energy() const;
  double chain_target_energy() const;
  virtual bool thermostat_chain_active() const;
  void update_chain0_acceleration(double);
  void propagate_chain_tail_halfstep(double);
  double propagate_chain0_halfstep(double, bool);
  void update_scaled_nuclear_kinetic(double &, double &) const;
  void advance_chain_positions(double);
  void complete_chain0_halfstep(double, double);
  void update_outer_chain_accelerations(double);
  void complete_chain_tail_halfstep(double, double);

  virtual void thermostat_step();
  virtual void force_half_step();
  virtual void centroid_position_half_step();
  virtual void nh_v_temp();
  double thermostat_work_delta(double) const;
  virtual double chain0_target_energy() const;

  void setup_subclass_state() override;
  int base_restart_size() const override;
  int pack_base_restart(double *) const override;
  int unpack_base_restart(const double *) override;
  int nuclear_vector_size() const override;
  double compute_nuclear_vector(int) const override;
};

}    // namespace LAMMPS_NS

#endif
#endif
