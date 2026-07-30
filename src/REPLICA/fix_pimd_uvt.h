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
FixStyle(pimd/uvt,FixPIMDUVT);
// clang-format on
#else

#ifndef LMP_FIX_PIMD_UVT_H
#define LMP_FIX_PIMD_UVT_H

#include "arg_info.h"
#include "fix_pimd_nvt.h"

namespace LAMMPS_NS {

class Fix;

class FixPIMDUVT : public FixPIMDNVT {
 public:
  FixPIMDUVT(class LAMMPS *, int, char **);
  ~FixPIMDUVT() override;
  double compute_scalar() override;
  void *extract(const char *, int &) override;

 protected:
  bool parse_uvt_keyword(int, char **, int &);
  void finish_uvt_constructor_setup();
  void setup_subclass_state() override;
  void after_force_transform_hook() override;
  void thermostat_step() override;
  void force_half_step() override;
  void centroid_position_half_step() override;

  bool thermostat_chain_active() const override;
  bool ne_thermostat_participates() const;
  double ne_thermostat_chain_count() const;
  double ne_target_current_share() const;
  double ne_kinetic_current_share() const;
  void scale_ne_velocity(double);
  double chain0_target_energy() const override;

  int subclass_restart_size() const override;
  int pack_subclass_restart(double *, int) const override;
  int unpack_subclass_restart(const double *, int) override;
  int subclass_vector_size() const override;
  double compute_subclass_vector(int) const override;

  void compute_mu_target();
  void nhc_mu_integrate();
  double evaluate_dedn();
  void refresh_dedn_cache();
  void parse_dedn_source(const char *);
  void resolve_dedn_source();

  int ustat_flag;
  int mu_flag;
  double mu;
  double *Ne;
  double *Ne_dot;
  double *Ne_mass;
  double u_start, u_stop;
  double u_current, u_target;
  double u_freq;
  double u_period;
  double ne_ecouple_work;
  char *dedn_name;
  int dedn_which;
  int dedn_index;
  int dedn_var;
  class Compute *dedn_compute;
  Fix *dedn_fix;
  double dedn_current;
};
}    // namespace LAMMPS_NS

#endif
#endif
