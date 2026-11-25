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
  void post_force(int) override;
  double compute_scalar() override;
  std::string get_thermo_colname(int) override;
  void *extract(const char *, int &) override;

 protected:
  bool parse_keyword(int, char **, int &) override;
  void setup_subclass_state() override;
  void b_step() override;
  void qc_step() override;

  bool thermostat_chain_active() const override;
  double thermostat_extra_kinetic_energy() const override;
  void scale_ne_velocity(double);
  double chain0_target_energy() const override;

  int subclass_restart_size() const override;
  int pack_subclass_restart(double *, int) const override;
  int unpack_subclass_restart(const double *, int) override;
  double compute_subclass_vector(int) const override;

  void thermostat_extra_velocity_step() override;
  double evaluate_dedn();
  void parse_dedn_source(const char *);
  void resolve_dedn_source();

  double mu;
  double *Ne;
  double *Ne_dot;
  double *Ne_mass;
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
  double dedn_bead;
};
}    // namespace LAMMPS_NS

#endif
#endif
