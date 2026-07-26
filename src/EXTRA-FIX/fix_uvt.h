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
FixStyle(uvt,FixUVT);
// clang-format on
#else

#ifndef LMP_FIX_UVT_H
#define LMP_FIX_UVT_H

#include "arg_info.h"
#include "fix_nh.h"

namespace LAMMPS_NS {

class Compute;
class Fix;

class FixUVT : public FixNH {
 public:
  FixUVT(class LAMMPS *, int, char **);
  ~FixUVT() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  std::string get_thermo_colname(int) override;
  int size_restart_global() override;
  int pack_restart_data(double *) override;
  void restart(char *) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
  void nve_v() override;
  void nve_x() override;

 private:
  void nhc_mu_integrate();
  void compute_mu_target();
  double evaluate_dedn();
  void parse_dedn_source(const char *);

  double u_start, u_stop;
  double u_current, u_target;
  double u_freq;
  int ustat_flag;
  double *Ne;
  double *Ne_dot;
  double *Ne_mass;
  char *dedn_name;
  int dedn_which;
  int dedn_index;
  int dedn_var;
  Compute *dedn_compute;
  Fix *dedn_fix;
  double dedn_current;
  int dedn_defer;
};

}    // namespace LAMMPS_NS

#endif
#endif
