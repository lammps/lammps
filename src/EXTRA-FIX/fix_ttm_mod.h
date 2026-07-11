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
FixStyle(ttm/mod,FixTTMMod);
// clang-format on
#else

#ifndef LMP_FIX_TTM_MOD_H
#define LMP_FIX_TTM_MOD_H

#include "fix.h"

namespace LAMMPS_NS {

struct el_heat_capacity_thermal_conductivity {
  double el_heat_capacity;
  double el_thermal_conductivity;
};

class FixTTMMod : public Fix {
 public:
  FixTTMMod(class LAMMPS *, int, char **);
  ~FixTTMMod() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void post_force_setup(int);
  void post_force_respa_setup(int, int, int);
  void end_of_step() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  double memory_usage() override;
  void grow_arrays(int) override;
  double compute_vector(int) override;

 private:
  int nlevels_respa;
  int seed;
  int outevery;
  double shift;
  char *infile, *outfile;

  class RanMars *random;

  int nxgrid, nygrid, nzgrid;
  int ngridtotal;

  double *gfactor1, *gfactor2, *ratio, **flangevin;
  double ***T_electron, ***T_electron_old, ***T_electron_first;
  double ***net_energy_transfer, ***net_energy_transfer_all;

  double gamma_p, gamma_s, v_0, v_0_sq;
  int skin_layer, surface_l, surface_r, t_surface_l, t_surface_r;
  int movsur;
  double esheat_0, esheat_1, esheat_2, esheat_3, esheat_4, C_limit, electronic_density;
  double el_th_diff, T_damp;
  double intensity, width, duration, surface_double;
  double mult_factor, ttm_dt;
  double pres_factor, free_path, ionic_density;
  double electron_temperature_min;
  el_heat_capacity_thermal_conductivity el_properties(double);
  double el_sp_heat_integral(double);

  void read_parameters(const std::string &);
  void read_electron_temperatures(const std::string &);
  void write_electron_temperatures(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif
