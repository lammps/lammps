/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include <exception>

namespace LAMMPS_NS {

struct el_heat_capacity_thermal_conductivity {
  double el_heat_capacity;
  double el_thermal_conductivity;
};

class FixTTMMod : public Fix {
 public:
  FixTTMMod(class LAMMPS *, int, char **);
  ~FixTTMMod();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void post_force_setup(int);
  void post_force_respa_setup(int, int, int);
  void end_of_step();
  void reset_dt();
  void write_restart(FILE *);
  void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  double memory_usage();
  void grow_arrays(int);
  double compute_vector(int);

 private:
  int nlevels_respa;
  int seed;
  int outevery;
  double shift;
  char *infile, *outfile;

  class RanMars *random;

  int nxgrid, nygrid, nzgrid;
  int ngridtotal;

  int ***nsum, ***nsum_all;
  double *gfactor1, *gfactor2, *ratio, **flangevin;
  double ***T_electron, ***T_electron_old, ***T_electron_first;
  double ***sum_vsq, ***sum_mass_vsq;
  double ***sum_vsq_all, ***sum_mass_vsq_all;
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

  class parser_error : public std::exception {
    std::string message;

   public:
    parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept { return message.c_str(); }
  };
};

}    // namespace LAMMPS_NS

#endif
#endif
