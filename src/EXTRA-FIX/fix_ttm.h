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
FixStyle(ttm,FixTTM);
// clang-format on
#else

#ifndef LMP_FIX_TTM_H
#define LMP_FIX_TTM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTTM : public Fix {
 public:
  FixTTM(class LAMMPS *, int, char **);
  ~FixTTM() override;
  void post_constructor() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force_setup(int);
  void post_force(int) override;
  void post_force_respa_setup(int, int, int);
  void post_force_respa(int, int, int) override;
  void end_of_step() override;
  void reset_dt() override;
  void grow_arrays(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  double compute_vector(int) override;
  double memory_usage() override;

 protected:
  int nlevels_respa;
  int seed;
  int nxgrid, nygrid, nzgrid;    // size of global grid
  int ngridtotal;                // total size of global grid
  int deallocate_flag;
  int outflag, outevery;
  double shift, tinit;
  double e_energy, transfer_energy;
  char *infile, *outfile;

  class RanMars *random;
  double electronic_specific_heat, electronic_density;
  double electronic_thermal_conductivity;
  double gamma_p, gamma_s, v_0, v_0_sq;

  double *gfactor1, *gfactor2, *ratio, **flangevin;
  double ***T_electron, ***T_electron_old;
  double ***net_energy_transfer, ***net_energy_transfer_all;
  double ***T_atomic;
  int ***nsum, ***nsum_all;
  double ***sum_vsq, ***sum_vsq_all;
  double ***sum_mass_vsq, ***sum_mass_vsq_all;

  virtual void allocate_grid();
  virtual void deallocate_grid();
  virtual void read_electron_temperatures(const std::string &);
  virtual void write_electron_temperatures(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif
