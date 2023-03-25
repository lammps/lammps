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
FixStyle(rheo/thermal,FixRHEOThermal)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_THERMAL_H
#define LMP_FIX_RHEO_THERMAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOThermal : public Fix {
 public:
  FixRHEOThermal(class LAMMPS *, int, char **);
  ~FixRHEOThermal() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void initial_integrate() override;
  void post_integrate() override;
  void post_neighbor() override;
  void pre_force(int) override;
  void final_integrate() override;
  void reset_dt() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  double *cv_type, cv;
  double *Tc_type, Tc;
  double *kappa_type, kappa;
  double dtf, dtv;
  int Tc_style;
  int cv_style;
  int conductivity_style;
  int first_flag, last_flag;
  int nmax_old, index_cond;

  class FixRHEO *fix_rheo;

  double calc_cv(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
