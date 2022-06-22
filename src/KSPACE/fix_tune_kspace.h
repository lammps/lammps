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
FixStyle(tune/kspace,FixTuneKspace);
// clang-format on
#else

#ifndef LMP_FIX_TUNE_KSPACE_H
#define LMP_FIX_TUNE_KSPACE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTuneKspace : public Fix {
 public:
  FixTuneKspace(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void pre_exchange() override;
  double get_timing_info();
  void store_old_kspace_settings();
  void update_pair_style(const std::string &, double);
  void update_kspace_style(const std::string &, const std::string &);
  void adjust_rcut(double);
  void mnbrak();
  void brent0();
  void brent1();
  void brent2();

 private:
  int nevery;

  int last_step;        // previous timestep when timing info was collected
  double last_spcpu;    // old elapsed CPU time value
  int firststep;        // 0 if this is the first time timing info is collected
  int niter;            // number of kspace switches

  double ewald_time, pppm_time, msm_time;
  double pair_cut_coul;
  std::string acc_str;
  std::string kspace_style;
  std::string pair_style;
  std::string base_pair_style;

  int old_differentiation_flag;
  int old_slabflag;
  double old_slab_volfactor;

  int niter_adjust_rcut;
  double ax_brent, bx_brent, cx_brent, dx_brent;
  double fa_brent, fb_brent, fc_brent, fd_brent;
  double v_brent, w_brent, x_brent;
  double fv_brent, fw_brent, fx_brent;
  double a_brent, b_brent;
  double fd2_brent;
  double dxlim;
  bool keep_bracketing, first_brent_pass;
  bool converged, need_fd2_brent;

  inline void shft3(double &a, double &b, double &c, const double d)
  {
    a = b;
    b = c;
    c = d;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
