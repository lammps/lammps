/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(tune/kspace,FixTuneKspace)

#else

#ifndef LMP_FIX_TUNE_KSPACE_H
#define LMP_FIX_TUNE_KSPACE_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixTuneKspace : public Fix {
 public:
  FixTuneKspace(class LAMMPS *, int, char **);
  ~FixTuneKspace();
  int setmask();
  void init();
  void pre_exchange();
  double get_timing_info();
  void store_old_kspace_settings();
  void update_pair_style(char *, double);
  void update_kspace_style(char *, char *);
  void adjust_rcut(double);
  void mnbrak();
  void brent0();
  void brent1();
  void brent2();

 private:
  int nevery;

  int last_step;      // previous timestep when timing info was collected
  double last_spcpu;  // old elapsed CPU time value
  int firststep;      // 0 if this is the first time timing info is collected
  int niter;          // number of kspace switches

  double ewald_time,pppm_time,msm_time;
  double pair_cut_coul;
  char new_acc_str[12];
  char new_kspace_style[20];
  char new_pair_style[20];
  char base_pair_style[20];

  int old_differentiation_flag;
  int old_slabflag;
  double old_slab_volfactor;

  int niter_adjust_rcut;
  double ax_brent,bx_brent,cx_brent,dx_brent;
  double fa_brent,fb_brent,fc_brent,fd_brent;
  double v_brent,w_brent,x_brent;
  double fv_brent,fw_brent,fx_brent;
  double a_brent,b_brent;
  double fd2_brent;
  double dxlim;
  bool keep_bracketing,first_brent_pass;
  bool converged,need_fd2_brent;

  inline void shft3(double &a, double &b, double &c, const double d)
  {
    a=b;
    b=c;
    c=d;
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix tune/kspace without kspace

This fix (tune/kspace) can only be used when a kspace style has been specified.

E: Cannot use fix tune/kspace without a pair style

This fix (tune/kspace) can only be used when a pair style has been specified.

E: Bad real space Coulomb cutoff in fix tune/kspace

Fix tune/kspace tried to find the optimal real space Coulomb cutoff using
the Newton-Rhaphson method, but found a non-positive or NaN cutoff

*/
