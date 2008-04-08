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

#ifndef FIX_NPT_H
#define FIX_NPT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNPT : public Fix {
 public:
  FixNPT(class LAMMPS *, int, char **);
  virtual ~FixNPT();
  int setmask();
  virtual void init();
  void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int);
  double compute_scalar();
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  void reset_dt();

 protected:
  int dimension,which;
  double dtv,dtf,dthalf;
  double boltz,nktv2p;
  double vol0;

  double t_start,t_stop;
  double t_current,t_target;
  double t_freq;
  double f_eta,eta_dot,eta;

  int press_couple,allremap;
  int p_flag[3];                   // 1 if control P on this dim, 0 if not
  double p_start[3],p_stop[3];
  double p_freq[3],p_target[3];
  double omega[3],omega_dot[3];
  double p_current[3],dilation[3];
  double drag,drag_factor;
  double factor[3];
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes

  int nlevels_respa;
  double *step_respa;

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tflag,pflag;

  void couple();
  void remap(int);
};

}

#endif
