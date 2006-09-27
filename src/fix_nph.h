/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_NPH_H
#define FIX_NPH_H

#include "fix.h"
class Temperature;
class Pressure;

class FixNPH : public Fix {
  friend class Pressure;              // accesses press_couple

 public:
  FixNPH(int, char **);
  ~FixNPH();
  int setmask();
  void init();
  void setup();
  void initial_integrate();
  void final_integrate();
  void initial_integrate_respa(int,int);
  void final_integrate_respa(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  int thermo_fields(int, int *, char **);
  int thermo_compute(double *);

 private:
  Temperature *temperature;
  Pressure *pressure;

  double dtv,dtf,dthalf;
  double boltz,nktv2p;
  double total_mass,vol0,nkt;

  int press_couple,dilate_partial;
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

  void couple();
  void box_dilate(int);
};

#endif
