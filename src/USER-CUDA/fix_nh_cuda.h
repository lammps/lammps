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

#ifndef LMP_FIX_NH_CUDA_H
#define LMP_FIX_NH_CUDA_H

#include "fix.h"
#include "cuda_precision.h"

namespace LAMMPS_NS {

class FixNHCuda : public Fix {
 public:
  FixNHCuda(class LAMMPS *, int, char **);
  virtual ~FixNHCuda();
  int setmask();
  virtual void init();
  void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  double compute_scalar();
  double compute_vector(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  void reset_dt();

 protected:
  class Cuda *cuda;
  int dimension,which;
  double dtv,dtf,dthalf,dt4,dt8,dto;
  double boltz,nktv2p,tdof;
  double vol0,t0;

  double t_start,t_stop;
  double t_current,t_target;
  double t_freq;

  int tstat_flag;                   // 1 if control T
  int pstat_flag;                   // 1 if control P

  int pstyle,pcouple,allremap;
  int p_flag[6];                   // 1 if control P on this dim, 0 if not
  double p_start[6],p_stop[6];
  double p_freq[6],p_target[6];
  double omega[6],omega_dot[6];
  double omega_mass[6];
  double p_current[6],dilation[6];
  double drag,tdrag_factor;        // drag factor on particle thermostat
  double pdrag_factor;             // drag factor on barostat
  double factor[6];                // velocity scaling due to barostat
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int *rfix;                       // indices of rigid fixes

  int nlevels_respa;
  double *step_respa;

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tflag,pflag;

  double *eta,*eta_dot;            // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;
  int mtchain;                     // length of chain
                                   
  double *etap;                    // chain thermostat for barostat
  double *etap_dot;
  double *etap_dotdot;
  double *etap_mass;
  int mpchain;                     // length of chain
                                   
  int mtk_flag;                    // 0 if using Hoover barostat
  double mtk_term1,mtk_term2;
  int mtchain_default_flag;
  int pdim;                        // number of barostatted dims
  double mvv_current[3];           // diagonal of KE tensor
  double mtk_factor;               // MTK factor
  double p_freq_max;               // maximum barostat frequency

  double p_hydro;                  // hydrostatic target pressure

  int nc_tchain,nc_pchain;
  double factor_eta;
  double sigma[6];                 // scaled target stress
  double fdev[6];                  // deviatoric force on barostat
  int deviatoric_flag;             // 0 if target stress tensor is hydrostatic
  double h0_inv[6];                // h_inv of reference (zero strain) box
  int nreset_h0;                   // interval for resetting h0

  void couple();
  void couple_ke();
  void remap();
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void nve_x();            // may be overwritten by child classes
  virtual void nve_v();
  virtual void nh_v_press();
  virtual void nh_v_temp();

  void compute_sigma();
  void compute_deviatoric();
  double compute_strain_energy();
  void compute_press_target();
  void nh_omega_dot();
  
  X_FLOAT triggerneighsq;
};

}

#endif
