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

#ifndef LMP_FIX_NH_H
#define LMP_FIX_NH_H

#include "fix.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class FixNH : public Fix {
 public:
  FixNH(class LAMMPS *, int, char **);
  ~FixNH() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void pre_force_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void pre_exchange() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  void write_restart(FILE *) override;
  virtual int pack_restart_data(double *);    // pack restart data
  void restart(char *) override;
  int modify_param(int, char **) override;
  void reset_target(double) override;
  void reset_dt() override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
  int dimension, which;
  double dtv, dtf, dthalf, dt4, dt8, dto;
  double boltz, nktv2p, tdof;
  double vol0;    // reference volume
  double t0;      // reference temperature
                  // used for barostat mass
  double t_start, t_stop;
  double t_current, t_target, ke_target;
  double t_freq;

  int tstat_flag;    // 1 if control T
  int pstat_flag;    // 1 if control P

  int pstyle, pcouple, allremap;
  int p_flag[6];    // 1 if control P on this dim, 0 if not
  double p_start[6], p_stop[6];
  double p_freq[6], p_target[6];
  double omega[6], omega_dot[6];
  double omega_mass[6];
  double p_current[6];
  double drag, tdrag_factor;     // drag factor on particle thermostat
  double pdrag_factor;           // drag factor on barostat
  int kspace_flag;               // 1 if KSpace invoked, 0 if not
  int nrigid;                    // number of rigid fixes
  int dilate_group_bit;          // mask for dilation group
  int *rfix;                     // indices of rigid fixes
  char *id_dilate;               // group name to dilate
  class Irregular *irregular;    // for migrating atoms after box flips

  double p_temp;    // target temperature for barostat
  int p_temp_flag;

  int nlevels_respa;
  double *step_respa;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tcomputeflag, pcomputeflag;    // 1 = compute was created by fix
                                     // 0 = created externally

  double *eta, *eta_dot;    // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;
  int mtchain;                 // length of chain
  int mtchain_default_flag;    // 1 = mtchain is default

  double *etap;    // chain thermostat for barostat
  double *etap_dot;
  double *etap_dotdot;
  double *etap_mass;
  int mpchain;    // length of chain

  int mtk_flag;         // 0 if using Hoover barostat
  int pdim;             // number of barostatted dims
  double p_freq_max;    // maximum barostat frequency

  double p_hydro;    // hydrostatic target pressure

  int nc_tchain, nc_pchain;
  double factor_eta;
  double sigma[6];        // scaled target stress
  double fdev[6];         // deviatoric force on barostat
  int deviatoric_flag;    // 0 if target stress tensor is hydrostatic
  double h0_inv[6];       // h_inv of reference (zero strain) box
  int nreset_h0;          // interval for resetting h0

  double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

  int eta_mass_flag;      // 1 if eta_mass updated, 0 if not.
  int omega_mass_flag;    // 1 if omega_mass updated, 0 if not.
  int etap_mass_flag;     // 1 if etap_mass updated, 0 if not.
  int dipole_flag;        // 1 if dipole is updated, 0 if not.
  int dlm_flag;           // 1 if using the DLM rotational integrator, 0 if not

  int scaleyz;     // 1 if yz scaled with lz
  int scalexz;     // 1 if xz scaled with lz
  int scalexy;     // 1 if xy scaled with ly
  int flipflag;    // 1 if box flips are invoked as needed

  int pre_exchange_flag;    // set if pre_exchange needed for box flips

  double fixedpoint[3];    // location of dilation fixed-point

  void couple();
  virtual void remap();
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void nve_x();    // may be overwritten by child classes
  virtual void nve_v();
  virtual void nh_v_press();
  virtual void nh_v_temp();
  virtual void compute_temp_target();
  virtual int size_restart_global();

  void compute_sigma();
  void compute_deviatoric();
  double compute_strain_energy();
  void compute_press_target();
  void nh_omega_dot();
};

}    // namespace LAMMPS_NS

#endif
