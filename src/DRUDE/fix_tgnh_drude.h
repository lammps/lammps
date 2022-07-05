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

#ifndef LMP_FIX_TGNH_DRUDE_H
#define LMP_FIX_TGNH_DRUDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTGNHDrude : public Fix {
 public:
  FixTGNHDrude(class LAMMPS *, int, char **);
  ~FixTGNHDrude() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void pre_force_respa(int, int, int) override;
  void initial_integrate_respa(int, int, int) override;
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
  double memory_usage() override;

 protected:
  int dimension, which;
  double dtv, dtf, dthalf, dt4, dt8, dto;
  double boltz, nktv2p, tdof;
  double vol0;    // reference volume
  double t0;      // reference temperature
                  // used for barostat mass
  double t_start, t_stop;
  double t_current, t_target;
  double t_freq;

  int tstat_flag;    // 1 if control T
  int pstat_flag;    // 1 if control P

  int pstyle, pcouple;
  int p_flag[6];    // 1 if control P on this dim, 0 if not
  double p_start[6], p_stop[6];
  double p_freq[6], p_target[6];
  double omega[6], omega_dot[6];
  double omega_mass[6];
  double p_current[6];
  int kspace_flag;               // 1 if KSpace invoked, 0 if not
  int nrigid;                    // number of rigid fixes
  int *rfix;                     // indices of rigid fixes
  class Irregular *irregular;    // for migrating atoms after box flips

  int nlevels_respa;
  double *step_respa;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tcomputeflag, pcomputeflag;    // 1 = compute was created by fix
                                     // 0 = created externally

  double *etamol;
  double *etamol_dot;    // chain thermostat for motion of whole molecules
  double *etamol_dotdot;
  double *etamol_mass;

  double *etaint;
  double *etaint_dot;    // chain thermostat for internal DOFs
  double *etaint_dotdot;
  double *etaint_mass;

  double *etadrude;
  double *etadrude_dot;    // chain thermostat for Drude relative motions
  double *etadrude_dotdot;
  double *etadrude_mass;

  double *etap;    // chain thermostat for barostat
  double *etap_dot;
  double *etap_dotdot;
  double *etap_mass;

  int mtchain;    // length of chain
  int mpchain;    // length of chain

  int mtk_flag;         // 0 if using Hoover barostat
  int pdim;             // number of barostatted dims
  double p_freq_max;    // maximum barostat frequency

  double p_hydro;    // hydrostatic target pressure

  int nc_tchain, nc_pchain;
  double sigma[6];        // scaled target stress
  double fdev[6];         // deviatoric force on barostat
  int deviatoric_flag;    // 0 if target stress tensor is hydrostatic
  double h0_inv[6];       // h_inv of reference (zero strain) box
  int nreset_h0;          // interval for resetting h0

  double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

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

  class FixDrude *fix_drude;
  int n_mol;    // number of molecules in the system
  double *mass_mol;
  double dof_mol, dof_int, dof_drude;    // DOFs of different modes in the fix group
  void setup_mol_mass_dof();
  double **v_mol, **v_mol_tmp;
  void compute_temp_mol_int_drude(bool);    // calculate the temperatures of three sets of DOFs
  bool temp_computed_end_of_step = false;
  double tdrude_target, tdrude_freq;
  double t_mol, t_int, t_drude;
  double ke2mol, ke2int, ke2drude;
  double ke2mol_target, ke2int_target, ke2drude_target;
  double factor_eta_mol, factor_eta_int, factor_eta_drude;
  double propagate(double *, double *, double *, const double *, const double &, const double &,
                   const double &) const;
};

}    // namespace LAMMPS_NS

#endif
