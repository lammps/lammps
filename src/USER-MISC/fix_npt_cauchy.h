/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(npt/cauchy,FixNPTCauchy)

#else

#ifndef LMP_FIX_NPT_CAUCHY_H
#define LMP_FIX_NPT_CAUCHY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNPTCauchy : public Fix {
 public:
  FixNPTCauchy(class LAMMPS *, int, char **);
  virtual void post_constructor();
  virtual ~FixNPTCauchy();
  int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  virtual void pre_exchange();
  double compute_scalar();
  virtual double compute_vector(int);
  void write_restart(FILE *);
  virtual int pack_restart_data(double *); // pack restart data
  virtual void restart(char *);
  int modify_param(int, char **);
  void reset_target(double);
  void reset_dt();
  virtual void *extract(const char*,int &);
  double memory_usage();

 protected:
  int dimension,which;
  double dtv,dtf,dthalf,dt4,dt8,dto;
  double boltz,nktv2p,tdof;
  double vol0;                      // reference volume
  double t0;                        // reference temperature
                                    // used for barostat mass
  double t_start,t_stop;
  double t_current,t_target,ke_target;
  double t_freq;

  int tstat_flag;                   // 1 if control T
  int pstat_flag;                   // 1 if control P

  int pstyle,pcouple,allremap;
  int p_flag[6];                   // 1 if control P on this dim, 0 if not
  double p_start[6],p_stop[6];
  double p_freq[6],p_target[6];
  double omega[6],omega_dot[6];
  double omega_mass[6];
  double p_current[6];
  double drag,tdrag_factor;        // drag factor on particle thermostat
  double pdrag_factor;             // drag factor on barostat
  int kspace_flag;                 // 1 if KSpace invoked, 0 if not
  int nrigid;                      // number of rigid fixes
  int dilate_group_bit;            // mask for dilation group
  int *rfix;                       // indices of rigid fixes
  char *id_dilate;                 // group name to dilate
  class Irregular *irregular;      // for migrating atoms after box flips

  int nlevels_respa;
  double *step_respa;

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tcomputeflag,pcomputeflag;   // 1 = compute was created by fix
                                   // 0 = created externally

  double *eta,*eta_dot;            // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;
  int mtchain;                     // length of chain
  int mtchain_default_flag;        // 1 = mtchain is default

  double *etap;                    // chain thermostat for barostat
  double *etap_dot;
  double *etap_dotdot;
  double *etap_mass;
  int mpchain;                     // length of chain

  int mtk_flag;                    // 0 if using Hoover barostat
  int pdim;                        // number of barostatted dims
  double p_freq_max;               // maximum barostat frequency

  double p_hydro;                  // hydrostatic target pressure

  int nc_tchain,nc_pchain;
  double factor_eta;
  double sigma[6];                 // scaled target stress
  double fdev[6];                  // deviatoric force on barostat
  int deviatoric_flag;             // 0 if target stress tensor is hydrostatic
  double h0_inv[6];                // h_inv of reference (zero strain) box
  int nreset_h0;                   // interval for resetting h0

  double mtk_term1,mtk_term2;      // Martyna-Tobias-Klein corrections

  int eta_mass_flag;               // 1 if eta_mass updated, 0 if not.
  int omega_mass_flag;             // 1 if omega_mass updated, 0 if not.
  int etap_mass_flag;              // 1 if etap_mass updated, 0 if not.
  int dipole_flag;                 // 1 if dipole is updated, 0 if not.
  int dlm_flag;                    // 1 if using the DLM rotational integrator, 0 if not

  int scaleyz;                     // 1 if yz scaled with lz
  int scalexz;                     // 1 if xz scaled with lz
  int scalexy;                     // 1 if xy scaled with ly
  int flipflag;                    // 1 if box flips are invoked as needed

  int pre_exchange_flag;           // set if pre_exchange needed for box flips

  double fixedpoint[3];            // location of dilation fixed-point

  void couple();
  virtual void remap();
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void nve_x();            // may be overwritten by child classes
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

  // Implementation of CauchyStat
  char *id_store;       // fix id of the STORE fix for retaining data
  class FixStore *init_store;  // fix pointer to STORE fix
  double H0[3][3];      // shape matrix for the undeformed cell
  double h_old[6];      // previous time step shape matrix for
                        // the undeformed cell
  double invH0[3][3];   // inverse of H0;
  double CSvol0;        // volume of undeformed cell
  double setPK[3][3];   // current set values of the PK stress
                        // (this is modified until the cauchy
                        // stress converges)
  double alpha;         // integration parameter for the cauchystat
  int initPK;           // 1 if setPK needs to be initialized either
                        // from cauchy or restart, else 0
  int restartPK;        // Read PK stress from the previous run
  int restart_stored;   // values of PK stress from the previous step stored
  int initRUN;          // 0 if run not initialized
                        // (pressure->vector not computed yet),
                        // else 1 (pressure->vector available)

  void CauchyStat_init();
  void CauchyStat_cleanup();
  void CauchyStat();
  void CauchyStat_Step(double (&Fi)[3][3], double (&Fdot)[3][3],
                       double (&cauchy)[3][3], double (&setcauchy)[3][3],
                       double (&setPK)[3][3], double volume, double volume0,
                       double deltat, double alpha);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Target temperature for fix npt/cauchy cannot be 0.0

Self-explanatory.

E: Invalid fix npt/cauchy command for a 2d simulation

Cannot control z dimension in a 2d model.

E: Fix npt/cauchy dilate group ID does not exist

Self-explanatory.

E: Invalid fix npt/cauchy command pressure settings

If multiple dimensions are coupled, those dimensions must be
specified.

E: Cannot use fix npt/cauchy on a non-periodic dimension

When specifying a diagonal pressure component, the dimension must be
periodic.

E: Cannot use fix npt/cauchy on a 2nd non-periodic dimension

When specifying an off-diagonal pressure component, the 2nd of the two
dimensions must be periodic.  E.g. if the xy component is specified,
then the y dimension must be periodic.

E: Cannot use fix npt/cauchy with yz scaling when z is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix npt/cauchy with xz scaling when z is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix npt/cauchy with xy scaling when y is non-periodic dimension

The 2nd dimension in the barostatted tilt factor must be periodic.

E: Cannot use fix npt/cauchy with both yz dynamics and yz scaling

Self-explanatory.

E: Cannot use fix npt/cauchy with both xz dynamics and xz scaling

Self-explanatory.

E: Cannot use fix npt/cauchy with both xy dynamics and xy scaling

Self-explanatory.

E: Can not specify Pxy/Pxz/Pyz in fix npt/cauchy with non-triclinic box

Only triclinic boxes can be used with off-diagonal pressure components.
See the region prism command for details.

E: Invalid fix npt/cauchy pressure settings

Settings for coupled dimensions must be the same.

E: Using update dipole flag requires atom style sphere

Self-explanatory.

E: Using update dipole flag requires atom attribute mu

Self-explanatory.

E: Fix npt/cauchy damping parameters must be > 0.0

Self-explanatory.

E: Cannot use fix npt/cauchy and fix deform on same component of stress tensor

This would be changing the same box dimension twice.

E: Temperature ID for fix npt/cauchy does not exist

Self-explanatory.

E: Pressure ID for fix npt/cauchy does not exist

Self-explanatory.

E: Non-numeric pressure - simulation unstable

UNDOCUMENTED

E: Fix npt/cauchy has tilted box too far in one step - periodic cell is too far from equilibrium state

Self-explanatory.  The change in the box tilt is too extreme
on a short timescale.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Temperature for fix modify is not for group all

The temperature compute is being used with a pressure calculation
which does operate on group all, so this may be inconsistent.

E: Pressure ID for fix modify does not exist

Self-explanatory.

E: Could not find fix_modify pressure ID

The compute ID for computing pressure does not exist.

E: Fix_modify pressure ID does not compute pressure

The compute ID assigned to the fix must compute pressure.

U: The dlm flag must be used with update dipole

Self-explanatory.

*/
