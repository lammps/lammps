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

#ifndef LMP_FIX_RIGID_NH_H
#define LMP_FIX_RIGID_NH_H

#include "fix_rigid.h"

namespace LAMMPS_NS {

class FixRigidNH : public FixRigid {
 public:
  FixRigidNH(class LAMMPS *, int, char **);
  virtual ~FixRigidNH();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual double compute_scalar();
  int modify_param(int, char **);
  void write_restart(FILE *);
  void restart(char *buf);
  void reset_target(double);
  
 protected:
  double **conjqm;                    // conjugate quaternion momentum
  double boltz,nktv2p,mvv2e;          // boltzman constant, conversion factors

  int nf_t,nf_r;                      // trans/rot degrees of freedom
  double onednft,onednfr;             // factors 1 + dimension/trans(rot) 
                                      //   degrees of freedom
  double *w,*wdti1,*wdti2,*wdti4;     // Yoshida-Suzuki coefficients
  double *q_t,*q_r;                   // trans/rot thermostat masses
  double *eta_t,*eta_r;               // trans/rot thermostat positions
  double *eta_dot_t,*eta_dot_r;       // trans/rot thermostat velocities
  double *f_eta_t,*f_eta_r;           // trans/rot thermostat forces
  
  double epsilon_mass[3], *q_b;       // baro/thermo masses
  double epsilon[3],*eta_b;           // baro/thermo positions
  double epsilon_dot[3],*eta_dot_b;   // baro/thermo velocities
  double *f_eta_b;                    // thermo forces
  double akin_t,akin_r;               // translational/rotational kinetic energies
  
  int kspace_flag;                    // 1 if KSpace invoked, 0 if not
  int nrigidfix;                      // number of rigid fixes
  int *rfix;                          // indicies of rigid fixes

  double vol0;                        // reference volume
  double t0;                          // reference temperature
  int pdim,g_f;                       // number of barostatted dims, total DoFs
  double p_hydro;                     // hydrostatic target pressure
  double p_freq_max;                  // maximum barostat frequency
  
  double mtk_term1,mtk_term2;         // Martyna-Tobias-Klein corrections
  
  double t_current,t_target;
  double p_current[3],p_target[3];

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tcomputeflag,pcomputeflag;

  void couple();
  void remap();  
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void compute_temp_target();
  void compute_press_target();
  void nh_epsilon_dot();

  void allocate_chain();
  void allocate_order();
  void deallocate_chain();
  void deallocate_order();

  inline double maclaurin_series(double);
};

inline double FixRigidNH::maclaurin_series(double x)
{
  double x2,x4;
  x2 = x * x;
  x4 = x2 * x2;
  return (1.0 + (1.0/6.0) * x2 + (1.0/120.0) * x4 + (1.0/5040.0) * x2 * x4 +
         (1.0/362880.0) * x4 * x4);
}

}

#endif

/* ERROR/WARNING messages:

E: Fix rigid npt/nph period must be > 0.0

Self-explanatory.

E: Invalid fix rigid npt/nph command for a 2d simulation

Cannot control z dimension in a 2d model.

E: Invalid fix rigid npt/nph command pressure settings

If multiple dimensions are coupled, those dimensions must be
specified.

E: Cannot use fix rigid npt/nph on a non-periodic dimension

When specifying a diagonal pressure component, the dimension must be
periodic.

E: Invalid fix rigid npt/nph pressure settings

Settings for coupled dimensions must be the same.

E: Fix rigid nvt/npt/nph damping parameters must be > 0.0

Self-explanatory.

E: Fix rigid npt/nph dilate group ID does not exist

Self-explanatory.

E: Temperature ID for fix rigid nvt/npt/nph does not exist

Self-explanatory.

E: Fix rigid npt/nph does not yet allow triclinic box

Self-explanatory.

E: Cannot use fix rigid npt/nph and fix deform on same component of stress tensor

This would be changing the same box dimension twice.

E: Pressure ID for fix rigid npt/nph does not exist

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

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

*/
