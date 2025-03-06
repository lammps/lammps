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

#ifndef LMP_FIX_RIGID_NH_SMALL_H
#define LMP_FIX_RIGID_NH_SMALL_H

#include "fix_rigid_small.h"

namespace LAMMPS_NS {

class FixRigidNHSmall : public FixRigidSmall {
 public:
  FixRigidNHSmall(class LAMMPS *, int, char **);
  ~FixRigidNHSmall() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override;
  int modify_param(int, char **) override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  void reset_target(double) override;

 protected:
  double boltz, nktv2p, mvv2e;    // boltzman constant, conversion factors

  int nf_t, nf_r;                       // trans/rot degrees of freedom
  double *w, *wdti1, *wdti2, *wdti4;    // Yoshida-Suzuki coefficients
  double *q_t, *q_r;                    // trans/rot thermostat masses
  double *eta_t, *eta_r;                // trans/rot thermostat positions
  double *eta_dot_t, *eta_dot_r;        // trans/rot thermostat velocities
  double *f_eta_t, *f_eta_r;            // trans/rot thermostat forces

  double epsilon_mass[3], *q_b;         // baro/thermo masses
  double epsilon[3], *eta_b;            // baro/thermo positions
  double epsilon_dot[3], *eta_dot_b;    // baro/thermo velocities
  double *f_eta_b;                      // thermo forces
  double akin_t, akin_r;                // translational/rotational kinetic energies

  int kspace_flag;            // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;    // indices of rigid fixes

  double vol0;          // reference volume
  double t0;            // reference temperature
  int pdim, g_f;        // number of barostatted dims, total DoFs
  double p_hydro;       // hydrostatic target pressure
  double p_freq_max;    // maximum barostat frequency

  double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

  double t_target, t_current;
  double t_freq;

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tcomputeflag, pcomputeflag;    // 1 = compute was created by fix. 0 = external

  void couple();
  void remap();
  void nhc_temp_integrate();
  void nhc_press_integrate();

  virtual void compute_temp_target();
  void compute_press_target();
  void nh_epsilon_dot();
  void compute_dof();

  void allocate_chain();
  void allocate_order();
  void deallocate_chain();
  void deallocate_order();

  inline double maclaurin_series(double);
};

inline double FixRigidNHSmall::maclaurin_series(double x)
{
  double x2, x4;
  x2 = x * x;
  x4 = x2 * x2;
  return (1.0 + (1.0 / 6.0) * x2 + (1.0 / 120.0) * x4 + (1.0 / 5040.0) * x2 * x4 +
          (1.0 / 362880.0) * x4 * x4);
}

}    // namespace LAMMPS_NS

#endif
