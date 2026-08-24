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

#ifndef LMP_FIX_NH_TSPIN_H
#define LMP_FIX_NH_TSPIN_H

#include "fix_nve_tspin.h"

namespace LAMMPS_NS {

class FixNHTSpin : public FixNVETSpin {
 public:
  FixNHTSpin(class LAMMPS *, int, char **);
  ~FixNHTSpin() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  double compute_scalar() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 protected:
  int tstat_flag;    // 1 if the temp keyword was used
  double t_start, t_stop, t_target, t_freq, t_period;
  double t_current, ke_target, tdof;
  double drag, tdrag_factor;

  int mtchain, nc_tchain;
  double *eta, *eta_dot, *eta_dotdot, *eta_mass;
  double dtq, dt4, dt8, dthalf;

  void compute_temp_target();
  double compute_spin_temp();
  void nhc_spin_integrate();
  void nh_vs_scale(double);
};

}    // namespace LAMMPS_NS

#endif
