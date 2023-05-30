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

#ifdef FIX_CLASS
// clang-format off
FixStyle(pafi,FixPAFI);
// clang-format on
#else

#ifndef LMP_FIX_PAFI_H
#define LMP_FIX_PAFI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPAFI : public Fix {
 public:
  FixPAFI(class LAMMPS *, int, char **);
  ~FixPAFI() override;
  int setmask() override;
  void init() override;

  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;

  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_vector(int) override;
  // nve
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void reset_dt() override;

  double memory_usage() override;

 protected:
  int varflag, icompute;
  char *computename;
  class Compute *PathCompute;
  double proj[6], proj_all[6];          // f,v,h, psi
  double results[5], results_all[5];    // f.n, (f.n)**2, psi, dx.n
  double c_v[10], c_v_all[10];
  double temperature, gamma, sqrtD, t_period, local_norm, mass_f;
  int force_flag, od_flag, com_flag;
  int nlevels_respa, ilevel_respa;
  int maxatom;
  class RanMars *random;
  int seed;
  double **h;
  // nve
  double dtv, dtf;
  double *step_respa;
  int mass_require;
};

}    // namespace LAMMPS_NS

#endif
#endif
