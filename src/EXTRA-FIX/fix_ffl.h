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
FixStyle(ffl,FixFFL);
// clang-format on
#else

#ifndef LMP_FIX_FFL_H
#define LMP_FIX_FFL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFFL : public Fix {
 public:
  FixFFL(class LAMMPS *, int, char **);
  ~FixFFL() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void ffl_integrate();
  void initial_integrate_respa(int vflag, int ilevel, int iloop) override;
  void final_integrate_respa(int ilevel, int iloop) override;
  void initial_integrate(int vflag) override;
  void final_integrate() override;
  double compute_scalar() override;
  void reset_target(double) override;
  void reset_dt() override;
  double memory_usage() override;
  void grow_arrays(int) override;

  void *extract(const char *, int &) override;

  void init_ffl();

 protected:
  double *ffl_tmp1, *ffl_tmp2;
  double t_start, t_stop, t_target;
  double dtv, dtf, c1, c2, gamma;
  char flip_type[10];

  int doffl, ffl_every, ffl_step, flip_int;
  class RanMars *random;
  double *sqrt_m;
  double *step_respa;
  double energy;
  int nlevels_respa;

  double **vaux;
};

}    // namespace LAMMPS_NS

#endif
#endif
