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
FixStyle(gle,FixGLE);
// clang-format on
#else

#ifndef LMP_FIX_GLE_H
#define LMP_FIX_GLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGLE : public Fix {
 public:
  FixGLE(class LAMMPS *, int, char **);
  ~FixGLE() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void gle_integrate();
  void initial_integrate_respa(int vflag, int ilevel, int iloop) override;
  void final_integrate_respa(int ilevel, int iloop) override;
  void initial_integrate(int vflag) override;
  void final_integrate() override;
  double compute_scalar() override;
  void reset_target(double) override;
  void reset_dt() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  void *extract(const char *, int &) override;

  void init_gle();
  void init_gles();

 protected:
  int ns, ns1sq;
  double *A, *C, *S, *T, *ST, *TT;
  double *gle_tmp1, *gle_tmp2;
  double t_start, t_stop, t_target;
  double dtv, dtf;

  int dogle, fnoneq, gle_every, gle_step;
  class RanMars *random;
  double *sqrt_m;
  double *step_respa;
  double energy;
  int nlevels_respa;

  double **gle_s;
  double **vaux;
};

}    // namespace LAMMPS_NS

#endif
#endif
