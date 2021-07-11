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
  virtual ~FixGLE();
  int setmask();
  void init();
  void setup(int);
  void gle_integrate();
  void initial_integrate_respa(int vflag, int ilevel, int iloop);
  void final_integrate_respa(int ilevel, int iloop);
  void initial_integrate(int vflag);
  void final_integrate();
  double compute_scalar();
  void reset_target(double);
  virtual void reset_dt();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

  virtual void *extract(const char *, int &);

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
