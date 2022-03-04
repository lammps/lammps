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
FixStyle(lb/pc,FixLbPC);
// clang-format on
#else

#ifndef LMP_FIX_LB_PC_H
#define LMP_FIX_LB_PC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLbPC : public Fix {
 public:
  FixLbPC(class LAMMPS *, int, char **);
  ~FixLbPC() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  //  void set_arrays(int);
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 private:
  double dtv, dtf;
  int me;
  double *Gamma_MD;
  double expminusdttimesgamma;
  double DMDcoeff;

  double **force_old;
  double **up;
  double **up_old;

  void compute_up();
  class FixLbFluid *fix_lb_fluid;
};

}    // namespace LAMMPS_NS

#endif
#endif
