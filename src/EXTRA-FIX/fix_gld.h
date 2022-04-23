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
FixStyle(gld,FixGLD);
// clang-format on
#else

#ifndef LMP_FIX_GLD_H
#define LMP_FIX_GLD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGLD : public Fix {
 public:
  FixGLD(class LAMMPS *, int, char **);
  ~FixGLD() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
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
  void init_s_gld();

 protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;
  int freezeflag, zeroflag;
  double t_start, t_stop, t_target;

  int prony_terms;
  int series_type;
  double *prony_c;
  double *prony_tau;

  double **s_gld;

  class RanMars *random;
};

}    // namespace LAMMPS_NS

#endif
#endif
