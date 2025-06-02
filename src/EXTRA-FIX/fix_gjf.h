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
FixStyle(gjf,FixGJF);
// clang-format on
#else

#ifndef LMP_FIX_GJF_H
#define LMP_FIX_GJF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGJF : public Fix {
 public:
  FixGJF(class LAMMPS *, int, char **);
  ~FixGJF() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;
  void reset_target(double) override;
  void reset_dt() override;
  int modify_param(int, char **) override;
  double memory_usage() override;
  void *extract(const char *, int &) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  int osflag, tbiasflag, GJmethod, maxatom, lv_allocated;
  double t_start, t_stop, t_period, t_target, tsqrt;
  double gjfc1, gjfc2, gjfc3;
  int tstyle, tvar;
  char *tstr;

  double *tforce;
  double **lv;    //half step velocity

  char *id_temp;
  class Compute *temperature;

  class RanMars *random;
  int seed;

  void compute_target();
};

}    // namespace LAMMPS_NS

#endif
#endif
