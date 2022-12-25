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
FixStyle(tmd,FixTMD);
// clang-format on
#else

#ifndef LMP_FIX_TMD_H
#define LMP_FIX_TMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTMD : public Fix {
 public:
  FixTMD(class LAMMPS *, int, char **);
  ~FixTMD() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void initial_integrate_respa(int, int, int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void reset_dt() override;

 private:
  int me;
  int nfileevery, compressed;
  bigint previous_stat;
  FILE *fp;
  double rho_start, rho_stop, rho_old, masstotal;
  double dtv, dtf;
  double *step_respa;
  double work_lambda, work_analytical;
  double **xf, **xold;

  void readfile(char *);
  void open(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif
