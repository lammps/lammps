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
FixStyle(SRP,FixSRP);
// clang-format on
#else

#ifndef LMP_FIX_SRP_H
#define LMP_FIX_SRP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSRP : public Fix {
 public:
  FixSRP(class LAMMPS *, int, char **);
  ~FixSRP() override;
  int setmask() override;
  void init() override;

  void pre_exchange() override;
  void setup_pre_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;
  void post_run() override;

  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int maxsize_restart() override;
  int size_restart(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int modify_param(int, char **) override;

  double **array;

 protected:
  int btype;
  int bptype;
  std::string pair_name;
};

}    // namespace LAMMPS_NS

#endif
#endif
