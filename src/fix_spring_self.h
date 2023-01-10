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
FixStyle(spring/self,FixSpringSelf);
// clang-format on
#else

#ifndef LMP_FIX_SPRING_SELF_H
#define LMP_FIX_SPRING_SELF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringSelf : public Fix {
 public:
  FixSpringSelf(class LAMMPS *, int, char **);
  ~FixSpringSelf() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  double k, espring;
  double **xoriginal;    // original coords of atoms
  int xflag, yflag, zflag;
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
