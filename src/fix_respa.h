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
FixStyle(RESPA,FixRespa);
// clang-format on
#else

#ifndef LMP_FIX_RESPA_H
#define LMP_FIX_RESPA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRespa : public Fix {
  friend class Respa;
  friend class FixShake;
  friend class FixRattle;

 public:
  FixRespa(class LAMMPS *, int, char **);
  ~FixRespa() override;
  int setmask() override;
  void init() override {}

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 private:
  int nlevels;
  int store_torque;     // 1 if torques should be stored in addition to forces
  double ***f_level;    // force at each rRESPA level
  double ***t_level;    // torque at each rRESPA level
};

}    // namespace LAMMPS_NS

#endif
#endif
