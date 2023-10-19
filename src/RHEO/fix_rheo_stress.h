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
FixStyle(rheo/stress,FixRHEOStress);
// clang-format on
#else

#ifndef LMP_FIX_RHEO_STRESS_H
#define LMP_FIX_RHEO_STRESS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOStress : public Fix {
 public:
  FixRHEOStress(class LAMMPS *, int, char **);
  ~FixRHEOStress() override;
  void post_constructor() override;
  int setmask() override;
  void init() override;
  void pre_force(int) override;
  void end_of_step() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  char *id_compute, *id_fix;
  class Compute *stress_compute;
  class FixStoreAtom *store_fix;
};

}    // namespace LAMMPS_NS

#endif
#endif
