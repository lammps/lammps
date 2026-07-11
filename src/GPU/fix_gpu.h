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
FixStyle(GPU,FixGPU);
// clang-format on
#else

#ifndef LMP_FIX_GPU_H
#define LMP_FIX_GPU_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGPU : public Fix {
 public:
  FixGPU(class LAMMPS *, int, char **);
  ~FixGPU() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  void post_force_respa(int, int, int) override;
  double memory_usage() override;

  double binsize(const double subx, const double suby, const double subz, const int nlocal,
                 const double cut);

 private:
  int _gpu_mode;
  int _nlevels_respa;
  double _particle_split;
  double _binsize;
};

}    // namespace LAMMPS_NS

#endif
#endif
