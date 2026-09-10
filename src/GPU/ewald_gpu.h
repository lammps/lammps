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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald/gpu,EwaldGPU);
// clang-format on
#else

#ifndef LMP_EWALD_GPU_H
#define LMP_EWALD_GPU_H

#include "ewald.h"

namespace LAMMPS_NS {

class EwaldGPU : public Ewald {
 public:
  EwaldGPU(class LAMMPS *);
  ~EwaldGPU() override;
  void init() override;
  void setup() override;
  void compute(int, int) override;
  double memory_usage() override;

 protected:
  double cpu_time;    // accumulated host-side (non-GPU) k-space time
};

}    // namespace LAMMPS_NS

#endif
#endif
