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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvt/gpu,FixNVTGPU);
// clang-format on
#else

#ifndef LMP_FIX_NVT_GPU_H
#define LMP_FIX_NVT_GPU_H

#include "fix_nh_gpu.h"

namespace LAMMPS_NS {

class FixNVTGPU : public FixNHGPU {
 public:
  FixNVTGPU(class LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
