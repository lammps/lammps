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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/gpu,PPPMGPU);
// clang-format on
#else

#ifndef LMP_PPPM_GPU_H
#define LMP_PPPM_GPU_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMGPU : public PPPM {
 public:
  PPPMGPU(class LAMMPS *);
  ~PPPMGPU() override;
  void init() override;
  void setup() override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;
  double memory_usage() override;

  void compute_group_group(int, int, int) override;

 protected:
  FFT_SCALAR ***density_brick_gpu, ***vd_brick;
  bool kspace_split, im_real_space;
  int old_nlocal;
  double poisson_time;

  void brick2fft_gpu();
  void poisson_ik() override;

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;

  FFT_SCALAR ***create_3d_offset(int, int, int, int, int, int, const char *, FFT_SCALAR *, int);
  void destroy_3d_offset(FFT_SCALAR ***, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
