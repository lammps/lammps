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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam/gpu,PairEAMGPU);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_GPU_H
#define LMP_PAIR_EAM_GPU_H

#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMGPU : public PairEAM {
 public:
  PairEAMGPU(class LAMMPS *);
  ~PairEAMGPU() override;
  void compute(int, int) override;
  void init_style() override;
  double single(int, int, int, int, double, double, double, double &) override;
  double memory_usage() override;
  void *extract(const char *, int &) override { return nullptr; }

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  // function-pointer types matching the per-style device instance
  // functions in the gpu library (lal_eam*_ext.cpp)

  typedef int (*EAMGPUInitFn)(const int, double, int **, int **, int *, double ***, double ***,
                              double ***, double **, double, double, double, int, int, int, int,
                              int, const int, const int, const int, const int, const double, int &,
                              FILE *, int &);
  typedef void (*EAMGPUClearFn)();
  typedef int **(*EAMGPUComputeNFn)(const int, const int, const int, double **, int *, double *,
                                    double *, tagint *, int **, tagint **, const bool, const bool,
                                    const bool, const bool, int &, int **, int **, const double,
                                    bool &, int &, void **, double *, int *);
  typedef void (*EAMGPUComputeFn)(const int, const int, const int, const int, double **, int *,
                                  int *, int *, int **, const bool, const bool, const bool,
                                  const bool, int &, const double, bool &, void **);
  typedef void (*EAMGPUComputeForceFn)(int *, const bool, const bool, const bool, const bool);
  typedef double (*EAMGPUBytesFn)();

  // bindings to the device instance of this style, set by the constructor

  EAMGPUInitFn gpu_init_fn;
  EAMGPUClearFn gpu_clear_fn;
  EAMGPUComputeNFn gpu_compute_n_fn;
  EAMGPUComputeFn gpu_compute_fn;
  EAMGPUComputeForceFn gpu_compute_force_fn;
  EAMGPUBytesFn gpu_bytes_fn;

  int gpu_mode;
  double cpu_time;
  void *fp_pinned;
  bool fp_single;
};

}    // namespace LAMMPS_NS

#endif
#endif
