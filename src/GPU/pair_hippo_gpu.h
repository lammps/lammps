/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hippo/gpu,PairHippoGPU);
// clang-format on
#else

#ifndef LMP_PAIR_HIPPO_GPU_H
#define LMP_PAIR_HIPPO_GPU_H

#include "pair_amoeba.h"

namespace LAMMPS_NS {

class PairHippoGPU : public PairAmoeba {
 public:
  PairHippoGPU(LAMMPS *lmp);
  ~PairHippoGPU() override;
  void compute(int, int) override;
  void init_style() override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

  void induce() override;

  void repulsion() override;
  void dispersion_real() override;
  void multipole_real() override;
  void udirect2b(double **, double **) override;
  void umutual1(double **, double **) override;
  void fphi_uind(FFT_SCALAR ****, double **, double **, double **) override;
  void umutual2b(double **, double **) override;
  void ufield0c(double **, double **) override;
  void polar_real() override;

 private:
  int gpu_mode;
  double cpu_time;
  void *tq_pinned;
  void *fieldp_pinned;
  bool acc_float;

  bool gpu_hal_ready;
  bool gpu_repulsion_ready;
  bool gpu_dispersion_real_ready;
  bool gpu_multipole_real_ready;
  bool gpu_udirect2b_ready;
  bool gpu_umutual1_ready;
  bool gpu_fphi_uind_ready;
  bool gpu_umutual2b_ready;
  bool gpu_polar_real_ready;

  void udirect2b_cpu();

  template<class numtyp>
  void compute_force_from_torque(const numtyp*, double**, double*);
};

}    // namespace LAMMPS_NS
#endif
#endif
