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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(rmsd/kk,ComputeRmsdKokkos<LMPDeviceType>);
ComputeStyle(rmsd/kk/device,ComputeRmsdKokkos<LMPDeviceType>);
ComputeStyle(rmsd/kk/host,ComputeRmsdKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_RMSD_KOKKOS_H
#define LMP_COMPUTE_RMSD_KOKKOS_H

#include "compute_rmsd.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeRmsd{};

template<class DeviceType>
class ComputeRmsdKokkos : public ComputeRmsd {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeRmsdKokkos(class LAMMPS *, int, char **);
  ~ComputeRmsdKokkos() override;
  void init() override;
  double compute_scalar() override;

  KOKKOS_INLINE_FUNCTION
  double rmsd(double *);

  KOKKOS_INLINE_FUNCTION
  double gpu_q_j(typename AT::t_x_array, typename AT::t_x_array, double*);

  KOKKOS_INLINE_FUNCTION
  double rmsd_grad_gpu(double*);

 private:
  typename AT::t_x_array d_x_group;

  typename AT::tdual_tagint_1d k_group_taglist;
  typename AT::t_tagint_1d_randomread d_group_taglist;

  typename AT::tdual_x_array k_ref_positions;
  typename AT::t_x_array d_ref_positions;

};

}

#endif
#endif

