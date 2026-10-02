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
ComputeStyle(entropy/atom/kk,ComputeEntropyAtomKokkos<LMPDeviceType>);
ComputeStyle(entropy/atom/kk/device,ComputeEntropyAtomKokkos<LMPDeviceType>);
ComputeStyle(entropy/atom/kk/host,ComputeEntropyAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_ENTROPY_ATOM_KOKKOS_H
#define LMP_COMPUTE_ENTROPY_ATOM_KOKKOS_H

#include "compute_entropy_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int LOCAL>
struct TagComputeEntropyAtom{};

struct TagComputeEntropyAtomAvg{};

template<class DeviceType>
class ComputeEntropyAtomKokkos : public ComputeEntropyAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeEntropyAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeEntropyAtomKokkos() override;
  void init() override;
  void compute_peratom() override;

  template<int LOCAL>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeEntropyAtom<LOCAL>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeEntropyAtomAvg, const int&) const;

 protected:
  int inum;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkfloat_1d k_pair_entropy;
  DAT::ttransform_kkfloat_1d k_pair_entropy_avg;
  typename AT::t_kkfloat_1d d_pair_entropy;
  typename AT::t_kkfloat_1d d_pair_entropy_avg;

  // per-atom scratch for the g(r) histogram and its integrand

  typename AT::t_kkfloat_2d d_gofr;
  typename AT::t_kkfloat_1d d_rbin;
  typename AT::t_kkfloat_1d d_rbinsq;

  KK_FLOAT sigmasq2_kk, density_kk, deltar_kk, cutsq_kk, cutsq2_kk;
  KK_FLOAT local_volume_kk;
  int groupbit_kk;
};

}    // namespace LAMMPS_NS

#endif
#endif
