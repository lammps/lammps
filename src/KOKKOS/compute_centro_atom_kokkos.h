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
ComputeStyle(centro/atom/kk,ComputeCentroAtomKokkos<LMPDeviceType>);
ComputeStyle(centro/atom/kk/device,ComputeCentroAtomKokkos<LMPDeviceType>);
ComputeStyle(centro/atom/kk/host,ComputeCentroAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_CENTRO_ATOM_KOKKOS_H
#define LMP_COMPUTE_CENTRO_ATOM_KOKKOS_H

#include "compute_centro_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int AXES>
struct TagComputeCentroAtom{};

template<class DeviceType>
class ComputeCentroAtomKokkos : public ComputeCentroAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeCentroAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeCentroAtomKokkos() override;
  void init() override;
  void compute_peratom() override;

  template<int AXES>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeCentroAtom<AXES>, const int&) const;

  // device versions of the two quickselect helpers of the CPU style

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void select_kk(int k, int n, int ii) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void select2_kk(int k, int n, int ii) const;

 protected:
  int inum;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkfloat_1d k_centro;
  typename AT::t_kkfloat_1d d_centro;
  DAT::ttransform_kkfloat_2d k_array_atom;
  typename AT::t_kkfloat_2d d_array_atom;

  // per-atom scratch, so that the kernel needs no dynamic allocation

  typename AT::t_kkfloat_2d d_distsq;
  typename AT::t_int_2d d_nearest;
  typename AT::t_kkfloat_2d d_pairs;

  KK_FLOAT cutsq_kk;
  int groupbit_kk, npairs, nhalf;
};

}    // namespace LAMMPS_NS

#endif
#endif
