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
ComputeStyle(hexorder/atom/kk,ComputeHexOrderAtomKokkos<LMPDeviceType>);
ComputeStyle(hexorder/atom/kk/device,ComputeHexOrderAtomKokkos<LMPDeviceType>);
ComputeStyle(hexorder/atom/kk/host,ComputeHexOrderAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_HEXORDER_ATOM_KOKKOS_H
#define LMP_COMPUTE_HEXORDER_ATOM_KOKKOS_H

#include "compute_hexorder_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeHexOrderAtom{};

template<class DeviceType>
class ComputeHexOrderAtomKokkos : public ComputeHexOrderAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeHexOrderAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeHexOrderAtomKokkos() override;
  void init() override;
  void compute_peratom() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeHexOrderAtom, const int&) const;

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

  DAT::ttransform_kkfloat_2d k_qnarray;
  typename AT::t_kkfloat_2d d_qnarray;

  // per-atom scratch, so that the kernel needs no dynamic allocation

  typename AT::t_kkfloat_2d d_distsq;
  typename AT::t_int_2d d_nearest;

  KK_FLOAT cutsq_kk;
  int groupbit_kk;
};

}    // namespace LAMMPS_NS

#endif
#endif
