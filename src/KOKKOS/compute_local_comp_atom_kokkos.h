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
ComputeStyle(local/comp/atom/kk,ComputeLocalCompAtomKokkos<LMPDeviceType>);
ComputeStyle(local/comp/atom/kk/device,ComputeLocalCompAtomKokkos<LMPDeviceType>);
ComputeStyle(local/comp/atom/kk/host,ComputeLocalCompAtomKokkos<LMPHostType>);
// clang-format on

#else

#ifndef LMP_COMPUTE_LOCAL_COMP_ATOM_KOKKOS_H
#define LMP_COMPUTE_LOCAL_COMP_ATOM_KOKKOS_H

#include "compute_local_comp_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

// clang-format off
struct TagComputeLocalCompAtom {};
// clang-format on

template <class DeviceType> class ComputeLocalCompAtomKokkos : public ComputeLocalCompAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeLocalCompAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeLocalCompAtomKokkos() override;
  void init() override;
  void compute_peratom() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeLocalCompAtom, const int &) const;

 private:

  typename AT::t_x_array x;
  typename ArrayTypes<DeviceType>::t_int_1d type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
  DAT::tdual_float_2d k_result;
  typename AT::t_float_2d d_result;
};

}    // namespace LAMMPS_NS

#endif
#endif
