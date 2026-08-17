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
ComputeStyle(ke/atom/kk,ComputeKEAtomKokkos<LMPDeviceType>);
ComputeStyle(ke/atom/kk/device,ComputeKEAtomKokkos<LMPDeviceType>);
ComputeStyle(ke/atom/kk/host,ComputeKEAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_KE_ATOM_KOKKOS_H
#define LMP_COMPUTE_KE_ATOM_KOKKOS_H

#include "compute_ke_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int RMASS>
struct TagComputeKEAtom{};

template<class DeviceType>
class ComputeKEAtomKokkos : public ComputeKEAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeKEAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeKEAtomKokkos() override;
  void compute_peratom() override;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeKEAtom<RMASS>, const int&) const;

 protected:
  double mvv2e;

  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;

  DAT::ttransform_kkfloat_1d k_ke;
  typename AT::t_kkfloat_1d d_ke;
};

}    // namespace LAMMPS_NS

#endif
#endif
