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
ComputeStyle(erotate/sphere/atom/kk,ComputeErotateSphereAtomKokkos<LMPDeviceType>);
ComputeStyle(erotate/sphere/atom/kk/device,ComputeErotateSphereAtomKokkos<LMPDeviceType>);
ComputeStyle(erotate/sphere/atom/kk/host,ComputeErotateSphereAtomKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_EROTATE_SPHERE_ATOM_KOKKOS_H
#define LMP_COMPUTE_EROTATE_SPHERE_ATOM_KOKKOS_H

#include "compute_erotate_sphere_atom.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagComputeErotateSphereAtom{};

template<class DeviceType>
class ComputeErotateSphereAtomKokkos : public ComputeErotateSphereAtom {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeErotateSphereAtomKokkos(class LAMMPS *, int, char **);
  ~ComputeErotateSphereAtomKokkos() override;
  void compute_peratom() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeErotateSphereAtom, const int&) const;

 protected:
  typename AT::t_kkfloat_1d_3_randomread omega;
  typename AT::t_kkfloat_1d_randomread radius;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_int_1d_randomread mask;

  DAT::ttransform_kkfloat_1d k_erot;
  typename AT::t_kkfloat_1d d_erot;
};

}    // namespace LAMMPS_NS

#endif
#endif
