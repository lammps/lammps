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
ComputeStyle(erotate/asphere/kk,ComputeERotateAsphereKokkos<LMPDeviceType>);
ComputeStyle(erotate/asphere/kk/device,ComputeERotateAsphereKokkos<LMPDeviceType>);
ComputeStyle(erotate/asphere/kk/host,ComputeERotateAsphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_EROTATE_ASPHERE_KOKKOS_H
#define LMP_COMPUTE_EROTATE_ASPHERE_KOKKOS_H

#include "compute_erotate_asphere.h"
#include "kokkos_type.h"

#include "atom_vec_ellipsoid_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeERotateAsphereKokkos : public ComputeERotateAsphere {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef AtomVecEllipsoidKokkosBonusArray<DeviceType> EllipBonusAT;

  ComputeERotateAsphereKokkos(class LAMMPS *, int, char **);
  void init() override;
  double compute_scalar() override;

 private:
  class AtomVecEllipsoidKokkos *avecEllipKK;
  typename AT::t_kkfloat_1d_3_randomread omega;
  typename AT::t_kkfloat_1d_3_randomread angmom;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_int_1d_randomread mask;

  typename AT::t_int_1d_randomread ellipsoid;
  typename EllipBonusAT::t_bonus_1d_randomread bonus;
};

}    // namespace LAMMPS_NS

#endif
#endif
