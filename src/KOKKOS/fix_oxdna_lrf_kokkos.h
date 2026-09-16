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

#ifdef FIX_CLASS
// clang-format off
FixStyle(OXDNA/LRF/kk,FixOxdnaLRFKokkos<LMPDeviceType>);
FixStyle(OXDNA/LRF/kk/device,FixOxdnaLRFKokkos<LMPDeviceType>);
FixStyle(OXDNA/LRF/kk/host,FixOxdnaLRFKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_OXDNA_LRF_KOKKOS_H
#define LMP_FIX_OXDNA_LRF_KOKKOS_H

#include "fix.h"
#include "kokkos_type.h"

#include "atom_vec_ellipsoid_kokkos.h"

namespace LAMMPS_NS {

struct TagFixOxdnaLRFComputeQuatToXYZ{};

template<class DeviceType>
class FixOxdnaLRFKokkos : public Fix {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixOxdnaLRFKokkos(class LAMMPS *, int, char **);
  ~FixOxdnaLRFKokkos() override;

  int setmask() override;
  void init() override;
  void min_setup_pre_force(int);
  void min_pre_force(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  // Unlike vanilla FixOxdnaLRF, we calc nlocal+nghost rather than
  // just nlocal and communicating ghost values via [un]pack routines.
  // So none of these routines are needed here.

  // per-atom arrays for local unit vectors in lab frame
  DAT::tdual_kkfloat_1d_3 k_nx, k_ny, k_nz;
  typename AT::t_kkfloat_1d_3 d_nx, d_ny, d_nz;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaLRFComputeQuatToXYZ, const int &) const;

 private:

  AtomVecEllipsoidKokkos *avecEllipKK;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread ellipsoid;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d bonus;

  void compute_lrf_kokkos();
};

}    // namespace LAMMPS_NS
#endif
#endif
