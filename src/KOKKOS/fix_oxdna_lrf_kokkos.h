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

// Read-only strided view into the packed per-atom record of fix OXDNA/LRF/kk
// (position and frame vectors of an atom in one 64-byte row); indexed (i,k)
// like the separate x and nx/ny/nz views it replaces in the force kernels.
template<class DeviceType>
using t_oxdna_packed_sub = Kokkos::View<const KK_FLOAT**, Kokkos::LayoutStride,
  typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr::device_type,
  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

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
  DAT::tdual_kkfloat_1d_3_lr k_nx, k_ny, k_nz;    // LayoutRight: the 3 components of an atom are adjacent
  typename AT::t_kkfloat_1d_3_lr d_nx, d_ny, d_nz;

  // Packed per-atom record for the force kernels, one 64-byte row per atom:
  // columns 0-2 position, 4-6 nx, 7-9 ny, 10-12 nz (3, 13-15 unused), so that
  // all data of a neighbor atom is fetched with one or two memory transactions.
  // Filled for all owned and ghost atoms every time the frames are computed.
  Kokkos::View<KK_FLOAT*[16], Kokkos::LayoutRight, typename AT::t_kkfloat_1d_3_lr::device_type> d_xn;
  t_oxdna_packed_sub<DeviceType> packed_x() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, Kokkos::make_pair(0,3)); }
  t_oxdna_packed_sub<DeviceType> packed_nx() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, Kokkos::make_pair(4,7)); }
  t_oxdna_packed_sub<DeviceType> packed_ny() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, Kokkos::make_pair(7,10)); }
  t_oxdna_packed_sub<DeviceType> packed_nz() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, Kokkos::make_pair(10,13)); }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixOxdnaLRFComputeQuatToXYZ, const int &) const;

 private:

  AtomVecEllipsoidKokkos *avecEllipKK;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkacc_1d_3 f, torque;
  int zero_forces;    // 1 if this fix zeroes f and torque in place of VerletKokkos::force_clear()
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread ellipsoid;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d bonus;

  void compute_lrf_kokkos(int zero_forces_flag);
};

}    // namespace LAMMPS_NS
#endif
#endif
