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

// one column of the packed per-atom record (e.g. the atom type)
template<class DeviceType>
using t_oxdna_packed_col = Kokkos::View<const KK_FLOAT*, Kokkos::LayoutStride,
  typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr::device_type,
  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

// the whole packed record, for kernels that load complete rows
template<class DeviceType>
using t_oxdna_packed = Kokkos::View<const KK_FLOAT*[16], Kokkos::LayoutRight,
  typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr::device_type,
  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

// one row of the packed record: x (0-2), type (3), nx (4-6), ny (7-9), nz (10-12), qeff (13)
struct OxdnaRow {
  KK_FLOAT v[16];
};

// load the first ncol (a multiple of 16 bytes) values of row i with 16-byte vector loads
template<int NCOL, class ViewType>
KOKKOS_INLINE_FUNCTION
void oxdna_load_row(const ViewType &xn, const int i, OxdnaRow &r)
{
  struct alignas(16) Chunk16 { KK_FLOAT v[16 / sizeof(KK_FLOAT)]; };
  constexpr int nper = 16 / sizeof(KK_FLOAT);
  static_assert(NCOL % nper == 0, "oxdna_load_row: NCOL must fill whole 16-byte chunks");
  const Chunk16 *src = reinterpret_cast<const Chunk16 *>(xn.data() + (size_t) i * 16);
  for (int k = 0; k < NCOL / nper; k++) {
    const Chunk16 c = src[k];
    for (int m = 0; m < nper; m++) r.v[nper * k + m] = c.v[m];
  }
}

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

  // Packed per-atom record for the force kernels, one 64-byte row per atom:
  // columns 0-2 position, 3 type, 4-6 nx, 7-9 ny, 10-12 nz, 13 qeff (14-15
  // unused), so that
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
  // atom type (exact as KK_FLOAT) and effective charge
  t_oxdna_packed_col<DeviceType> packed_type() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, 3); }
  t_oxdna_packed_col<DeviceType> packed_qeff() const
    { return Kokkos::subview(d_xn, Kokkos::ALL, 13); }
  t_oxdna_packed<DeviceType> packed() const { return d_xn; }

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
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread qeff;
  typename AtomVecEllipsoidKokkosBonusArray<DeviceType>::t_bonus_1d bonus;

  void compute_lrf_kokkos(int zero_forces_flag);

  // write one 16-value row of d_xn; rows are 64 B (single and mixed precision)
  // or 128 B (double) and aligned, so they are written as 16-byte chunks
  struct alignas(16) Chunk16 { KK_FLOAT v[16 / sizeof(KK_FLOAT)]; };

  KOKKOS_INLINE_FUNCTION
  void store_row(const int i, const KK_FLOAT *row) const
  {
    constexpr int nper = 16 / sizeof(KK_FLOAT);
    Chunk16 *dst = reinterpret_cast<Chunk16 *>(&d_xn(i, 0));
    for (int k = 0; k < 16 / nper; k++) {
      Chunk16 c;
      for (int m = 0; m < nper; m++) c.v[m] = row[nper * k + m];
      dst[k] = c;
    }
  }
};

}    // namespace LAMMPS_NS
#endif
#endif
