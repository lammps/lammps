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

#ifdef BOND_CLASS
// clang-format off
BondStyle(table/kk,BondTableKokkos<LMPDeviceType>);
BondStyle(table/kk/device,BondTableKokkos<LMPDeviceType>);
BondStyle(table/kk/host,BondTableKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_BOND_TABLE_KOKKOS_H
#define LMP_BOND_TABLE_KOKKOS_H

#include "bond_table.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagBondTableCompute{};

template<class DeviceType>
class BondTableKokkos : public BondTable {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  BondTableKokkos(class LAMMPS *);
  ~BondTableKokkos() override;
  void compute(int, int) override;
  void init_style() override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondTableCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondTableCompute<NEWTON_BOND,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void uf_lookup_kk(const int &type, const KK_FLOAT &x, KK_FLOAT &u, KK_FLOAT &mdu) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &ebond, const KK_FLOAT &fbond, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename Kokkos::View<KK_ACC_FLOAT*[3],DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic> > f;
  typename AT::t_int_2d_lr bondlist;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  TransformView<KK_ACC_FLOAT*,double*,Kokkos::LayoutRight,KKDeviceType> k_eatom;
  TransformView<KK_ACC_FLOAT*[6],double*[6],LMPDeviceLayout,KKDeviceType> k_vatom;
  Kokkos::View<KK_ACC_FLOAT*,Kokkos::LayoutRight,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic>> d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6],LMPDeviceLayout,KKDeviceType,Kokkos::MemoryTraits<Kokkos::Atomic>> d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  // device copies of the tabulated data.  the per-table scalars and the
  // tablength-long arrays are packed into views indexed by table number,
  // so the kernel only needs tabindex[type] to reach them

  DAT::tdual_int_1d k_tabindex;
  DAT::tdual_kkfloat_1d k_lo, k_invdelta, k_deltasq6;
  DAT::tdual_kkfloat_2d k_r, k_e, k_de, k_f, k_df, k_e2, k_f2;

  typename AT::t_int_1d d_tabindex;
  typename AT::t_kkfloat_1d d_lo, d_invdelta, d_deltasq6;
  typename AT::t_kkfloat_2d d_r, d_e, d_de, d_f, d_df, d_e2, d_f2;

  DAT::tdual_int_scalar k_error_flag;
  typename AT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  void setup_tables();
};

}

#endif
#endif
