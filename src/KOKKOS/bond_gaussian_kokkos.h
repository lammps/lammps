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
BondStyle(gaussian/kk,BondGaussianKokkos<LMPDeviceType>);
BondStyle(gaussian/kk/device,BondGaussianKokkos<LMPDeviceType>);
BondStyle(gaussian/kk/host,BondGaussianKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_BOND_GAUSSIAN_KOKKOS_H
#define LMP_BOND_GAUSSIAN_KOKKOS_H

#include "bond_gaussian.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagBondGaussianCompute{};

template<class DeviceType>
class BondGaussianKokkos : public BondGaussian {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  BondGaussianKokkos(class LAMMPS *);
  ~BondGaussianKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondGaussianCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagBondGaussianCompute<NEWTON_BOND,EVFLAG>, const int&) const;

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

  // per-type 1D parameters
  DAT::tdual_kkfloat_1d k_bond_temperature;
  DAT::tdual_int_1d k_nterms;
  typename AT::t_kkfloat_1d d_bond_temperature;
  typename AT::t_int_1d d_nterms;

  // per-type 2D parameters [ntypes+1][max_nterms]
  DAT::tdual_kkfloat_2d k_alpha;
  DAT::tdual_kkfloat_2d k_width;
  DAT::tdual_kkfloat_2d k_r0;
  typename AT::t_kkfloat_2d d_alpha;
  typename AT::t_kkfloat_2d d_width;
  typename AT::t_kkfloat_2d d_r0;

  int kk_max_nterms;

  void allocate() override;
  void update_2d_views(int ntypes);
};

}

#endif
#endif
