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

/* ----------------------------------------------------------------------
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mesomem/dipole/kk,PairMesomemDipoleKokkos<LMPDeviceType>);
PairStyle(mesomem/dipole/kk/device,PairMesomemDipoleKokkos<LMPDeviceType>);
PairStyle(mesomem/dipole/kk/host,PairMesomemDipoleKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MESOMEM_DIPOLE_KOKKOS_H
#define LMP_PAIR_MESOMEM_DIPOLE_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_mesomem_dipole.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairMesomemDipoleKernel{};

template<class DeviceType>
class PairMesomemDipoleKokkos : public PairMesomemDipole {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  PairMesomemDipoleKokkos(class LAMMPS *);
  ~PairMesomemDipoleKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_mesomem {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_mesomem() {sigma=eps=zeta=gscale=ktilt=ksplay=c0=wc=wc_inv=wc_half2inv=0;zpow=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_mesomem(int /*i*/) {sigma=eps=zeta=gscale=ktilt=ksplay=c0=wc=wc_inv=wc_half2inv=0;zpow=0;};
    KK_FLOAT sigma,eps,zeta,gscale,ktilt,ksplay,c0,wc,wc_inv,wc_half2inv;
    int zpow;
  };

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int, EV_FLOAT &ev) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j, const KK_FLOAT &epair,
                    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:
  Kokkos::DualView<params_mesomem**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_mesomem**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_4_randomread mu;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag,inum;

  KK_FLOAT special_lj[4];

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  void launch(EV_FLOAT &ev);

  void allocate() override;
  friend void pair_virial_fdotr_compute<PairMesomemDipoleKokkos>(PairMesomemDipoleKokkos*);
};

}

#endif
#endif
