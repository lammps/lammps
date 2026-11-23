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
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/sf/dipole/sf/kk,PairLJSFDipoleSFKokkos<LMPDeviceType>);
PairStyle(lj/sf/dipole/sf/kk/device,PairLJSFDipoleSFKokkos<LMPDeviceType>);
PairStyle(lj/sf/dipole/sf/kk/host,PairLJSFDipoleSFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_SF_DIPOLE_SF_KOKKOS_H
#define LMP_PAIR_LJ_SF_DIPOLE_SF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_sf_dipole_sf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairLJSFDipoleSFKernel{};

template<class DeviceType>
class PairLJSFDipoleSFKokkos : public PairLJSFDipoleSF {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  PairLJSFDipoleSFKokkos(class LAMMPS *);
  ~PairLJSFDipoleSFKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj_sf_dipole {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_sf_dipole() {cut_ljsq=cut_coulsq=lj1=lj2=lj3=lj4=scale=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_sf_dipole(int /*i*/) {cut_ljsq=cut_coulsq=lj1=lj2=lj3=lj4=scale=0;};
    KK_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,scale;
  };

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int, EV_FLOAT &ev) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int) const;

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
  Kokkos::DualView<params_lj_sf_dipole**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_sf_dipole**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;
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
  KK_FLOAT special_coul[4];
  KK_FLOAT qqrd2e;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  void launch(EV_FLOAT &ev);

  void allocate() override;
  friend void pair_virial_fdotr_compute<PairLJSFDipoleSFKokkos>(PairLJSFDipoleSFKokkos*);
};

}

#endif
#endif
