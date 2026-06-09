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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ylz/kk,PairYLZKokkos<LMPDeviceType>);
PairStyle(ylz/kk/device,PairYLZKokkos<LMPDeviceType>);
PairStyle(ylz/kk/host,PairYLZKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_YLZ_KOKKOS_H
#define LMP_PAIR_YLZ_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_ylz.h"
#include "neigh_list_kokkos.h"
#include "atom_vec_ellipsoid_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
struct TagPairYLZKernel{};

template<class DeviceType>
class PairYLZKokkos : public PairYLZ {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef AtomVecEllipsoidKokkosBonusArray<DeviceType> EllipBonusAT;
  typedef EV_FLOAT value_type;

  PairYLZKokkos(class LAMMPS *);
  ~PairYLZKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairYLZKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairYLZKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j, const KK_FLOAT &epair,
                    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  struct params_ylz {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_ylz() {cutsq=0;epsilon=0;sigma=0;zeta=0;mu=0;beta=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_ylz(int /*i*/) {cutsq=0;epsilon=0;sigma=0;zeta=0;mu=0;beta=0;};
    KK_FLOAT cutsq,epsilon,sigma,zeta,mu,beta;
  };

 protected:
  Kokkos::DualView<params_ylz**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_ylz**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um d_params;
  params_ylz m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread d_ellipsoid;
  typename EllipBonusAT::t_bonus_1d_randomread bonus;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int newton_pair;
  int nlocal,nall,eflag,vflag;
  KK_FLOAT special_lj[4];

  int neighflag;

  class AtomVecEllipsoidKokkos *avecKK;

  void allocate() override;
  friend void pair_virial_fdotr_compute<PairYLZKokkos>(PairYLZKokkos*);
};

}    // namespace LAMMPS_NS

#endif
#endif
