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
PairStyle(lj/cut/dipole/cut/kk,PairLJCutDipoleCutKokkos<LMPDeviceType>);
PairStyle(lj/cut/dipole/cut/kk/device,PairLJCutDipoleCutKokkos<LMPDeviceType>);
PairStyle(lj/cut/dipole/cut/kk/host,PairLJCutDipoleCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_DIPOLE_CUT_KOKKOS_H
#define LMP_PAIR_LJ_CUT_DIPOLE_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_dipole_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
struct TagPairLJCutDipoleCutKernel{};

template<class DeviceType>
class PairLJCutDipoleCutKokkos : public PairLJCutDipoleCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  PairLJCutDipoleCutKokkos(class LAMMPS *);
  ~PairLJCutDipoleCutKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutDipoleCutKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int, EV_FLOAT &ev) const;
  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairLJCutDipoleCutKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j, const KK_FLOAT &epair,
                    KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:
  Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_kkfloat_1d_4_randomread mu;
  typename AT::t_kkfloat_1d_4 c_mu;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  DAT::ttransform_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  DAT::ttransform_kkfloat_2d k_cut_coulsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;


  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  void allocate() override;
  friend void pair_virial_fdotr_compute<PairLJCutDipoleCutKokkos>(PairLJCutDipoleCutKokkos*);
};

}

#endif
#endif

