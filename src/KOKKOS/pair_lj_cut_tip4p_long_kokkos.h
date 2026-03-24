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
PairStyle(lj/cut/tip4p/long/kk,PairLJCutTIP4PLongKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/long/kk/device,PairLJCutTIP4PLongKokkos<LMPDeviceType>);
PairStyle(lj/cut/tip4p/long/kk/host,PairLJCutTIP4PLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_KOKKOS_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_tip4p_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType, unsigned NEIGHFLAG, int ZEROFLAG>
struct PairTIP4PLongComputeFunctor;

template<class DeviceType>
class PairLJCutTIP4PLongKokkos : public PairLJCutTIP4PLong {
 public:
  enum { EnabledNeighFlags = FULL | HALFTHREAD | HALF };
  enum { COUL_FLAG = 1 };
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairLJCutTIP4PLongKokkos(class LAMMPS *);
  ~PairLJCutTIP4PLongKokkos() override;

  void compute(int, int) override;

  void init_tables(double cut_coul, double *cut_respa) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  // Host Kokkos views for TIP4P charge sites and H neighbors (device used in compute)
  Kokkos::DualView<int **, Kokkos::LayoutRight, DeviceType> k_tip4p_hneigh;
  Kokkos::DualView<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> k_tip4p_newsite;
  int nmax_tip4p;

  Kokkos::View<int **, Kokkos::LayoutRight, DeviceType> d_tip4p_hneigh;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> d_tip4p_newsite;

  KK_FLOAT tip4p_cut_coulsqplus;
  KK_FLOAT tip4p_alpha;
  int tip4p_typeO, tip4p_typeH;
  int kokkos_ntypes;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT eval_fpair(bool stack, const KK_FLOAT &rsq, const int &itype, const int &jtype) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT eval_fcoul(bool stack, const KK_FLOAT &rsq, const int &j, const KK_FLOAT &factor_coul,
                      const KK_FLOAT &qtmp) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT eval_ecoul(bool stack, const KK_FLOAT &rsq, const int &j, const KK_FLOAT &factor_coul,
                      const KK_FLOAT &qtmp) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT eval_evdwl(bool stack, const KK_FLOAT &rsq, const int &itype, const int &jtype) const;

  Kokkos::DualView<params_lj_coul **, Kokkos::LayoutRight, DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul **, Kokkos::LayoutRight, DeviceType>::t_dev_const_um params;
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS + 1][MAX_TYPES_STACKPARAMS + 1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS + 1][MAX_TYPES_STACKPARAMS + 1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS + 1][MAX_TYPES_STACKPARAMS + 1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS + 1][MAX_TYPES_STACKPARAMS + 1];

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  DAT::ttransform_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;

  typename AT::t_kkfloat_1d_randomread d_rtable, d_drtable, d_ftable, d_dftable, d_ctable, d_dctable,
      d_etable, d_detable;

  int neighflag;
  int nlocal, nall, eflag, vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT g_ewald_kk;
  KK_FLOAT tabinnersq_kk;

  void allocate() override;

  template<class D, unsigned NF, int ZF>
  friend struct PairTIP4PLongComputeFunctor;

  template<class D>
  friend EV_FLOAT pair_compute_tip4p(PairLJCutTIP4PLongKokkos<D> *, NeighListKokkos<D> *);
};

}    // namespace LAMMPS_NS

#endif
#endif
