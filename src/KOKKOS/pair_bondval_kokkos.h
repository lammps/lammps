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
PairStyle(bondval/kk,PairBondValKokkos<LMPDeviceType>);
PairStyle(bondval/kk/device,PairBondValKokkos<LMPDeviceType>);
PairStyle(bondval/kk/host,PairBondValKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BONDVAL_KOKKOS_H
#define LMP_PAIR_BONDVAL_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_bondval.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairBondValPackForwardComm{};
struct TagPairBondValUnpackForwardComm{};
struct TagPairBondValInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairBondValKernelA{};

template<int EFLAG>
struct TagPairBondValKernelB{};

template<int EFLAG>
struct TagPairBondValKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairBondValKernelC{};

template<class DeviceType>
class PairBondValKokkos : public PairBondVal, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairBondValKokkos(class LAMMPS *);
  ~PairBondValKokkos() override;
  void compute(int, int) override;
  void allocate() override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;


  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelAB<EFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  struct params_bv{
    KOKKOS_INLINE_FUNCTION
    params_bv() {cut=0;cutsq=0;r0=0;alpha=0;sparam=0;v0=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_bv(int /*i*/) {cut=0;cutsq=0;r0=0;alpha=0;sparam=0;v0=0;offset=0;};
    KK_FLOAT cut,cutsq,r0,alpha,sparam,v0,offset;
  };

 protected:

  Kokkos::DualView<params_bv**,Kokkos::LayoutRight,DeviceType> k_params;
  //params is the unmanaged/device view of the dual view
  typename Kokkos::DualView<params_bv**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  //m_params is an instance of params_bv stucture
  params_bv m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  DAT::tdual_kkfloat_1d k_energy0;
  typename AT::t_kkfloat_1d d_energy0;
  KK_FLOAT m_energy0[MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::tdual_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  typename AT::t_kkfloat_1d_3 x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d type;

  DAT::tdual_kkacc_1d k_eatom;
  DAT::tdual_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int need_dup,inum;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_s0;
  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_s0;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  DAT::tdual_kkacc_1d k_s0;
  DAT::tdual_kkfloat_1d k_fp;
  typename AT::t_kkacc_1d d_s0;
  typename AT::t_kkfloat_1d d_fp;
  HAT::t_kkacc_1d h_s0;
  HAT::t_kkfloat_1d h_fp;

  template<class TAG>
  struct policyInstance;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairBondValKokkos>(PairBondValKokkos*);
};

}
#endif
#endif

