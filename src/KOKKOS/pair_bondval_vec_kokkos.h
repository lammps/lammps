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
PairStyle(bondval/vec/kk,PairBondValVecKokkos<LMPDeviceType>);
PairStyle(bondval/vec/kk/device,PairBondValVecKokkos<LMPDeviceType>);
PairStyle(bondval/vec/kk/host,PairBondValVecKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BONDVAL_VEC_KOKKOS_H
#define LMP_PAIR_BONDVAL_VEC_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_bondval_vec.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairBondValVecPackForwardComm{};
struct TagPairBondValVecUnpackForwardComm{};
struct TagPairBondValVecInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairBondValVecKernelA{};

template<int EFLAG>
struct TagPairBondValVecKernelB{};

template<int EFLAG>
struct TagPairBondValVecKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairBondValVecKernelC{};

template<class DeviceType>
class PairBondValVecKokkos : public PairBondValVec, public KokkosBase {
 private:
  int datom;
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairBondValVecKokkos(class LAMMPS *);
  ~PairBondValVecKokkos() override;
  void compute(int, int) override;
  void allocate() override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;


  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelAB<EFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBondValVecKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fx, const KK_FLOAT &fy, const KK_FLOAT &fz, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  struct params_bvv{
    KOKKOS_INLINE_FUNCTION
    params_bvv() {cut=0;cutsq=0;r0=0;alpha=0;bvvsparam=0;bvvv0=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_bvv(int /*i*/) {cut=0;cutsq=0;r0=0;alpha=0;bvvsparam=0;bvvv0=0;offset=0;};
    KK_FLOAT cut,cutsq,r0,alpha,bvvsparam,bvvv0,offset;
  };

 protected:

  Kokkos::DualView<params_bvv**,Kokkos::LayoutRight,DeviceType> k_params;
  //params is the unmanaged/device view of the dual view
  typename Kokkos::DualView<params_bvv**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  //m_params is an instance of params_bv stucture
  params_bvv m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_s0;
  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_s0;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  DAT::tdual_kkacc_1d_3 k_s0;
  DAT::tdual_kkfloat_1d_3 k_Di;

  typename AT::t_kkacc_1d_3 d_s0;
  typename AT::t_kkfloat_1d_3 d_Di;

  HAT::t_kkfloat_1d_3 h_s0;
  HAT::t_kkfloat_1d_3 h_Di;

  template<class TAG>
  struct policyInstance;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_kkfloat_1d_um v_buf;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairBondValVecKokkos>(PairBondValVecKokkos*);
};

}
#endif
#endif

