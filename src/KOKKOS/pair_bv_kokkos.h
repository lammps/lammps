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
PairStyle(bv/kk,PairBVKokkos<LMPDeviceType>);
PairStyle(bv/kk/device,PairBVKokkos<LMPDeviceType>);
PairStyle(bv/kk/host,PairBVKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BV_KOKKOS_H
#define LMP_PAIR_BV_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_bv.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagPairBVPackForwardComm{};
struct TagPairBVUnpackForwardComm{};
struct TagPairBVInitialize{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct TagPairBVKernelA{};

template<int EFLAG>
struct TagPairBVKernelB{};

template<int EFLAG>
struct TagPairBVKernelAB{};

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairBVKernelC{};

template<class DeviceType>
class PairBVKokkos : public PairBV, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairBVKokkos(class LAMMPS *);
  ~PairBVKokkos() override;
  void compute(int, int) override;
  void allocate() override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;


  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVInitialize, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelA<NEIGHFLAG,NEWTON_PAIR>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelB<EFLAG>, const int&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelAB<EFLAG>, const int&, EV_FLOAT&) const;

  template<int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelAB<EFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairBVKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_xfloat_1d&, int, int *);
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  struct params_bv{
    KOKKOS_INLINE_FUNCTION
    params_bv() {cut=0;cutsq=0;r0=0;alpha=0;sparam=0;v0=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_bv(int /*i*/) {cut=0;cutsq=0;r0=0;alpha=0;sparam=0;v0=0;offset=0;};
    F_FLOAT cut,cutsq,r0,alpha,sparam,v0,offset;
  };

 protected:

  Kokkos::DualView<params_bv**,Kokkos::LayoutRight,DeviceType> k_params;
  //params is the unmanaged/device view of the dual view
  typename Kokkos::DualView<params_bv**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  //m_params is an instance of params_bv stucture 
  params_bv m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];  
  
  DAT::tdual_ffloat_1d k_energy0;
  typename AT::t_ffloat_1d d_energy0;
  F_FLOAT m_energy0[MAX_TYPES_STACKPARAMS+1];
  
  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;

  typename AT::t_x_array x;
  typename AT::t_f_array f;
  typename AT::t_int_1d type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int need_dup,inum;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout> dup_s0;
  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_f;
  DupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> dup_eatom;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_vatom;
  NonDupScatterView<F_FLOAT*, typename DAT::t_ffloat_1d::array_layout> ndup_s0;
  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_f;
  NonDupScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout> ndup_eatom;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_vatom;

  DAT::tdual_ffloat_1d k_s0;
  DAT::tdual_ffloat_1d k_fp;
  typename AT::t_ffloat_1d d_s0;
  typename AT::t_ffloat_1d d_fp;
  HAT::t_ffloat_1d h_s0;
  HAT::t_ffloat_1d h_fp;
  
  template<class TAG>
  struct policyInstance;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  friend void pair_virial_fdotr_compute<PairBVKokkos>(PairBVKokkos*);
};

}
#endif
#endif

