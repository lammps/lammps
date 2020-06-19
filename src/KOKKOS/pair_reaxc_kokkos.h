/* -*- c++ -*- ----------------------------------------------------------

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(reax/c/kk,PairReaxCKokkos<LMPDeviceType>)
PairStyle(reax/c/kk/device,PairReaxCKokkos<LMPDeviceType>)
PairStyle(reax/c/kk/host,PairReaxCKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_REAXC_KOKKOS_H
#define LMP_PAIR_REAXC_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_reaxc.h"
#include "neigh_list_kokkos.h"
#include "reaxc_types.h"

#define C_ele 332.06371
#define SMALL 0.0001
#define KCALpMOL_to_EV 23.02
#define HB_THRESHOLD   1e-2  // 0.01
#define MAX_BONDS      30

#define SQR(x)        ((x)*(x))

namespace LAMMPS_NS {

template<class DeviceType>
struct LR_lookup_table_kk
{
  typedef Kokkos::DualView<cubic_spline_coef*,Kokkos::LayoutRight,DeviceType> tdual_cubic_spline_coef_1d;
  typedef typename tdual_cubic_spline_coef_1d::t_dev t_cubic_spline_coef_1d;

  double dx, inv_dx;

  t_cubic_spline_coef_1d d_vdW, d_CEvd;
  t_cubic_spline_coef_1d d_ele, d_CEclmb;
};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputePolar{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeLJCoulomb{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeTabulatedLJCoulomb{};

struct PairReaxBuildListsFull{};

template<int NEIGHFLAG>
struct PairReaxBuildListsHalf{};

struct PairReaxZero{};

struct PairReaxZeroEAtom{};

struct PairReaxZeroVAtom{};

struct PairReaxBondOrder1{};

struct PairReaxBondOrder2{};

struct PairReaxBondOrder3{};

template<int NEIGHFLAG>
struct PairReaxUpdateBond{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeBond1{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeBond2{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeMulti1{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeMulti2{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeAngular{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeTorsion{};

template<int NEIGHFLAG, int EVFLAG>
struct PairReaxComputeHydrogen{};

struct PairReaxFindBondZero{};

struct PairReaxFindBondSpeciesZero{};

struct PairReaxFindBondSpecies{};


template<class DeviceType>
class PairReaxCKokkos : public PairReaxC {
 public:
  enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT_REAX value_type;

  PairReaxCKokkos(class LAMMPS *);
  virtual ~PairReaxCKokkos();

  void ev_setup(int, int, int alloc = 1);
  void compute(int, int);
  void *extract(const char *, int &);
  void init_style();
  double memory_usage();
  void FindBond(int &);
  void PackBondBuffer(DAT::tdual_ffloat_1d, int &);
  void FindBondSpecies();

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputePolar<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputePolar<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxBuildListsFull, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxBuildListsHalf<NEIGHFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxZeroEAtom, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxZeroVAtom, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxBondOrder1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxBondOrder2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxBondOrder3, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxUpdateBond<NEIGHFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeBond1<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeBond1<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeBond2<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeBond2<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeMulti1<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeMulti2<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeMulti2<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeAngular<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeAngular<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeTorsion<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeTorsion<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxFindBondZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void calculate_find_bond_item(int, int&) const;

  KOKKOS_INLINE_FUNCTION
  void pack_bond_buffer_item(int, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxFindBondSpeciesZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(PairReaxFindBondSpecies, const int&) const;

  struct params_sing{
    KOKKOS_INLINE_FUNCTION
    params_sing(){mass=0;chi=0;eta=0;r_s=0;r_pi=0;r_pi2=0;valency=0;valency_val=0;valency_e=0;valency_boc=0;nlp_opt=0;
      p_lp2=0;p_ovun2=0;p_ovun5=0;p_val3=0;p_val5=0;p_hbond=0;};
    KOKKOS_INLINE_FUNCTION
    params_sing(int i){mass=0;chi=0;eta=0;r_s=0;r_pi=0;r_pi2=0;valency=0;valency_val=0;valency_e=0;valency_boc=0;nlp_opt=0;
      p_lp2=0;p_ovun2=0;p_ovun5=0;p_val3=0;p_val5=0;p_hbond=0;};
    F_FLOAT mass,chi,eta,r_s,r_pi,r_pi2,valency,valency_val,valency_e,valency_boc,nlp_opt,
      p_lp2,p_ovun2,p_ovun5, p_val3, p_val5, p_hbond;
  };

  struct params_twbp{
    KOKKOS_INLINE_FUNCTION
    params_twbp(){gamma=0;gamma_w=0;alpha=0;r_vdw=0;epsilon=0;acore=0;ecore=0;rcore=0;lgre=0;lgcij=0;
      r_s=0;r_pi=0;r_pi2=0;p_bo1=0;p_bo2=0;p_bo3=0;p_bo4=0;p_bo5=0;p_bo6=0;ovc=0;v13cor=0;
      p_boc3=0;p_boc4=0;p_boc5=0;p_be1=0,p_be2=0,De_s=0,De_p=0;De_pp=0;
          p_ovun1=0;};
    KOKKOS_INLINE_FUNCTION
    params_twbp(int i){gamma=0;gamma_w=0;alpha=0;r_vdw=0;epsilon=0;acore=0;ecore=0;rcore=0;lgre=0;lgcij=0;
      r_s=0;r_pi=0;r_pi2=0;p_bo1=0;p_bo2=0;p_bo3=0;p_bo4=0;p_bo5=0;p_bo6=0;ovc=0;v13cor=0;
      p_boc3=0;p_boc4=0;p_boc5=0;p_be1=0,p_be2=0,De_s=0,De_p=0;De_pp=0;
          p_ovun1=0;};
    F_FLOAT gamma,gamma_w,alpha,r_vdw,epsilon,acore,ecore,rcore,lgre,lgcij,
      r_s,r_pi,r_pi2,p_bo1,p_bo2,p_bo3,p_bo4,p_bo5,p_bo6,ovc,v13cor,
      p_boc3,p_boc4,p_boc5,p_be1,p_be2,De_s,De_p,De_pp,
      p_ovun1;
  };

  struct params_thbp{
    KOKKOS_INLINE_FUNCTION
    params_thbp(){cnt=0;theta_00=0;p_val1=0;p_val2=0;p_val4=0;p_val7=0;p_pen1=0;p_coa1=0;};
    KOKKOS_INLINE_FUNCTION
    params_thbp(int i){cnt=0;theta_00=0;p_val1=0;p_val2=0;p_val4=0;p_val7=0;p_pen1=0;p_coa1=0;};
    F_FLOAT cnt, theta_00, p_val1, p_val2, p_val4, p_val7, p_pen1, p_coa1;
  };

  struct params_fbp{
    KOKKOS_INLINE_FUNCTION
    params_fbp(){p_tor1=0;p_cot1=0;V1=0;V2=0;V3=0;};
    KOKKOS_INLINE_FUNCTION
    params_fbp(int i){p_tor1=0;p_cot1=0;V1=0;V2=0;V3=0;};
    F_FLOAT p_tor1, p_cot1, V1, V2, V3;
  };

  struct params_hbp{
    KOKKOS_INLINE_FUNCTION
    params_hbp(){p_hb1=0;p_hb2=0;p_hb3=0;r0_hb=0;};
    KOKKOS_INLINE_FUNCTION
    params_hbp(int i){p_hb1=0;p_hb2=0;p_hb3=0;r0_hb=0;};
    F_FLOAT p_hb1, p_hb2, p_hb3, r0_hb;
  };

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT_REAX &ev, const int &i, const int &j, const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void e_tally(EV_FLOAT_REAX &ev, const int &i, const int &j, const F_FLOAT &epair) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void e_tally_single(EV_FLOAT_REAX &ev, const int &i, const F_FLOAT &epair) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT_REAX &ev, const int &i, F_FLOAT *fi, F_FLOAT *drij) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
    F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally3_atom(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
    F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drjk) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally4(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k, const int &l,
    F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *dril, F_FLOAT *drjl, F_FLOAT *drkl) const;

 protected:
  void cleanup_copy();
  void allocate();
  void allocate_array();
  void setup();
  void init_md();
  int Init_Lookup_Tables();
  void Deallocate_Lookup_Tables();
  void LR_vdW_Coulomb( int i, int j, double r_ij, LR_data *lr );

  typedef Kokkos::DualView<int*,DeviceType> tdual_int_1d;
  Kokkos::DualView<params_sing*,typename DeviceType::array_layout,DeviceType> k_params_sing;
  typename Kokkos::DualView<params_sing*,typename DeviceType::array_layout,DeviceType>::t_dev_const paramssing;

  typedef Kokkos::DualView<int**,DeviceType> tdual_int_2d;
  Kokkos::DualView<params_twbp**,typename DeviceType::array_layout,DeviceType> k_params_twbp;
  typename Kokkos::DualView<params_twbp**,typename DeviceType::array_layout,DeviceType>::t_dev_const paramstwbp;

  typedef Kokkos::DualView<int***,DeviceType> tdual_int_3d;
  Kokkos::DualView<params_thbp***,typename DeviceType::array_layout,DeviceType> k_params_thbp;
  typename Kokkos::DualView<params_thbp***,typename DeviceType::array_layout,DeviceType>::t_dev_const paramsthbp;
  Kokkos::DualView<params_hbp***,typename DeviceType::array_layout,DeviceType> k_params_hbp;
  typename Kokkos::DualView<params_hbp***,typename DeviceType::array_layout,DeviceType>::t_dev_const paramshbp;

  typedef Kokkos::DualView<int****,DeviceType> tdual_int_4d;
  Kokkos::DualView<params_fbp****,typename DeviceType::array_layout,DeviceType> k_params_fbp;
  typename Kokkos::DualView<params_fbp****,typename DeviceType::array_layout,DeviceType>::t_dev_const paramsfbp;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d_randomread tag;
  typename AT::t_float_1d_randomread q;
  typename AT::t_tagint_1d_randomread molecule;

  DAT::tdual_efloat_1d k_eatom;
  typename AT::t_efloat_1d d_eatom;

  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;
  HAT::t_virial_array h_vatom;

  DAT::tdual_float_1d k_tap;
  typename AT::t_float_1d d_tap;
  HAT::t_float_1d h_tap;

  typename AT::t_float_1d d_bo_rij, d_hb_rsq, d_Deltap, d_Deltap_boc, d_total_bo;
  typename AT::t_float_1d d_Delta, d_Delta_boc, d_Delta_lp, d_dDelta_lp, d_Delta_lp_temp, d_CdDelta;
  typename AT::t_ffloat_2d_dl d_BO, d_BO_s, d_BO_pi, d_BO_pi2, d_dBOp;
  typename AT::t_ffloat_2d_dl d_dln_BOp_pix, d_dln_BOp_piy, d_dln_BOp_piz;
  typename AT::t_ffloat_2d_dl d_dln_BOp_pi2x, d_dln_BOp_pi2y, d_dln_BOp_pi2z;
  typename AT::t_ffloat_2d_dl d_C1dbo, d_C2dbo, d_C3dbo;
  typename AT::t_ffloat_2d_dl d_C1dbopi, d_C2dbopi, d_C3dbopi, d_C4dbopi;
  typename AT::t_ffloat_2d_dl d_C1dbopi2, d_C2dbopi2, d_C3dbopi2, d_C4dbopi2;
  typename AT::t_ffloat_2d_dl d_Cdbo, d_Cdbopi, d_Cdbopi2, d_dDeltap_self;

  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_total_bo;
  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_CdDelta;
  Kokkos::Experimental::ScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_eatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_f;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_vatom;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_dDeltap_self;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_Cdbo;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_Cdbopi;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_Cdbopi2;

  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_total_bo;
  Kokkos::Experimental::ScatterView<F_FLOAT*, typename DAT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_CdDelta;
  Kokkos::Experimental::ScatterView<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_eatom;
  Kokkos::Experimental::ScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_f;
  Kokkos::Experimental::ScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_vatom;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_dDeltap_self;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_Cdbo;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_Cdbopi;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_ffloat_2d_dl::array_layout,typename KKDevice<DeviceType>::value,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_Cdbopi2;

  int need_dup;

  typedef Kokkos::DualView<F_FLOAT**[7],typename DeviceType::array_layout,DeviceType> tdual_ffloat_2d_n7;
  typedef typename tdual_ffloat_2d_n7::t_dev_const_randomread t_ffloat_2d_n7_randomread;
  typedef typename tdual_ffloat_2d_n7::t_host t_host_ffloat_2d_n7;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_1d d_bo_first, d_bo_num, d_bo_list, d_hb_first, d_hb_num, d_hb_list;

  DAT::tdual_int_scalar k_resize_bo, k_resize_hb;
  typename AT::t_int_scalar d_resize_bo, d_resize_hb;

  typename AT::t_ffloat_2d_dl d_sum_ovun;
  typename AT::t_ffloat_2d_dl d_dBOpx, d_dBOpy, d_dBOpz;

  int neighflag, newton_pair, maxnumneigh, maxhb, maxbo;
  int nlocal,nall,eflag,vflag;
  F_FLOAT cut_nbsq, cut_hbsq, cut_bosq, bo_cut, thb_cut, thb_cutsq;
  F_FLOAT bo_cut_bond;

  int vdwflag, lgflag;
  F_FLOAT gp[39], p_boc1, p_boc2;

  friend void pair_virial_fdotr_compute<PairReaxCKokkos>(PairReaxCKokkos*);

  int bocnt,hbcnt,enobondsflag;

  typedef LR_lookup_table_kk<DeviceType> LR_lookup_table_kk_DT;

  typedef Kokkos::DualView<LR_lookup_table_kk_DT**,LMPDeviceType::array_layout,DeviceType> tdual_LR_lookup_table_kk_2d;
  typedef typename tdual_LR_lookup_table_kk_2d::t_dev t_LR_lookup_table_kk_2d;

  tdual_LR_lookup_table_kk_2d k_LR;
  t_LR_lookup_table_kk_2d d_LR;

  DAT::tdual_int_2d k_tmpid;
  DAT::tdual_ffloat_2d k_tmpbo;
  DAT::tdual_int_scalar k_error_flag;

  typename AT::t_int_1d d_numneigh_bonds;
  typename AT::t_tagint_2d d_neighid;
  typename AT::t_ffloat_2d d_abo;

  typename AT::t_ffloat_1d d_buf;
  DAT::tdual_int_scalar k_nbuf_local;
};

template <class DeviceType>
struct PairReaxCKokkosFindBondFunctor  {
  typedef DeviceType device_type;
  typedef int value_type;
  PairReaxCKokkos<DeviceType> c;
  PairReaxCKokkosFindBondFunctor(PairReaxCKokkos<DeviceType>* c_ptr):c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void join(volatile int &dst,
             const volatile int &src) const {
    dst = MAX(dst,src);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &numbonds) const {
    c.calculate_find_bond_item(ii,numbonds);
  }
};

template <class DeviceType>
struct PairReaxCKokkosPackBondBufferFunctor  {
  typedef DeviceType device_type;
  typedef int value_type;
  PairReaxCKokkos<DeviceType> c;
  PairReaxCKokkosPackBondBufferFunctor(PairReaxCKokkos<DeviceType>* c_ptr):c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &j, const bool &final) const {
    c.pack_bond_buffer_item(ii,j,final);
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
