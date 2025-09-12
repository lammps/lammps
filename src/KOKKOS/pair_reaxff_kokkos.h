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
PairStyle(reaxff/kk,PairReaxFFKokkos<LMPDeviceType>);
PairStyle(reaxff/kk/device,PairReaxFFKokkos<LMPDeviceType>);
PairStyle(reaxff/kk/host,PairReaxFFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_REAXC_KOKKOS_H
#define LMP_PAIR_REAXC_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_reaxff.h"
#include "neigh_list_kokkos.h"

#include "reaxff_inline.h"

namespace LAMMPS_NS {

template<class DeviceType>
struct LR_lookup_table_kk
{
  typedef Kokkos::DualView<ReaxFF::cubic_spline_coef*,Kokkos::LayoutRight,DeviceType> tdual_cubic_spline_coef_1d;
  typedef typename tdual_cubic_spline_coef_1d::t_dev t_cubic_spline_coef_1d;

  double dx, inv_dx;

  t_cubic_spline_coef_1d d_vdW, d_CEvd;
  t_cubic_spline_coef_1d d_ele, d_CEclmb;
};

template<int NEIGHFLAG>
struct TagPairReaxComputePolar{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairReaxComputeLJCoulomb{};

struct TagPairReaxComputeLJCoulombShortNeigh{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairReaxComputeTabulatedLJCoulomb{};

template<int NEIGHFLAG>
struct TagPairReaxBuildListsHalfBlocking{};

template<int NEIGHFLAG>
struct TagPairReaxBuildListsHalfBlockingPreview{};

template<int NEIGHFLAG>
struct TagPairReaxBuildListsHalfPreview{};

struct TagPairReaxBuildListsFull{};

struct TagPairReaxZero{};

struct TagPairReaxBondOrder1{};

struct TagPairReaxBondOrder2{};

struct TagPairReaxBondOrder3{};

template<int NEIGHFLAG>
struct TagPairReaxUpdateBond{};

template<int NEIGHFLAG, int EFLAG>
struct TagPairReaxComputeBond1{};

template<int NEIGHFLAG, int VFLAG>
struct TagPairReaxComputeBond2{};

struct TagPairReaxComputeMulti1{};

template<int NEIGHFLAG, int EFLAG>
struct TagPairReaxComputeMulti2{};

template<bool POPULATE>
struct TagPairReaxCountAngularTorsion{};
template<int NEIGHFLAG, int EVFLAG>
struct TagPairReaxComputeAngularPreprocessed{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairReaxComputeTorsionPreprocessed{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairReaxComputeHydrogen{};

struct TagPairReaxFindBondZero{};

struct TagPairReaxFindBondSpeciesZero{};

struct TagPairReaxFindBondSpecies{};


template<class DeviceType>
class PairReaxFFKokkos : public PairReaxFF {
 public:
  enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT_REAX value_type;

  // "Blocking" factors to reduce thread divergence within some kernels
  using blocking_t = unsigned short int;

  // "PairReaxFFComputeTorsion"
  static constexpr int compute_torsion_blocksize = 8;

  // "PairReaxBuildListsHalfBlocking"
  static constexpr int build_lists_half_blocksize = 64;

  PairReaxFFKokkos(class LAMMPS *);
  virtual ~PairReaxFFKokkos();

  void ev_setup(int, int, int alloc = 1);
  void compute(int, int);
  void init_style();
  double memory_usage();
  void FindBond(int &, int groupbit = 1);
  void PackBondBuffer(DAT::tdual_double_1d, int &);
  void PackReducedBondBuffer(DAT::tdual_double_1d, int &, bool);
  void FindBondSpecies();

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeLJCoulombShortNeigh, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBuildListsHalfBlocking<NEIGHFLAG>, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBuildListsHalfBlockingPreview<NEIGHFLAG>, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBuildListsHalfPreview<NEIGHFLAG>, const int&) const;

  // Isolated function that builds the hbond list, reused across
  // TagPairReaxBuildListsHalfBlocking, HalfBlockingPreview, HalfPreview
  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void build_hb_list(KK_FLOAT, int, int, int, int) const;

  // Isolated function that builds the bond order list, reused across
  // TagPairReaxBuildListsHalfBlocking, HalfBlockingPreview, HalfPreview
  // Returns if we need to populate d_d* functions or not
  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  bool build_bo_list(int, int, int&, int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBuildListsFull, const int&) const;

  // Isolated function that computes bond order parameters
  // Returns BO_s, BO_pi, BO_pi2, C12, C34, C56 by reference
  KOKKOS_INLINE_FUNCTION
  void compute_bo(KK_FLOAT, int, int, KK_FLOAT, KK_FLOAT, KK_FLOAT,
    KK_FLOAT&, KK_FLOAT&, KK_FLOAT&, KK_FLOAT&, KK_FLOAT&, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBondOrder1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBondOrder2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxBondOrder3, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxUpdateBond<NEIGHFLAG>, const int&) const;

  template<int NEIGHFLAG, int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int&) const;

  template<int NEIGHFLAG, int VFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int VFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeMulti1, const int&) const;

  template<int NEIGHFLAG, int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int&) const;

  template<bool POPULATE>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxCountAngularTorsion<POPULATE>, const int&) const;

  // Abstraction for computing SBSO2, CSBO2, dSBO1, dsBO2
  KOKKOS_INLINE_FUNCTION
  void compute_angular_sbo(int, int, int) const;

  // Abstraction for counting and populating angular intermediates
  template<bool POPULATE>
  KOKKOS_INLINE_FUNCTION
  int preprocess_angular(int, int, int, int) const;

  // Abstraction for counting and populating torsion intermediated
  template<bool POPULATE>
  KOKKOS_INLINE_FUNCTION
  int preprocess_torsion(int, int, tagint, KK_FLOAT, KK_FLOAT, KK_FLOAT, int, int) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT_REAX&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxFindBondZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void calculate_find_bond_item(int, int&, int) const;

  KOKKOS_INLINE_FUNCTION
  void pack_bond_buffer_item(int, int&, const bool&) const;

  template<bool STORE_BONDS>
  KOKKOS_INLINE_FUNCTION
  void pack_reduced_bond_buffer_item(int, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxFindBondSpeciesZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairReaxFindBondSpecies, const int&) const;

  struct params_sing {
    KOKKOS_INLINE_FUNCTION
    params_sing() {mass=0;chi=0;eta=0;r_s=0;r_pi=0;r_pi2=0;valency=0;valency_val=0;valency_e=0;valency_boc=0;nlp_opt=0;
      p_lp2=0;p_ovun2=0;p_ovun5=0;p_val3=0;p_val5=0;p_hbond=0;bcut_acks2=0;};
    KOKKOS_INLINE_FUNCTION
    params_sing(int /*i*/) {mass=0;chi=0;eta=0;r_s=0;r_pi=0;r_pi2=0;valency=0;valency_val=0;valency_e=0;valency_boc=0;nlp_opt=0;
      p_lp2=0;p_ovun2=0;p_ovun5=0;p_val3=0;p_val5=0;p_hbond=0;bcut_acks2=0;};
    KK_FLOAT mass,chi,eta,r_s,r_pi,r_pi2,valency,valency_val,valency_e,valency_boc,nlp_opt,
      p_lp2,p_ovun2,p_ovun5, p_val3, p_val5, p_hbond, bcut_acks2;
  };

  struct params_twbp {
    KOKKOS_INLINE_FUNCTION
    params_twbp() {gamma=0;gamma_w=0;alpha=0;r_vdw=0;epsilon=0;acore=0;ecore=0;rcore=0;lgre=0;lgcij=0;
      r_s=0;r_pi=0;r_pi2=0;p_bo1=0;p_bo2=0;p_bo3=0;p_bo4=0;p_bo5=0;p_bo6=0;ovc=0;v13cor=0;
      p_boc3=0;p_boc4=0;p_boc5=0;p_be1=0,p_be2=0,De_s=0,De_p=0;De_pp=0;
          p_ovun1=0;};
    KOKKOS_INLINE_FUNCTION
    params_twbp(int /*i*/) {gamma=0;gamma_w=0;alpha=0;r_vdw=0;epsilon=0;acore=0;ecore=0;rcore=0;lgre=0;lgcij=0;
      r_s=0;r_pi=0;r_pi2=0;p_bo1=0;p_bo2=0;p_bo3=0;p_bo4=0;p_bo5=0;p_bo6=0;ovc=0;v13cor=0;
      p_boc3=0;p_boc4=0;p_boc5=0;p_be1=0,p_be2=0,De_s=0,De_p=0;De_pp=0;
          p_ovun1=0;};
    KK_FLOAT gamma,gamma_w,alpha,r_vdw,epsilon,acore,ecore,rcore,lgre,lgcij,
      r_s,r_pi,r_pi2,p_bo1,p_bo2,p_bo3,p_bo4,p_bo5,p_bo6,ovc,v13cor,
      p_boc3,p_boc4,p_boc5,p_be1,p_be2,De_s,De_p,De_pp,
      p_ovun1;
  };

  struct params_thbp {
    KOKKOS_INLINE_FUNCTION
    params_thbp() {cnt=0;theta_00=0;p_val1=0;p_val2=0;p_val4=0;p_val7=0;p_pen1=0;p_coa1=0;};
    KOKKOS_INLINE_FUNCTION
    params_thbp(int /*i*/) {cnt=0;theta_00=0;p_val1=0;p_val2=0;p_val4=0;p_val7=0;p_pen1=0;p_coa1=0;};
    KK_FLOAT cnt, theta_00, p_val1, p_val2, p_val4, p_val7, p_pen1, p_coa1;
  };

  struct params_fbp {
    KOKKOS_INLINE_FUNCTION
    params_fbp() {p_tor1=0;p_cot1=0;V1=0;V2=0;V3=0;};
    KOKKOS_INLINE_FUNCTION
    params_fbp(int /*i*/) {p_tor1=0;p_cot1=0;V1=0;V2=0;V3=0;};
    KK_FLOAT p_tor1, p_cot1, V1, V2, V3;
  };

  struct params_hbp {
    KOKKOS_INLINE_FUNCTION
    params_hbp() {p_hb1=0;p_hb2=0;p_hb3=0;r0_hb=0;};
    KOKKOS_INLINE_FUNCTION
    params_hbp(int /*i*/) {p_hb1=0;p_hb2=0;p_hb3=0;r0_hb=0;};
    KK_FLOAT p_hb1, p_hb2, p_hb3, r0_hb;
  };

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT_REAX &ev, const int &i, const int &j, const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void e_tally(EV_FLOAT_REAX &ev, const int &i, const int &j, const KK_FLOAT &epair) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void e_tally_single(EV_FLOAT_REAX &ev, const int &i, const KK_FLOAT &epair) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally(EV_FLOAT_REAX &ev, const int &i, KK_ACC_FLOAT *fi, KK_FLOAT *drij) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
    KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drij, KK_FLOAT *drik) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally3_atom(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
    KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drjk) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally4(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k, const int &l,
    KK_ACC_FLOAT *fi, KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *dril, KK_FLOAT *drjl, KK_FLOAT *drkl) const;

 protected:
  void deallocate_views_of_views();
  void allocate();
  void allocate_array();
  void setup();
  void init_md();
  int Init_Lookup_Tables();
  void Deallocate_Lookup_Tables();
  void LR_vdW_Coulomb(int i, int j, double r_ij, ReaxFF::LR_data *lr);

  Kokkos::DualView<params_sing*,LMPDeviceLayout,DeviceType> k_params_sing;
  typename Kokkos::DualView<params_sing*,LMPDeviceLayout,DeviceType>::t_dev_const paramssing;

  Kokkos::DualView<params_twbp**,LMPDeviceLayout,DeviceType> k_params_twbp;
  typename Kokkos::DualView<params_twbp**,LMPDeviceLayout,DeviceType>::t_dev_const paramstwbp;

  Kokkos::DualView<params_thbp***,LMPDeviceLayout,DeviceType> k_params_thbp;
  typename Kokkos::DualView<params_thbp***,LMPDeviceLayout,DeviceType>::t_dev_const paramsthbp;
  Kokkos::DualView<params_hbp***,LMPDeviceLayout,DeviceType> k_params_hbp;
  typename Kokkos::DualView<params_hbp***,LMPDeviceLayout,DeviceType>::t_dev_const paramshbp;

  Kokkos::DualView<params_fbp****,LMPDeviceLayout,DeviceType> k_params_fbp;
  typename Kokkos::DualView<params_fbp****,LMPDeviceLayout,DeviceType>::t_dev_const paramsfbp;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d_randomread tag;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_tagint_1d_randomread molecule;

  DAT::ttransform_kkacc_1d k_eatom;
  typename AT::t_kkacc_1d d_eatom;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  DAT::tdual_kkfloat_1d k_tap;
  typename AT::t_kkfloat_1d d_tap;

  typename AT::t_kkfloat_1d d_bo_rij, d_hb_rsq, d_Deltap, d_Deltap_boc, d_total_bo, d_s;
  typename AT::t_kkfloat_1d d_Delta, d_Delta_boc, d_Delta_lp, d_dDelta_lp, d_Delta_lp_temp, d_CdDelta;
  typename AT::t_kkfloat_2d_dl d_BO, d_BO_s, d_BO_pi, d_BO_pi2;
  typename AT::t_kkfloat_2d_dl d_dln_BOp_pi, d_dln_BOp_pi2;
  typename AT::t_kkfloat_2d_dl d_C1dbo, d_C2dbo, d_C3dbo;
  typename AT::t_kkfloat_2d_dl d_C1dbopi, d_C2dbopi, d_C3dbopi, d_C4dbopi;
  typename AT::t_kkfloat_2d_dl d_C1dbopi2, d_C2dbopi2, d_C3dbopi2, d_C4dbopi2;
  typename AT::t_kkfloat_2d_dl d_dDeltap_self, d_Cdbo, d_Cdbopi, d_Cdbopi2;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  DupScatterView<KK_FLOAT**, typename DAT::t_kkfloat_2d_dl::array_layout> dup_dDeltap_self;
  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_total_bo;
  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_CdDelta;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;
  NonDupScatterView<KK_FLOAT**, typename DAT::t_kkfloat_2d_dl::array_layout> ndup_dDeltap_self;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_total_bo;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_CdDelta;

  typedef Kokkos::DualView<KK_FLOAT**[7],LMPDeviceLayout,DeviceType> tdual_kkfloat_2d_n7;
  typedef typename tdual_kkfloat_2d_n7::t_dev_const_randomread t_kkfloat_2d_n7_randomread;
  typedef typename tdual_kkfloat_2d_n7::t_host t_hostkkfloat_2d_n7;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_1d d_bo_num, d_hb_num;
  typename AT::t_int_2d d_bo_list, d_hb_list;

  DAT::tdual_int_scalar k_resize_bo, k_resize_hb;
  typename AT::t_int_scalar d_resize_bo, d_resize_hb;

  typename AT::t_kkfloat_2d_dl d_sum_ovun;
  typename AT::t_kkfloat_2d_dl d_dBOp;

  int neighflag, newton_pair, maxnumneigh, maxhb, maxbo;
  int nlocal,nn,NN,eflag,vflag,acks2_flag;
  KK_FLOAT cut_nbsq, cut_hbsq, cut_bosq, bo_cut, thb_cut, thb_cutsq;
  KK_FLOAT bo_cut_bond;

  int vdwflag, lgflag;
  KK_FLOAT gp[39], p_boc1, p_boc2;

  friend void pair_virial_fdotr_compute<PairReaxFFKokkos>(PairReaxFFKokkos*);

  int bocnt,hbcnt,enobondsflag;

  typedef LR_lookup_table_kk<DeviceType> LR_lookup_table_kk_DT;

  typedef Kokkos::DualView<LR_lookup_table_kk_DT**,LMPDeviceLayout,DeviceType> tdual_LR_lookup_table_kk_2d;
  typedef typename tdual_LR_lookup_table_kk_2d::t_dev t_LR_lookup_table_kk_2d;

  tdual_LR_lookup_table_kk_2d k_LR;

  DAT::ttransform_int_2d k_tmpid;
  typename AT::t_int_2d d_tmpid;
  DAT::ttransform_kkfloat_2d k_tmpbo;
  typename AT::t_kkfloat_2d d_tmpbo;
  DAT::tdual_int_scalar k_error_flag;

  typename AT::t_int_1d d_numneigh_bonds;
  typename AT::t_tagint_2d d_neighid;
  typename AT::t_kkfloat_2d d_abo;

  typename AT::t_double_1d d_buf;
  DAT::tdual_int_scalar k_nbuf_local;

  typedef Kokkos::View<reax_int4**, LMPDeviceLayout, DeviceType> t_reax_int4_2d;

  t_reax_int4_2d d_angular_pack, d_torsion_pack;

  typename AT::t_kkfloat_2d d_angular_intermediates;

  DAT::tdual_int_1d k_count_angular_torsion;
  typename AT::t_int_1d d_count_angular_torsion;

};

template <class DeviceType>
struct PairReaxKokkosFindBondFunctor  {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;
  int groupbit;
  PairReaxFFKokkos<DeviceType> c;
  PairReaxKokkosFindBondFunctor(PairReaxFFKokkos<DeviceType>* c_ptr, int groupbit):groupbit(groupbit),c(*c_ptr){};

  KOKKOS_INLINE_FUNCTION
  void join(int &dst,
             const int &src) const {
    dst = MAX(dst,src);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &numbonds) const {
    c.calculate_find_bond_item(ii,numbonds,groupbit);
  }
};

template <class DeviceType>
struct PairReaxKokkosPackBondBufferFunctor  {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;
  PairReaxFFKokkos<DeviceType> c;
  PairReaxKokkosPackBondBufferFunctor(PairReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &j, const bool &final) const {
    c.pack_bond_buffer_item(ii,j,final);
  }
};

template <class DeviceType, bool STORE_BONDS>
struct PairReaxKokkosPackReducedBondBufferFunctor  {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef int value_type;
  PairReaxFFKokkos<DeviceType> c;
  PairReaxKokkosPackReducedBondBufferFunctor(PairReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &j, const bool &final) const {
    c.template pack_reduced_bond_buffer_item<STORE_BONDS>(ii,j,final);
  }
};

}

#endif
#endif

