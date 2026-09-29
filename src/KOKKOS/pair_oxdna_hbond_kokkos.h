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
PairStyle(oxdna/hbond/kk,PairOxdnaHbondKokkos<LMPDeviceType>);
PairStyle(oxdna/hbond/kk/device,PairOxdnaHbondKokkos<LMPDeviceType>);
PairStyle(oxdna/hbond/kk/host,PairOxdnaHbondKokkos<LMPHostType>);
PairStyle(oxdna2/hbond/kk,PairOxdnaHbondKokkos<LMPDeviceType>);
PairStyle(oxdna2/hbond/kk/device,PairOxdnaHbondKokkos<LMPDeviceType>);
PairStyle(oxdna2/hbond/kk/host,PairOxdnaHbondKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA_HBOND_KOKKOS_H
#define LMP_PAIR_OXDNA_HBOND_KOKKOS_H

#include "kokkos_base.h"
#include "fix_oxdna_lrf_kokkos.h"
#include "pair_kokkos.h"
#include "pair_oxdna_hbond.h"
#include "nucleotide_oxdna.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<class DeviceType>
class FixOxdnaNpairKokkos;  // forward declaration

template<class DeviceType>
class PairOxdna3XstkKokkos;  // forward declaration

template<class DeviceType, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct PairOxdna3HbXstkFused;    // fused hbond + oxdna3/xstk kernel

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaHbondCompute{};

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaHbondComputeGPUPair{};

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaHbondComputeGPURadial{};

// packed per-type-pair coefficients of PairOxdnaHbondKokkos
struct ParamsHbond {
  KK_FLOAT epsilon_hb, a_hb, cut_hb_0, cut_hb_c, cut_hb_lo, cut_hb_hi;
  KK_FLOAT cut_hb_lc, cut_hb_hc, b_hb_lo, b_hb_hi, shift_hb, cutsq_hb_hc;
  KK_FLOAT a_hb1, theta_hb1_0, dtheta_hb1_ast, b_hb1, dtheta_hb1_c, a_hb2;
  KK_FLOAT theta_hb2_0, dtheta_hb2_ast, b_hb2, dtheta_hb2_c, a_hb3, theta_hb3_0;
  KK_FLOAT dtheta_hb3_ast, b_hb3, dtheta_hb3_c, a_hb4, theta_hb4_0, dtheta_hb4_ast;
  KK_FLOAT b_hb4, dtheta_hb4_c, a_hb7, theta_hb7_0, dtheta_hb7_ast, b_hb7;
  KK_FLOAT dtheta_hb7_c, a_hb8, theta_hb8_0, dtheta_hb8_ast, b_hb8, dtheta_hb8_c;
};

template<class DeviceType>
class PairOxdnaHbondKokkos : public PairOxdnaHbond, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdnaHbondKokkos(class LAMMPS *);
  ~PairOxdnaHbondKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  // interaction sites of the model, used by init_one() of the CPU base class
  // for the cutoffs; the kernels compute the sites themselves
  void compute_base_site(int type, double e1[3], double /*e2*/[3], double /*e3*/[3],
                         double rbs[3]) const override
  {
    if (oxdnaflag == OXDNA3) {
      NucleotideOxdna3 oxdna3;
      switch (type) {
        case 0:
          oxdna3.base_site<0>(e1, nullptr, nullptr, rbs);
          break;
        case 1:
          oxdna3.base_site<1>(e1, nullptr, nullptr, rbs);
          break;
        case 2:
          oxdna3.base_site<2>(e1, nullptr, nullptr, rbs);
          break;
        case 3:
          oxdna3.base_site<3>(e1, nullptr, nullptr, rbs);
          break;
      }
    } else {
      NucleotideOxdna1 oxdna1;
      oxdna1.base_site<0>(e1, nullptr, nullptr, rbs);
    }
  }

  // Standard non-GPU Compute Functor(s). 1 with EV_FLOAT, 1 without.

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaHbondCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaHbondCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

// GPU ComputeGPUPair Functor(s). 1 with EV_FLOAT, 1 without.

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int RADIAL_ONLY = 0>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool screened_pair_body(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ipair,
                          KK_ACC_FLOAT (&fa)[3], KK_ACC_FLOAT (&ta)[3], EV_FLOAT &ev) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  // first phase of the two-phase evaluation (OXDNA_KK_TWO_PHASE)
  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaHbondComputeGPURadial<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int PAIRWISE = 0>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  // fused hbond + oxdna3/xstk kernel (OXDNA_KK_FUSE_HBXSTK), set up by
  // pair oxdna3/xstk/kk in its init_style()
  PairOxdna3XstkKokkos<DeviceType> *fuse_partner;
  bigint fuse_ncompute;    // # of compute() calls that used the fused kernel
  bool fuse_supported() const;

 protected:

  int oxdnaflag;
  enum EnabledOXDNAFlag{OXDNA=1,OXDNA3=2};

  t_oxdna_packed_sub<DeviceType> x;    // positions in the packed record of fix OXDNA/LRF/kk
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d_randomread tag;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;
  double special_lj[4];

  typename AT::tdual_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int neighflag;
  int nlocal, eflag, vflag;
  int anum;

  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;
  // GPU-specific: screened neighbor arrays for npair fix
  DAT::tdual_uint64_1d k_pairs_screened;
  typename AT::t_uint64_1d d_pairs_screened;
  typename AT::t_int_1d d_screened_offsets;  // per-atom segments of d_pairs_screened
  int screened_launch_count;   // number of threads of the screened-pair kernels
  int screened_pair_count;
  // compact list of the screened pairs that pass the radial test (OXDNA_KK_TWO_PHASE)
  typename AT::t_int_1d d_radial_pairs;
  typename AT::t_int_scalar d_radial_count;

  DAT::tdual_int_1d k_idc;
  typename AT::t_int_1d_randomread d_idc;
  int unique_basepair_enabled;
  bigint last_idc_nbuild;
  int last_idc_nall;

  // hydrogen-bonding interaction parameters
  // all per-type-pair coefficients of a pair packed in one struct
  Kokkos::DualView<ParamsHbond **, Kokkos::LayoutRight, DeviceType> k_params_hb;
  typename Kokkos::DualView<ParamsHbond **, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params_hb;
  // per-atom arrays for local unit vectors
  t_oxdna_packed_sub<DeviceType> d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> dup_torque;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> ndup_torque;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  void allocate() override;

  friend void pair_virial_fdotr_compute<PairOxdnaHbondKokkos>(PairOxdnaHbondKokkos*);
  template<class, int, int, int> friend struct PairOxdna3HbXstkFused;
  template<class> friend class PairOxdna3XstkKokkos;

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
  FixOxdnaNpairKokkos<DeviceType> *fix_oxdna_npairKK;    // ptr to OXDNA/NPAIR/kk fix

 private:

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_radial_terms(const int &atype, const int &btype, const KK_FLOAT &r_hb,
    KK_FLOAT &f1, KK_FLOAT &df1) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta1_terms(const int &atype, const int &btype,
    const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
    KK_FLOAT &theta1, KK_FLOAT &f4t1, KK_FLOAT &df4t1) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta2_terms(const int &atype, const int &btype,
    const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
    KK_FLOAT &theta2, KK_FLOAT &cost2, KK_FLOAT &f4t2, KK_FLOAT &df4t2) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta3_terms(const int &atype, const int &btype,
    const KK_FLOAT (&b_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
    KK_FLOAT &theta3, KK_FLOAT &cost3, KK_FLOAT &f4t3, KK_FLOAT &df4t3) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta4_terms(const int &atype, const int &btype,
    const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
    KK_FLOAT &theta4, KK_FLOAT &f4t4, KK_FLOAT &df4t4) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta7_terms(const int &atype, const int &btype,
    const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
    KK_FLOAT &theta7, KK_FLOAT &cost7, KK_FLOAT &f4t7, KK_FLOAT &df4t7) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  bool hbond_theta8_terms(const int &atype, const int &btype,
    const KK_FLOAT (&b_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
    KK_FLOAT &theta8, KK_FLOAT &cost8, KK_FLOAT &f4t8, KK_FLOAT &df4t8) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void hbond_force_contrib(const KK_FLOAT &f1,
    const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
    const KK_FLOAT &f4t4, const KK_FLOAT &f4t7, const KK_FLOAT &f4t8,
    const KK_FLOAT &df1, const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
    const KK_FLOAT &df4t7, const KK_FLOAT &df4t8,
    const KK_FLOAT &rinv_hb, const KK_FLOAT &factor_lj,
    const KK_FLOAT &theta2, const KK_FLOAT &theta3, const KK_FLOAT &theta7, const KK_FLOAT &theta8,
    const KK_FLOAT &cost2, const KK_FLOAT &cost3, const KK_FLOAT &cost7, const KK_FLOAT &cost8,
    const KK_FLOAT (&delr_hb)[3], const KK_FLOAT (&delr_hb_norm)[3],
    const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
    const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
    const KK_FLOAT (&ra_chb)[3], const KK_FLOAT (&rb_chb)[3],
    KK_ACC_FLOAT (&delf)[3], KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void hbond_torque_contrib(const KK_FLOAT &f1,
    const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
    const KK_FLOAT &f4t4, const KK_FLOAT &f4t7, const KK_FLOAT &f4t8,
    const KK_FLOAT &df4t1, const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
    const KK_FLOAT &df4t4, const KK_FLOAT &df4t7, const KK_FLOAT &df4t8,
    const KK_FLOAT &factor_lj,
    const KK_FLOAT &theta1, const KK_FLOAT &theta2, const KK_FLOAT &theta3,
    const KK_FLOAT &theta4, const KK_FLOAT &theta7, const KK_FLOAT &theta8,
    const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
    const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
    const KK_FLOAT (&delr_hb_norm)[3],
    KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const;
};

}

#endif
#endif

