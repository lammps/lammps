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
PairStyle(oxdna/excv/kk,PairOxdnaExcvKokkos<LMPDeviceType>);
PairStyle(oxdna/excv/kk/device,PairOxdnaExcvKokkos<LMPDeviceType>);
PairStyle(oxdna/excv/kk/host,PairOxdnaExcvKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA_EXCV_KOKKOS_H
#define LMP_PAIR_OXDNA_EXCV_KOKKOS_H

#include "kokkos_base.h"
#include "fix_oxdna_lrf_kokkos.h"
#include "pair_kokkos.h"
#include "pair_oxdna_excv.h"
#include "nucleotide_oxdna.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<class DeviceType>
class FixOxdnaNpairKokkos;  // forward declaration

template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos;  // forward declaration

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaExcvCompute{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct ev_tally_xyz{};

// packed per-type-pair coefficients of PairOxdnaExcvKokkos
struct ParamsOxdnaExcv2 {
  KK_FLOAT epsilon_bkbk, sigma_bkbk, cut_bkbk_ast, cutsq_bkbk_ast, lj1_bkbk, lj2_bkbk;
  KK_FLOAT b_bkbk, cut_bkbk_c, cutsq_bkbk_c, epsilon_bkbs, sigma_bkbs, cut_bkbs_ast;
  KK_FLOAT cutsq_bkbs_ast, lj1_bkbs, lj2_bkbs, b_bkbs, cut_bkbs_c, cutsq_bkbs_c;
  KK_FLOAT epsilon_bsbs, sigma_bsbs, cut_bsbs_ast, cutsq_bsbs_ast, lj1_bsbs, lj2_bsbs;
  KK_FLOAT b_bsbs, cut_bsbs_c, cutsq_bsbs_c;
};

// packed per-tetramer coefficients of PairOxdnaExcvKokkos
struct ParamsOxdnaExcv4 {
  KK_FLOAT sigma4_bsbs, cut4_bsbs_ast, cut4sq_bsbs_ast, lj14_bsbs, lj24_bsbs, b4_bsbs;
  KK_FLOAT cut4_bsbs_c, cut4sq_bsbs_c;
};

template<class DeviceType>
class PairOxdnaExcvKokkos : public PairOxdnaExcv, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdnaExcvKokkos(class LAMMPS *);
  ~PairOxdnaExcvKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  // interaction sites of the model, used by init_one() of the CPU base class
  // for the cutoffs; the kernels compute the sites themselves
  void compute_backbone_site(double e1[3], double e2[3], double e3[3],
                             double rbk[3]) const override
  {
    if ((oxdnaflag == OXDNA2) || (oxdnaflag == OXDNA3)) {
      NucleotideOxdna2 oxdna2;
      oxdna2.backbone_site(e1, e2, nullptr, rbk);
    } else if (oxdnaflag == OXRNA2) {
      NucleotideOxrna2 oxrna2;
      oxrna2.backbone_site(e1, nullptr, e3, rbk);
    } else {
      NucleotideOxdna1 oxdna1;
      oxdna1.backbone_site(e1, nullptr, nullptr, rbk);
    }
  }
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
  void coeff(int, char **) override;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:

  void coeff_set_tetramers_kokkos(int narg, char **arg);

  int oxdnaflag;
  enum EnabledOXDNAFlag{OXDNA=1,OXDNA2=2,OXRNA2=4,OXDNA3=8};

  t_oxdna_packed_sub<DeviceType> x;    // positions in the packed record of fix OXDNA/LRF/kk
  t_oxdna_packed_col<DeviceType> xn_type;    // atom types in the packed record
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;

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

  // s=sugar-phosphate backbone site, b=base site, st=stacking site
  // excluded volume interaction parameters
  // all per-type-pair coefficients of a pair packed in one struct
  Kokkos::DualView<ParamsOxdnaExcv2 **, Kokkos::LayoutRight, DeviceType> k_params2_excv;
  typename Kokkos::DualView<ParamsOxdnaExcv2 **, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params2_excv;
  // tetramer-dependent coefficients
  // all per-tetramer coefficients of a pair packed in one struct
  Kokkos::DualView<ParamsOxdnaExcv4 ****, Kokkos::LayoutRight, DeviceType> k_params4_excv;
  typename Kokkos::DualView<ParamsOxdnaExcv4 ****, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params4_excv;
  // per-atom arrays for local unit vectors
  t_oxdna_packed_sub<DeviceType> d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread tag;
  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread id5p;
  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread id3p;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  // local index of the atom with a given tag (atom->map()), -1 if the tag is -1 or not present
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int map_tag(const tagint &itag) const;

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

  friend void pair_virial_fdotr_compute<PairOxdnaExcvKokkos>(PairOxdnaExcvKokkos*);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
  FixOxdnaNpairKokkos<DeviceType> *fix_oxdna_npairKK;    // ptr to OXDNA/NPAIR/kk fix
  FixOxdnaPrimeNeighsKokkos<DeviceType> *fix_oxdna_prime_neighsKK;    // ptr to OXDNA/PRIME_NEIGHS/kk fix
};

}

#endif
#endif

