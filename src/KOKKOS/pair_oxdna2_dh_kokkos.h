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
PairStyle(oxdna2/dh/kk,PairOxdna2DhKokkos<LMPDeviceType>);
PairStyle(oxdna2/dh/kk/device,PairOxdna2DhKokkos<LMPDeviceType>);
PairStyle(oxdna2/dh/kk/host,PairOxdna2DhKokkos<LMPHostType>);
PairStyle(oxdna3/dh/kk,PairOxdna2DhKokkos<LMPDeviceType>);
PairStyle(oxdna3/dh/kk/device,PairOxdna2DhKokkos<LMPDeviceType>);
PairStyle(oxdna3/dh/kk/host,PairOxdna2DhKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA2_DH_KOKKOS_H
#define LMP_PAIR_OXDNA2_DH_KOKKOS_H

#include "kokkos_base.h"
#include "fix_oxdna_lrf_kokkos.h"
#include "pair_kokkos.h"
#include "pair_oxdna2_dh.h"
#include "nucleotide_oxdna.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdna2DhCompute{};

// packed per-type-pair coefficients of PairOxdna2DhKokkos
struct ParamsOxdnaDh {
  KK_FLOAT qeff_dh_pf, kappa_dh, b_dh, cut_dh_ast, cutsq_dh_ast, cut_dh_c, cutsq_dh_c;
};

template<class DeviceType>
class PairOxdna2DhKokkos : public PairOxdna2Dh, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdna2DhKokkos(class LAMMPS *);
  ~PairOxdna2DhKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  // interaction sites of the model, used by init_one() of the CPU base class
  // for the cutoffs; the kernels compute the sites themselves
  void compute_backbone_site(double e1[3], double e2[3], double e3[3],
                             double rbk[3]) const override
  {
    if (oxdnaflag == OXRNA2) {
      // same as PairOxrna2Dh::compute_backbone_site()
      const double dx_cbk = ConstantsOxdna::get_dx_cbk_oxdna1();
      const double dz_cbk = ConstantsOxdna::get_dz_cbk_oxrna2();
      rbk[0] = dx_cbk * e1[0] + dz_cbk * e3[0];
      rbk[1] = dx_cbk * e1[1] + dz_cbk * e3[1];
      rbk[2] = dx_cbk * e1[2] + dz_cbk * e3[2];
    } else {
      NucleotideOxdna2 oxdna2;
      oxdna2.backbone_site(e1, e2, nullptr, rbk);
    }
  }

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdna2DhCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdna2DhCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:

  int oxdnaflag;
  enum EnabledOXDNAFlag{OXDNA2=1,OXRNA2=2};

  t_oxdna_packed_sub<DeviceType> x;    // positions in the packed record of fix OXDNA/LRF/kk
  t_oxdna_packed_col<DeviceType> xn_type;    // atom types in the packed record
  t_oxdna_packed_col<DeviceType> xn_qeff;    // effective charges in the packed record
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread qeff;

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
  int nsplit;    // threads per atom of the compute kernel

  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;

  // debye huckel interaction parameters
  // all coefficients of a type pair packed in one struct
  Kokkos::DualView<ParamsOxdnaDh **, Kokkos::LayoutRight, DeviceType> k_params_dh;
  typename Kokkos::DualView<ParamsOxdnaDh **, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params_dh;
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

  friend void pair_virial_fdotr_compute<PairOxdna2DhKokkos>(PairOxdna2DhKokkos*);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
};

}

#endif
#endif

