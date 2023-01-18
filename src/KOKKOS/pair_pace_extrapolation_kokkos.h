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
PairStyle(pace/extrapolation/kk,PairPACEExtrapolationKokkos<LMPDeviceType>);
PairStyle(pace/extrapolation/kk/device,PairPACEExtrapolationKokkos<LMPDeviceType>);
PairStyle(pace/extrapolation/kk/host,PairPACEExtrapolationKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_PACE_EXTRAPOLATION_KOKKOS_H
#define LMP_PAIR_PACE_EXTRAPOLATION_KOKKOS_H

#include "pair_pace_extrapolation.h"
#include "kokkos_type.h"
#include "pair_kokkos.h"

class SplineInterpolator;

namespace LAMMPS_NS {

template<class DeviceType>
class PairPACEExtrapolationKokkos : public PairPACEExtrapolation {
 public:
  struct TagPairPACEComputeNeigh{};
  struct TagPairPACEComputeRadial{};
  struct TagPairPACEComputeYlm{};
  struct TagPairPACEComputeAi{};
  struct TagPairPACEConjugateAi{};
  struct TagPairPACEComputeRho{};
  struct TagPairPACEComputeFS{};
  struct TagPairPACEComputeGamma{};
  struct TagPairPACEComputeWeights{};
  struct TagPairPACEComputeDerivative{};

  template<int NEIGHFLAG, int EVFLAG>
  struct TagPairPACEComputeForce{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  using complex = SNAComplex<double>;

  PairPACEExtrapolationKokkos(class LAMMPS *);
  ~PairPACEExtrapolationKokkos() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRadial,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeRadial>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeYlm,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeYlm>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeAi,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEConjugateAi,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeRho,const int& iter) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeFS,const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeGamma, const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeWeights,const int& iter) const;

  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeDerivative,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeDerivative>::member_type& team) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii, EV_FLOAT&) const;


  void *extract(const char *str, int &dim);
  void *extract_peratom(const char *str, int &ncol);

 protected:
  int inum, maxneigh, chunk_size, chunk_offset, idx_ms_combs_max, total_num_functions_max;
  int host_flag;
  int gamma_flag;

  int eflag, vflag;

  int neighflag, max_ndensity;
  int nelements, lmax, nradmax, nradbase;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  typedef Kokkos::DualView<F_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq, k_scale;
  typedef Kokkos::View<F_FLOAT**, DeviceType> t_fparams;
  t_fparams d_cutsq, d_scale;

  typename AT::t_int_1d d_map;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> dup_f;
  DupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> dup_vatom;

  NonDupScatterView<F_FLOAT*[3], typename DAT::t_f_array::array_layout> ndup_f;
  NonDupScatterView<F_FLOAT*[6], typename DAT::t_virial_array::array_layout> ndup_vatom;

  friend void pair_virial_fdotr_compute<PairPACEExtrapolationKokkos>(PairPACEExtrapolationKokkos*);

  void grow(int, int);
  void copy_pertype();
  void copy_splines();
  void copy_tilde();
  void allocate() override;
  void precompute_harmonics();
  double memory_usage() override;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  void compute_barplm(int, int, double, int) const;

  KOKKOS_INLINE_FUNCTION
  void compute_ylm(int, int, double, double, double, int) const;

  KOKKOS_INLINE_FUNCTION
  void cutoff_func_poly(const double, const double, const double, double &, double &) const;

  KOKKOS_INLINE_FUNCTION
  void Fexp(const double, const double, double &, double &) const;

  KOKKOS_INLINE_FUNCTION
  void FexpShiftedScaled(const double, const double, double &, double &) const;

  KOKKOS_INLINE_FUNCTION
  void inner_cutoff(const double, const double, const double, double &, double &) const;

  KOKKOS_INLINE_FUNCTION
  void FS_values_and_derivatives(const int, double&, const int) const;

  KOKKOS_INLINE_FUNCTION
  void evaluate_splines(const int, const int, double, int, int, int, int) const;

  template<class TagStyle>
  void check_team_size_for(int, int&, int);

  template<class TagStyle>
  void check_team_size_reduce(int, int&, int);

  // Utility routine which wraps computing per-team scratch size requirements for
  // ComputeNeigh, ComputeUi, and ComputeFusedDeidrj
  template <typename scratch_type>
  int scratch_size_helper(int values_per_team);

  typedef Kokkos::View<int*, DeviceType> t_ace_1i;
  typedef Kokkos::View<int**, DeviceType> t_ace_2i;
  typedef Kokkos::View<int***, DeviceType> t_ace_3i;
  typedef Kokkos::View<int****, DeviceType> t_ace_4i;
  typedef Kokkos::View<double*, DeviceType> t_ace_1d;
  typedef Kokkos::View<double**, DeviceType> t_ace_2d;
  typedef Kokkos::View<double*[3], DeviceType> t_ace_2d3;
  typedef Kokkos::View<double***, DeviceType> t_ace_3d;
  typedef Kokkos::View<const double***, DeviceType> tc_ace_3d;
  typedef Kokkos::View<double**[3], DeviceType> t_ace_3d3;
  typedef Kokkos::View<double**[4], DeviceType> t_ace_3d4;
  typedef Kokkos::View<double****, DeviceType> t_ace_4d;
  typedef Kokkos::View<complex*, DeviceType> t_ace_1c;
  typedef Kokkos::View<complex**, DeviceType> t_ace_2c;
  typedef Kokkos::View<complex***, DeviceType> t_ace_3c;
  typedef Kokkos::View<complex**[3], DeviceType> t_ace_3c3;
  typedef Kokkos::View<complex****, DeviceType> t_ace_4c;
  typedef Kokkos::View<complex***[3], DeviceType> t_ace_4c3;

  typedef Kokkos::View<double*>::HostMirror th_ace_1d;

  t_ace_3d A_rank1;
  t_ace_4c A;

  t_ace_3c A_list;
  t_ace_3c A_forward_prod;

  t_ace_3d weights_rank1;
  t_ace_4c weights;

  t_ace_1d e_atom;
  t_ace_2d rhos;
  t_ace_2d dF_drho;

  // hard-core repulsion
  t_ace_1d rho_core;
  t_ace_3c dB_flatten;
  t_ace_2d cr;
  t_ace_2d dcr;
  t_ace_1d dF_drho_core;

  // radial functions
  t_ace_4d fr;
  t_ace_4d dfr;
  t_ace_3d gr;
  t_ace_3d dgr;
  t_ace_3d d_values;
  t_ace_3d d_derivatives;

  // inverted active set
  tc_ace_3d d_ASI;
  t_ace_2d projections;
  t_ace_1d d_gamma;
  th_ace_1d h_gamma;

  // Spherical Harmonics
  void pre_compute_harmonics(int);

  KOKKOS_INLINE_FUNCTION
  void compute_barplm(double rz, int lmaxi);

  KOKKOS_INLINE_FUNCTION
  void compute_ylm(double rx, double ry, double rz, int lmaxi);

  t_ace_1d alm;
  t_ace_1d blm;
  t_ace_1d cl;
  t_ace_1d dl;

  t_ace_3d plm;
  t_ace_3d dplm;

  t_ace_3c ylm;
  t_ace_4c3 dylm;

  // short neigh list
  t_ace_1i d_ncount;
  t_ace_2d d_mu;
  t_ace_2d d_rnorms;
  t_ace_3d3 d_rhats;
  t_ace_2i d_nearest;

  // per-type
  t_ace_1i d_ndensity;
  t_ace_1i d_npoti;
  t_ace_1d d_rho_core_cutoff;
  t_ace_1d d_drho_core_cutoff;
  t_ace_1d d_E0vals;
  t_ace_2d d_wpre;
  t_ace_2d d_mexp;

  // tilde
  t_ace_1i d_idx_ms_combs_count;
  t_ace_1i d_total_basis_size;
  t_ace_2i d_rank;
  t_ace_2i d_num_ms_combs;
  t_ace_2i d_func_inds;
  t_ace_3i d_mus;
  t_ace_3i d_ns;
  t_ace_3i d_ls;
  t_ace_3i d_ms_combs;
//  t_ace_3d d_ctildes;
  t_ace_2d d_gen_cgs;
  t_ace_3d d_coeffs;

  t_ace_3d3 f_ij;

 public:
  struct SplineInterpolatorKokkos {
    int ntot, nlut, num_of_functions;
    double cutoff, deltaSplineBins, invrscalelookup, rscalelookup;

    t_ace_3d4 lookupTable;

    void operator=(const SplineInterpolator &spline);

    void deallocate() {
      lookupTable = t_ace_3d4();
    }

    double memory_usage() {
      return lookupTable.span() * sizeof(typename decltype(lookupTable)::value_type);
    }

    KOKKOS_INLINE_FUNCTION
    void calcSplines(const int ii, const int jj, const double r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const;
  };

  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_gk;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_rnl;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_hc;

};
}    // namespace LAMMPS_NS

#endif
#endif
