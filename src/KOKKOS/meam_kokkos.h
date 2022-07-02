/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MEAMKOKKOS_H
#define LMP_MEAMKOKKOS_H

#include "kokkos.h"
#include "meam.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include <cmath>
#include <cstdlib>

namespace LAMMPS_NS {

struct TagMEAMDensFinal {};
template <int NEIGHFLAG> struct TagMEAMDensInit {
};
struct TagMEAMZero {};
template <int NEIGHFLAG> struct TagMEAMForce {
};

template <class DeviceType> class MEAMKokkos : public MEAM {
 public:
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  MEAMKokkos(Memory *mem);
  ~MEAMKokkos() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMDensFinal, const int &, EV_FLOAT &) const;

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void operator()(TagMEAMDensInit<NEIGHFLAG>, const int &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagMEAMZero, const int &) const;

  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void operator()(TagMEAMForce<NEIGHFLAG>, const int &, EV_FLOAT &) const;

 private:
  // parameters to meam_dens_init

  int ntype, nlocal;
  typename AT::t_int_1d type;
  typename AT::t_int_1d d_offset;
  typename AT::t_int_1d d_map;
  typename AT::t_int_2d d_scale;
  typename AT::t_x_array x;
  typename AT::t_int_1d d_numneigh_half;
  typename AT::t_int_1d d_numneigh_full;
  typename AT::t_neighbors_2d d_neighbors_half;
  typename AT::t_neighbors_2d d_neighbors_full;
  typename AT::t_int_1d d_ilist_half;
  typename AT::t_f_array f;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  // parameters to meam_dens_final

  typename AT::t_int_scalar d_errorflag;
  int eflag_either, eflag_global, eflag_atom, vflag_either, vflag_global, vflag_atom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;

 public:
  void meam_dens_setup(int, int, int) override;
  void meam_setup_done(double *) override;
  void meam_dens_init(int, int, typename AT::t_int_1d, typename AT::t_int_1d,
                      typename AT::t_x_array, typename AT::t_int_1d, typename AT::t_int_1d,
                      typename AT::t_int_1d, typename AT::t_neighbors_2d,
                      typename AT::t_neighbors_2d, typename AT::t_int_1d, int, int);
  void meam_dens_final(int, int, int, int, typename ArrayTypes<DeviceType>::t_efloat_1d, int,
                       typename AT::t_int_1d, typename AT::t_int_1d, typename AT::t_int_2d, int &,
                       EV_FLOAT &);
  void meam_force(int, int, int, int, int, typename ArrayTypes<DeviceType>::t_efloat_1d, int,
                  typename AT::t_int_1d, typename AT::t_int_1d, typename AT::t_x_array,
                  typename AT::t_int_1d, typename AT::t_int_1d, typename AT::t_f_array,
                  typename ArrayTypes<DeviceType>::t_virial_array, typename AT::t_int_1d,
                  typename AT::t_int_1d, typename AT::t_neighbors_2d, typename AT::t_neighbors_2d,
                  int, int, EV_FLOAT &);
  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void getscreen(int, int, typename AT::t_x_array, typename AT::t_int_1d,
                                        typename AT::t_int_1d, int, typename AT::t_int_1d,
                                        typename AT::t_int_1d) const;
  template <int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION void calc_rho1(int, int, typename AT::t_int_1d, typename AT::t_int_1d,
                                        typename AT::t_x_array, typename AT::t_int_1d, int) const;
  KOKKOS_INLINE_FUNCTION
  double fcut(const double xi) const;
  KOKKOS_INLINE_FUNCTION
  double dfcut(const double xi, double &dfc) const;
  KOKKOS_INLINE_FUNCTION
  double dCfunc(const double, const double, const double) const;
  KOKKOS_INLINE_FUNCTION
  void dCfunc2(const double, const double, const double, double &, double &) const;
  KOKKOS_INLINE_FUNCTION
  double G_gam(const double, const int, int &) const;
  KOKKOS_INLINE_FUNCTION
  double dG_gam(const double, const int, double &) const;
  KOKKOS_INLINE_FUNCTION
  double zbl(const double, const int, const int) const;
  KOKKOS_INLINE_FUNCTION
  double embedding(const double, const double, const double, double &) const;
  KOKKOS_INLINE_FUNCTION
  double erose(const double, const double, const double, const double, const double, const double,
               const int) const;
  KOKKOS_INLINE_FUNCTION
  void get_shpfcn(const lattice_t latt, const double sthe, const double cthe, double (&s)[3]) const;
  KOKKOS_INLINE_FUNCTION
  int get_Zij(const lattice_t) const;

 public:
  DAT::tdual_ffloat_1d k_rho, k_rho0, k_rho1, k_rho2, k_rho3, k_frhop;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_rho, d_rho0, d_rho1, d_rho2, d_rho3, d_frhop;
  HAT::t_ffloat_1d h_rho, h_rho0, h_rho1, h_rho2, h_rho3, h_frhop;
  DAT::tdual_ffloat_1d k_gamma, k_dgamma1, k_dgamma2, k_dgamma3, k_arho2b;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_gamma, d_dgamma1, d_dgamma2, d_dgamma3, d_arho2b;
  HAT::t_ffloat_1d h_gamma, h_dgamma1, h_dgamma2, h_dgamma3, h_arho2b;
  DAT::tdual_ffloat_2d k_arho1, k_arho2, k_arho3, k_arho3b, k_t_ave, k_tsq_ave;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_arho1, d_arho2, d_arho3, d_arho3b, d_t_ave,
      d_tsq_ave;
  HAT::t_ffloat_2d h_arho1, h_arho2, h_arho3, h_arho3b, h_t_ave, h_tsq_ave;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_phir, d_phirar, d_phirar1, d_phirar2, d_phirar3,
      d_phirar4, d_phirar5, d_phirar6;
  DAT::tdual_ffloat_1d k_scrfcn, k_dscrfcn, k_fcpair;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_scrfcn, d_dscrfcn, d_fcpair;
  HAT::t_ffloat_1d h_scrfcn, h_dscrfcn, h_fcpair;

 protected:
  int need_dup;
  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template <typename DataType, typename Layout>
  using DupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template <typename DataType, typename Layout>
  using NonDupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<typename decltype(d_rho0)::data_type, typename decltype(d_rho0)::array_layout>
      dup_rho0;
  NonDupScatterView<typename decltype(d_rho0)::data_type, typename decltype(d_rho0)::array_layout>
      ndup_rho0;
  DupScatterView<typename decltype(d_arho2b)::data_type, typename decltype(d_arho2b)::array_layout>
      dup_arho2b;
  NonDupScatterView<typename decltype(d_arho2b)::data_type,
                    typename decltype(d_arho2b)::array_layout>
      ndup_arho2b;
  DupScatterView<typename decltype(d_arho1)::data_type, typename decltype(d_arho1)::array_layout>
      dup_arho1;
  NonDupScatterView<typename decltype(d_arho1)::data_type, typename decltype(d_arho1)::array_layout>
      ndup_arho1;
  DupScatterView<typename decltype(d_arho2)::data_type, typename decltype(d_arho2)::array_layout>
      dup_arho2;
  NonDupScatterView<typename decltype(d_arho2)::data_type, typename decltype(d_arho2)::array_layout>
      ndup_arho2;
  DupScatterView<typename decltype(d_arho3)::data_type, typename decltype(d_arho3)::array_layout>
      dup_arho3;
  NonDupScatterView<typename decltype(d_arho3)::data_type, typename decltype(d_arho3)::array_layout>
      ndup_arho3;
  DupScatterView<typename decltype(d_arho3b)::data_type, typename decltype(d_arho3b)::array_layout>
      dup_arho3b;
  NonDupScatterView<typename decltype(d_arho3b)::data_type,
                    typename decltype(d_arho3b)::array_layout>
      ndup_arho3b;
  DupScatterView<typename decltype(d_t_ave)::data_type, typename decltype(d_t_ave)::array_layout>
      dup_t_ave;
  NonDupScatterView<typename decltype(d_t_ave)::data_type, typename decltype(d_t_ave)::array_layout>
      ndup_t_ave;
  DupScatterView<typename decltype(d_tsq_ave)::data_type,
                 typename decltype(d_tsq_ave)::array_layout>
      dup_tsq_ave;
  NonDupScatterView<typename decltype(d_tsq_ave)::data_type,
                    typename decltype(d_tsq_ave)::array_layout>
      ndup_tsq_ave;
  DupScatterView<typename decltype(f)::data_type, typename decltype(f)::array_layout> dup_f;
  NonDupScatterView<typename decltype(f)::data_type, typename decltype(f)::array_layout> ndup_f;
  DupScatterView<typename decltype(d_eatom)::data_type, typename decltype(d_eatom)::array_layout>
      dup_eatom;
  NonDupScatterView<typename decltype(d_eatom)::data_type, typename decltype(d_eatom)::array_layout>
      ndup_eatom;
  DupScatterView<typename decltype(d_vatom)::data_type, typename decltype(d_vatom)::array_layout>
      dup_vatom;
  NonDupScatterView<typename decltype(d_vatom)::data_type, typename decltype(d_vatom)::array_layout>
      ndup_vatom;
};

KOKKOS_INLINE_FUNCTION
static bool iszero_kk(const double f)
{
  return fabs(f) < 1e-20;
}

KOKKOS_INLINE_FUNCTION
static bool isone_kk(const double f)
{
  return fabs(f - 1.0) < 1e-20;
}

KOKKOS_INLINE_FUNCTION
static double fdiv_zero_kk(const double n, const double d)
{
  if (iszero_kk(d)) return 0.0;
  return n / d;
}

// Functions we need for compat

}    // namespace LAMMPS_NS
#include "meam_impl_kokkos.h"

#endif
