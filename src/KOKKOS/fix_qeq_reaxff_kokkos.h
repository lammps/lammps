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

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/reaxff/kk,FixQEqReaxFFKokkos<LMPDeviceType>);
FixStyle(qeq/reaxff/kk/device,FixQEqReaxFFKokkos<LMPDeviceType>);
FixStyle(qeq/reaxff/kk/host,FixQEqReaxFFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_QEQ_REAXFF_KOKKOS_H
#define LMP_FIX_QEQ_REAXFF_KOKKOS_H

#include "fix_qeq_reaxff.h"

#include "kokkos_base.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "Kokkos_ScatterView.hpp"

// Include Kokkos-Kernels
#include "KokkosBlas.hpp"
#include "KokkosSparse.hpp"
#include "KokkosKernels_default_types.hpp"

namespace LAMMPS_NS {

template<class DeviceType>
class FixQEqReaxFFKokkos : public FixQEqReaxFF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixQEqReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixQEqReaxFFKokkos() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;

  // -------- Exchange --------

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  int pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &k_buf,
                           DAT::tdual_int_1d k_exchange_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space);

  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &k_indices,
                              int nrecv, int nrecv1, int nextrarecv1,
                              ExecutionSpace space);

  // -------- Forward comm --------

  int pack_forward_comm(int, int *, double *, int, int *) override;

  int pack_forward_comm_kokkos(int,
                               DAT::tdual_int_1d,
                               DAT::tdual_xfloat_1d &,
                               int, int *) override;

  void unpack_forward_comm(int, int, double *) override;

  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) override;

  // -------- Reverse comm --------

  int pack_reverse_comm(int, int, double *) override;

  int pack_reverse_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) override;

  void unpack_reverse_comm(int, int *, double *) override;
  
  void unpack_reverse_comm_kokkos(int,
                                  DAT::tdual_int_1d,
                                  DAT::tdual_xfloat_1d &) override;

 private:
  int neighflag, allocated_flag, last_allocate;
  double cutsq;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist, d_numneigh;

  DAT::tdual_float_1d k_chi_field;
  typename AT::t_float_1d d_chi_field;

  void resize_views();
  void get_chi_field();

  // -------- Mixed precision --------

  typedef float kk_compute;
  typedef double kk_accumulation;

  typedef Kokkos::DualView<kk_compute*, typename DeviceType::array_layout, DeviceType> tdual_compute_1d;
  typedef Kokkos::DualView<kk_accumulation*, typename DeviceType::array_layout, DeviceType> tdual_accumulation_1d;
  typedef typename tdual_compute_1d::t_dev t_compute_1d;
  typedef typename tdual_accumulation_1d::t_dev t_accumulation_1d;

  // -------- Extended Lagrangian --------

  tdual_compute_1d k_theta, k_theta_dot;
  t_compute_1d d_theta, d_theta_dot;
  double dt, dt_half, dt2_half, omega, omega2, K;

  void update_extended_lagrangian();

  // -------- QEQ parameters --------

  struct params_qeq {
    double chi, eta, gamma;
    KOKKOS_INLINE_FUNCTION
    params_qeq() {chi = 0; eta = 0; gamma = 0;};
    KOKKOS_INLINE_FUNCTION
    params_qeq(const params_qeq &p) { chi = p.chi; eta = p.eta; gamma = p.gamma; }
  };

  typedef Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType> tdual_qeq_1d;
  typedef typename tdual_qeq_1d::t_dev_const t_qeq_1d;
  tdual_qeq_1d k_params;
  t_qeq_1d d_params;

  DAT::tdual_float_2d k_shield;
  typename AT::t_float_2d d_shield;

  double tap0, tap1, tap2, tap3, tap4, tap5, tap6, tap7;

  // -------- CRS matrix --------

  bool matrix_sparsity_initialized, crs_matrix_allocated;

  typedef KokkosSparse::CrsMatrix<kk_compute, int, DeviceType> CRSMatrixType;
  CRSMatrixType crs_matrix;

  KOKKOS_INLINE_FUNCTION
  double calculate_H_k(const double &, const double &) const;

  void init_shielding_k();
  void build_crs_matrix();
  void update_crs_matrix_values();

  // -------- Schur CG --------

  tdual_compute_1d k_o, k_s;
  t_compute_1d d_o, d_s, d_r, d_p, d_d, d_Hdia_inv, d_b;

  void init_matvec();
  void sparse_matvec(t_compute_1d &in, t_compute_1d &out);
  void cg_solve();
  void calculate_q();

};

}    // namespace LAMMPS_NS

#endif
#endif
