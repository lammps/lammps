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
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "Kokkos_ScatterView.hpp"

// Include Kokkos-Kernels
#include "KokkosBlas.hpp"
#include "KokkosSparse.hpp"
#include "KokkosKernels_default_types.hpp"

namespace LAMMPS_NS {

struct TagQEqZero{};
struct TagQEqInitMatvec{};
struct TagQEqUpdate{};
struct TagQEqCalculateQ{};
struct TagQEqPackForwardComm{};
struct TagQEqUnpackForwardComm{};
struct TagQEqPackExchange{};
struct TagQEqUnpackExchange{};

template<class DeviceType>
class FixQEqReaxFFKokkos : public FixQEqReaxFF {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typedef KeyValueType<int, value_type, Kokkos::DefaultBinOp<int> > KeyValueType;
  typedef BinSort<KeyValueType, Kokkos::DefaultBinOp<int> > BinOp;

  // Mixed precision support
  // Use single precision for matrix and vectors, double precision for critical operations
  typedef float compute_float;      
  typedef double precision_float;   

  // Define views for different precision types
  typedef Kokkos::View<compute_float*, typename AT::t_float_1d::array_layout, 
          typename KKDevice<DeviceType>::value> t_compute_1d;
  typedef Kokkos::View<precision_float*, typename AT::t_float_1d::array_layout, 
          typename KKDevice<DeviceType>::value> t_precision_1d;
  
  // Use CSR matrix type for better performance
  typedef KokkosSparse::CrsMatrix<compute_float, int, DeviceType> CSRMatrixType;

  FixQEqReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixQEqReaxFFKokkos() override;
  void cleanup_copy();
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqInitMatvec, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqUpdate, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqCalculateQ, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqPackExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqUnpackExchange, const int&) const;
  
  struct params_qeq {
    KOKKOS_INLINE_FUNCTION
    params_qeq() {chi = 0; eta = 0; gamma = 0;};
    KOKKOS_INLINE_FUNCTION
    params_qeq(const params_qeq &p) {
      chi = p.chi;
      eta = p.eta;
      gamma = p.gamma;
    }
    double chi, eta, gamma;
  };

  struct FixQEqReaxFFKokkosNumNeighFunctor {
    typedef DeviceType device_type;

    FixQEqReaxFFKokkos<DeviceType> *fix;

    FixQEqReaxFFKokkosNumNeighFunctor(FixQEqReaxFFKokkos<DeviceType> *_fix) {
      fix = _fix;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int &ii, bigint &totneigh) const {
      fix->num_neigh_item(ii, totneigh);
    }

    KOKKOS_INLINE_FUNCTION
    void init(bigint &dst) const { dst = 0; }

    KOKKOS_INLINE_FUNCTION
    void join(volatile bigint &dst, const volatile bigint &src) const { dst += src; }
  };

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyValueType, BinOp> &Sorter);
  int pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &k_buf,
                          DAT::tdual_int_1d k_exchange_sendlist,
                          DAT::tdual_int_1d k_copylist,
                          ExecutionSpace space);
  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                             DAT::tdual_int_1d &k_indices, int nrecv,
                             int nrecv1, int nextrarecv1,
                             ExecutionSpace space);
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d,
                              DAT::tdual_xfloat_1d &,
                              int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int nn, NN, nlocal_last_allocate;
  int pack_flag;
  double last_allocate;

  // Extended Lagrangian variables
  t_compute_1d d_theta;      // Auxiliary charges (single precision)
  t_compute_1d d_theta_dot;  // Auxiliary charge velocities (single precision)
  double omega;             // Frequency parameter
  double K;                 // K = omega²δt²

  DAT::tdual_float_1d k_theta;
  DAT::tdual_float_1d k_theta_dot;

  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_f_array f;
  typename AT::t_int_1d mask;
  typename AT::t_int_1d type;
  typename AT::t_tagint_1d tag;

  DAT::tdual_float_1d k_q;
  typename AT::t_float_1d d_q;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  // Parameters for ReaxFF
  typedef Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType> tdual_params_qeq_1d;
  typedef typename tdual_params_qeq_1d::t_dev t_params_qeq_1d;
  t_params_qeq_1d params;
  DAT::tdual_float_2d k_shield;
  typename AT::t_float_2d d_shield;
  DAT::tdual_float_1d k_tap;
  typename AT::t_float_1d d_tap;

  // QEq vectors - use single precision for intermediate calculations
  DAT::tdual_float_1d k_o;
  t_compute_1d d_o;
  DAT::t_float_1d h_o;
  
  t_compute_1d d_r;         // Residual vector
  t_compute_1d d_p;         // Preconditioned vector
  t_compute_1d d_d;         // Direction vector
  t_compute_1d d_Hdia_inv;  // Diagonal preconditioner
  t_compute_1d d_b;         // RHS vector

  DAT::tdual_float_1d k_s;
  t_compute_1d d_s;
  DAT::t_float_1d h_s;

  DAT::tdual_float_1d k_chi_field;
  typename AT::t_float_1d d_chi_field;

  // CSR matrix
  CSRMatrixType csr_matrix;
  
  // Build matrix and methods for QEq
  void build_csr_matrix();
  void sparse_matvec(t_compute_1d &in, t_compute_1d &out);
  void one_cg_iter();
  void calculate_q();
  void get_chi_field();
  void init_shielding_k();
  void allocate_array();

  KOKKOS_INLINE_FUNCTION
  double calculate_H_k(const double &, const double &) const;

  KOKKOS_INLINE_FUNCTION
  void num_neigh_item(int, bigint &) const;

  int neighflag, nthreads;
};

template <class DeviceType>
struct FixQEqReaxFFKokkos<DeviceType>::KeyValueType {
  using key_type = int;
  using value_type = EV_FLOAT;
  key_type key;
  value_type value;
  KOKKOS_INLINE_FUNCTION KeyValueType() : key(0), value(0.0){};
  KOKKOS_INLINE_FUNCTION KeyValueType(key_type key_in, value_type value_in) :
    key(key_in), value(value_in){};
};

}    // namespace LAMMPS_NS

#endif
#endif
