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
FixStyle(acks2/reaxff/kk,FixACKS2ReaxFFKokkos<LMPDeviceType>);
FixStyle(acks2/reaxff/kk/device,FixACKS2ReaxFFKokkos<LMPDeviceType>);
FixStyle(acks2/reaxff/kk/host,FixACKS2ReaxFFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_ACKS2_REAXFF_KOKKOS_H
#define LMP_FIX_ACKS2_REAXFF_KOKKOS_H

#include "fix_acks2_reaxff.h"

#include "kokkos_type.h"
#include "kokkos_base.h"
#include "neigh_list_kokkos.h"
#include "KokkosSparse_BsrMatrix.hpp"
#include "pair_reaxff_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixACKS2ReaxFFKokkos : public FixACKS2ReaxFF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  using BSRMatrixType = KokkosSparse::Experimental::BsrMatrix<
      double, int, DeviceType, void, int>;

  FixACKS2ReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixACKS2ReaxFFKokkos() override;

  void post_constructor() override;
  void init() override;
  void init_storage() override {}
  void setup_pre_force(int) override; 
  void pre_force(int) override;

  void init_shielding_k();
  void init_bsr();
  void allocate_bsr_matrix();
  void build_bsr_matrix();
  void schur_cg_solve();
  void apply_block_diagonal_preconditioner(
      const typename AT::t_float_1d& in_vec,
      typename AT::t_float_1d& out_vec);
  void compute_B_times_vector(
      const typename AT::t_float_1d& in_vec,  
      typename AT::t_float_1d& out_vec);
  void compute_C_times_vector(
      const typename AT::t_float_1d& in_vec,
      typename AT::t_float_1d& out_vec);
  void sparse_matvec_bsr(
      const typename AT::t_float_1d& x_in,
      typename AT::t_float_1d& b_out);

  void compute_X_diagonal(); 
  void init_matvec();
  void calculate_Q();
  void allocate_array();
  void deallocate_array();
  void allocate_matrix();
  void deallocate_matrix();
  void get_chi_field() override;
  void compute_H() override { /* BSR matrix handles H internally */ }

  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  struct params_acks2 {
    KOKKOS_INLINE_FUNCTION
    params_acks2() {
      chi = 0.0;
      eta = 0.0;
      gamma = 0.0;
      bcut_acks2 = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    params_acks2(int /*i*/) {
      chi = 0.0;
      eta = 0.0;
      gamma = 0.0;
      bcut_acks2 = 0.0;
    }
    F_FLOAT chi, eta, gamma, bcut_acks2;
  };

  KOKKOS_INLINE_FUNCTION
  double calculate_H_k(const F_FLOAT&, const F_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  double calculate_X_k(const F_FLOAT&, const F_FLOAT&) const;

  DAT::tdual_ffloat_1d get_s() {return k_s;}

 private:
  BSRMatrixType bsr_matrix;
  typename AT::t_int_1d d_row_block_counts;  // Pre-allocated for reuse
  
  typedef Kokkos::DualView<int**, DeviceType> tdual_int_2d;
  using KKDeviceType = typename KKDevice<DeviceType>::value;
  
  typedef typename ArrayTypes<DeviceType>::t_float_1d float_1d;
  using d_float_2d = typename ArrayTypes<DeviceType>::t_float_2d;
  using d_float_1d = typename ArrayTypes<DeviceType>::t_float_1d;

  DAT::tdual_float_1d k_s;
  DAT::tdual_float_1d k_X_diag;
  DAT::tdual_float_1d k_chi_field;
  DAT::tdual_float_1d k_b_s;
  typename AT::t_float_1d d_s;
  typename AT::t_float_1d d_X_diag;
  typename AT::t_float_1d d_chi_field;
  typename AT::t_float_1d d_b_s;

  Kokkos::DualView<params_acks2*, Kokkos::LayoutRight, DeviceType> k_params;
  typename Kokkos::DualView<params_acks2*, Kokkos::LayoutRight, DeviceType>::t_dev_const params;

  typename ArrayTypes<DeviceType>::t_neighbors_2d d_neighbors;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread d_ilist;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread d_numneigh;

  typename ArrayTypes<DeviceType>::t_x_array_randomread d_x;
  typename ArrayTypes<DeviceType>::t_int_1d_const d_mask;
  typename ArrayTypes<DeviceType>::t_tagint_1d d_tag;
  typename ArrayTypes<DeviceType>::t_int_1d_const d_type;
  typename ArrayTypes<DeviceType>::t_ffloat_1d d_q;

  // FIXME double check everything again !
  //typename AT::t_x_array_randomread x;
  //typename AT::t_v_array v;
  //typename AT::t_f_array_const f;
  //typename AT::t_ffloat_1d_randomread mass;
  //typename AT::t_ffloat_1d q;
  //typename AT::t_int_1d type, mask;
  //typename AT::t_tagint_1d tag;

  DAT::tdual_int_scalar k_mfill_offset;
  typename AT::t_bigint_scalar d_mfill_offset;

  int NN, nn;
  int allocated_flag, last_allocate;


  DAT::tdual_float_2d k_shield;
  typename AT::t_float_2d d_shield;

  DAT::tdual_float_2d k_bcut;
  typename AT::t_float_2d d_bcut;
  
  DAT::tdual_float_1d k_tap;
  typename AT::t_float_1d d_tap;
  
  int neighflag, nlocal, nall, newton_pair;
  double cutsq;

};

}    // namespace LAMMPS_NS

#endif
#endif
