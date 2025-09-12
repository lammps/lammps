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
FixStyle(qeq/reax/kk,FixQEqReaxFFKokkos<LMPDeviceType>);
FixStyle(qeq/reax/kk/device,FixQEqReaxFFKokkos<LMPDeviceType>);
FixStyle(qeq/reax/kk/host,FixQEqReaxFFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_QEQ_REAXFF_KOKKOS_H
#define LMP_FIX_QEQ_REAXFF_KOKKOS_H

#include "fix_qeq_reaxff.h"
#include "kokkos_type.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

struct TagQEqZero{};
struct TagQEqInitMatvec{};
struct TagQEqSparseMatvec1{};
struct TagQEqZeroQGhosts{};

template<int NEIGHFLAG>
struct TagQEqSparseMatvec2_Half{};

struct TagQEqSparseMatvec2_Full{};
struct TagQEqNorm1{};
struct TagQEqDot1{};
struct TagQEqDot2{};
struct TagQEqDot3{};
struct TagQEqSum1{};
struct TagQEqSum2{};
struct TagQEqCalculateQ{};
struct TagQEqPackForwardComm{};
struct TagQEqUnpackForwardComm{};
struct TagQEqPackExchange{};
struct TagQEqUnpackExchange{};

template<class DeviceType>
class FixQEqReaxFFKokkos : public FixQEqReaxFF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef KK_double2 value_type;
  FixQEqReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixQEqReaxFFKokkos() override;

  void cleanup_copy();
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  KOKKOS_INLINE_FUNCTION
  void num_neigh_item(int, bigint&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqInitMatvec, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_item(int, bigint &, const bool &) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_team(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team, int, int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqSparseMatvec1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqZeroQGhosts, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqSparseMatvec2_Half<NEIGHFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagQEqSparseMatvec2_Half<NEIGHFLAG>>::member_type &team) const;

  typedef typename Kokkos::TeamPolicy<DeviceType, TagQEqSparseMatvec2_Full>::member_type membertype_vec;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqSparseMatvec2_Full, const membertype_vec &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqNorm1, const int&, KK_double2&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqDot1, const int&, KK_double2&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqDot2, const int&, KK_double2&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqDot3, const int&, KK_double2&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqSum1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqSum2, const int&, KK_double2&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqCalculateQ, const int&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT calculate_H_k(const KK_FLOAT &r, const KK_FLOAT &shld) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqUnpackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqPackExchange, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagQEqUnpackExchange, const int&) const;

  int pack_exchange_kokkos(const int &nsend,DAT::tdual_double_2d_lr &buf,
                           DAT::tdual_int_1d k_sendlist,
                           DAT::tdual_int_1d k_copylist,
                           ExecutionSpace space) override;

  void unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nrecv1,int nextrarecv1,
                              ExecutionSpace space) override;

  struct params_qeq{
    KOKKOS_INLINE_FUNCTION
    params_qeq() {chi=0;eta=0;gamma=0;};
    KOKKOS_INLINE_FUNCTION
    params_qeq(int /*i*/) {chi=0;eta=0;gamma=0;};
    KK_FLOAT chi, eta, gamma;
  };

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void sparse_matvec_kokkos(typename AT::t_kkfloat_1d_2 &);

  // There should be a better way to do this for other backends
#if defined(KOKKOS_ENABLE_CUDA)

  // warp length
  static constexpr int vectorsize = 32;

  // team size for sparse mat-vec operations
  static constexpr int spmv_teamsize = 8;

  // custom values for FixQEqReaxFFKokkosComputeHFunctor
  static constexpr int compute_h_vectorsize = vectorsize;
  static constexpr int compute_h_teamsize = 4;
#elif defined(KOKKOS_ENABLE_HIP)

  // wavefront length
  static constexpr int vectorsize = 64;

  // team size for sparse mat-vec operations
  static constexpr int spmv_teamsize = 16;

  // custom values for FixQEqReaxFFKokkosComputeHFunctor
  static constexpr int compute_h_vectorsize = 8; // not a typo, intentionally sub-wavefront
  static constexpr int compute_h_teamsize = 64;
#else
  // dummy values, to be updated
  static constexpr int spmv_teamsize = 4;
  static constexpr int vectorsize = 32;

  static constexpr int compute_h_vectorsize = 1;
  static constexpr int compute_h_teamsize = 32;
#endif

 private:
  int inum,ignum;
  int allocated_flag, last_allocate;
  int need_dup;
  int converged;
  bigint m_cap_big;

  typename AT::t_bigint_scalar d_mfill_offset;

  Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_qeq*, Kokkos::LayoutRight,DeviceType>::t_dev_const params;

  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_kkacc_1d_3_const f;
  //typename AT::t_kkfloat_1d_randomread mass, q;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_kkfloat_1d q;
  typename AT::t_int_1d type, mask;
  typename AT::t_tagint_1d tag;

  DAT::tdual_kkfloat_1d k_q;
  typename AT::t_kkfloat_1d d_q;
  HAT::t_double_1d h_q;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist, d_numneigh;

  DAT::tdual_kkfloat_1d k_tap;
  typename AT::t_kkfloat_1d d_tap;

  typename AT::t_bigint_1d d_firstnbr;
  typename AT::t_int_1d d_numnbrs;
  typename AT::t_int_1d d_jlist;
  typename AT::t_kkfloat_1d d_val;

  DAT::ttransform_kkfloat_1d k_chi_field;
  typename AT::t_kkfloat_1d d_Hdia_inv, d_chi_field;

  DAT::tdual_kkfloat_1d_2 k_o, k_d, k_st;
  typename AT::t_kkfloat_1d_2 d_p, d_o, d_r, d_d, d_b_st, d_st, d_xx;
  HAT::t_kkfloat_1d_2 h_o, h_d, h_st;

  DAT::ttransform_kkfloat_2d k_s_hist, k_t_hist;
  typename AT::t_kkfloat_2d d_s_hist, d_t_hist;
  HAT::t_double_2d_lr h_s_hist, h_t_hist;
  typename AT::t_kkfloat_2d_randomread r_s_hist, r_t_hist;

  DAT::tdual_kkfloat_2d k_shield;
  typename AT::t_kkfloat_2d d_shield;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_FLOAT**, DAT::t_kkfloat_1d_2::array_layout> dup_o;
  NonDupScatterView<KK_FLOAT**, DAT::t_kkfloat_1d_2::array_layout> ndup_o;

  int nsend;
  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d d_buf;
  typename AT::t_int_1d d_copylist;
  typename AT::t_int_1d d_indices;
  typename AT::t_int_1d d_exchange_sendlist;

  void init_shielding_k();
  void init_hist();
  void allocate_matrix() override;
  void allocate_array();

  int cg_solve();
  void calculate_q();

  int neighflag, pack_flag;
  int nlocal,nlocal_last_allocate,nall,nmax,newton_pair;
  int count, isuccess;
  double alpha[2];
  double beta[2];

  double delta, cutsq;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void get_chi_field() override;
};

template <class DeviceType>
struct FixQEqReaxFFKokkosNumNeighFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef bigint value_type;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosNumNeighFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, bigint &totneigh) const {
    c.num_neigh_item(ii, totneigh);
  }
};

template <class DeviceType, int NEIGHFLAG>
struct FixQEqReaxFFKokkosComputeHFunctor {
  int atoms_per_team, vector_length;
  typedef bigint value_type;
  typedef Kokkos::ScratchMemorySpace<DeviceType> scratch_space;
  FixQEqReaxFFKokkos<DeviceType> c;

  FixQEqReaxFFKokkosComputeHFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };

  FixQEqReaxFFKokkosComputeHFunctor(FixQEqReaxFFKokkos<DeviceType> *c_ptr,
                                  int _atoms_per_team, int _vector_length)
    : atoms_per_team(_atoms_per_team), vector_length(_vector_length), c(*c_ptr)  {
    c.cleanup_copy();
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, bigint &m_fill, const bool &final) const {
    c.template compute_h_item<NEIGHFLAG>(ii,m_fill,final);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(
      const typename Kokkos::TeamPolicy<DeviceType>::member_type &team) const {
    c.template compute_h_team<NEIGHFLAG>(team, atoms_per_team, vector_length);
  }

  size_t team_shmem_size(int /*team_size*/) const {
    size_t shmem_size =
        Kokkos::View<int *, scratch_space, Kokkos::MemoryUnmanaged>::shmem_size(
            atoms_per_team) + // s_ilist
        Kokkos::View<int *, scratch_space, Kokkos::MemoryUnmanaged>::shmem_size(
            atoms_per_team) + // s_numnbrs
        Kokkos::View<int *, scratch_space, Kokkos::MemoryUnmanaged>::shmem_size(
            atoms_per_team) + // s_firstnbr
        Kokkos::View<int **, scratch_space, Kokkos::MemoryUnmanaged>::
            shmem_size(atoms_per_team, vector_length) + // s_jtype
        Kokkos::View<int **, scratch_space, Kokkos::MemoryUnmanaged>::
            shmem_size(atoms_per_team, vector_length) + // s_j
        Kokkos::View<KK_FLOAT **, scratch_space,
                     Kokkos::MemoryUnmanaged>::shmem_size(atoms_per_team,
                                                          vector_length); // s_r
    return shmem_size;
  }
};

}

namespace Kokkos {
  // reduction identity must be defined in Kokkos namespace
  template<>
  struct reduction_identity<KK_double2> {
    KOKKOS_FORCEINLINE_FUNCTION static KK_double2 sum() {
      return KK_double2();
    }
  };
}

#endif
#endif
