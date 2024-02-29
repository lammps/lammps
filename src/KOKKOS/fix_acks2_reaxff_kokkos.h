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
FixStyle(acks2/reax/kk,FixACKS2ReaxFFKokkos<LMPDeviceType>);
FixStyle(acks2/reax/kk/device,FixACKS2ReaxFFKokkos<LMPDeviceType>);
FixStyle(acks2/reax/kk/host,FixACKS2ReaxFFKokkos<LMPHostType>);
// clang-format on
#else
// clang-format off
#ifndef LMP_FIX_ACKS2_REAXFF_KOKKOS_H
#define LMP_FIX_ACKS2_REAXFF_KOKKOS_H

#include "fix_acks2_reaxff.h"
#include "kokkos_type.h"
#include "kokkos_base.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagACKS2Zero{};
struct TagACKS2InitMatvec{};
struct TagACKS2SparseMatvec1{};
struct TagACKS2SparseMatvec2{};

template<int NEIGHFLAG>
struct TagACKS2SparseMatvec3_Half{};

struct TagACKS2SparseMatvec3_Full{};
struct TagACKS2Norm1{};
struct TagACKS2Norm2{};
struct TagACKS2Norm3{};
struct TagACKS2Dot1{};
struct TagACKS2Dot2{};
struct TagACKS2Dot3{};
struct TagACKS2Dot4{};
struct TagACKS2Dot5{};
struct TagACKS2Precon1A{};
struct TagACKS2Precon1B{};
struct TagACKS2Precon2{};
struct TagACKS2Add{};
struct TagACKS2ZeroQGhosts{};
struct TagACKS2CalculateQ{};

template<class DeviceType>
class FixACKS2ReaxFFKokkos : public FixACKS2ReaxFF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef double value_type;
  typedef ArrayTypes<DeviceType> AT;
  FixACKS2ReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixACKS2ReaxFFKokkos();

  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void cleanup_copy();

  DAT::tdual_ffloat_1d get_s() {return k_s;}

  KOKKOS_INLINE_FUNCTION
  void num_neigh_item(int, int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Zero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2InitMatvec, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_item(int, int &, const bool &) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_team(const typename Kokkos::TeamPolicy <DeviceType> ::member_type &team, int, int) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_x_item(int, int &, const bool &) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_x_team(const typename Kokkos::TeamPolicy <DeviceType> ::member_type &team, int, int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2SparseMatvec1, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2SparseMatvec2, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2SparseMatvec3_Half<NEIGHFLAG>, const int&) const;

  typedef typename Kokkos::TeamPolicy<DeviceType, TagACKS2SparseMatvec3_Full>::member_type membertype;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagACKS2SparseMatvec3_Full, const membertype &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Norm1, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Norm2, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Norm3, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Dot1, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Dot2, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Dot3, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Dot4, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Dot5, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Precon1A, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Precon1B, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Precon2, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2Add, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2ZeroQGhosts, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagACKS2CalculateQ, const int&) const;

  KOKKOS_INLINE_FUNCTION
  double calculate_H_k(const F_FLOAT &r, const F_FLOAT &shld) const;

  KOKKOS_INLINE_FUNCTION
  double calculate_X_k(const F_FLOAT &r, const F_FLOAT &bcut) const;

  struct params_acks2{
    KOKKOS_INLINE_FUNCTION
    params_acks2(){chi=0;eta=0;gamma=0;bcut_acks2=0;};
    KOKKOS_INLINE_FUNCTION
    params_acks2(int){chi=0;eta=0;gamma=0;bcut_acks2=0;};
    F_FLOAT chi, eta, gamma, bcut_acks2;
  };

 private:
  int inum;
  int allocated_flag, last_allocate;
  int need_dup,prev_last_rows_rank;
  double* buf;

  typename AT::t_int_scalar d_mfill_offset;

  typedef Kokkos::DualView<int***,DeviceType> tdual_int_1d;
  Kokkos::DualView<params_acks2*,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_acks2*, Kokkos::LayoutRight,DeviceType>::t_dev_const params;

  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_f_array_const f;
  typename AT::t_ffloat_1d_randomread mass;
  typename AT::t_ffloat_1d q;
  typename AT::t_int_1d type, mask;
  typename AT::t_tagint_1d tag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist, d_numneigh;

  DAT::tdual_ffloat_1d k_tap;
  typename AT::t_ffloat_1d d_tap;

  DAT::tdual_ffloat_2d k_bcut;
  typename AT::t_ffloat_2d d_bcut;

  typename AT::t_int_1d d_firstnbr;
  typename AT::t_int_1d d_numnbrs;
  typename AT::t_int_1d d_jlist;
  typename AT::t_ffloat_1d d_val;

  typename AT::t_int_1d d_firstnbr_X;
  typename AT::t_int_1d d_numnbrs_X;
  typename AT::t_int_1d d_jlist_X;
  typename AT::t_ffloat_1d d_val_X;

  DAT::tdual_ffloat_1d k_s, k_X_diag, k_chi_field;
  typename AT::t_ffloat_1d d_Hdia_inv,d_Xdia_inv, d_X_diag, d_chi_field, d_b_s,  d_s;
  typename AT::t_ffloat_1d_randomread r_b_s, r_s;

  DAT::tdual_ffloat_1d k_d, k_q_hat, k_z, k_y;
  typename AT::t_ffloat_1d d_p, d_q, d_r, d_d, d_g, d_q_hat, d_r_hat, d_y, d_z, d_bb, d_xx;
  typename AT::t_ffloat_1d_randomread r_p, r_r, r_d;

  DAT::tdual_ffloat_2d k_shield, k_s_hist, k_s_hist_X, k_s_hist_last;
  typename AT::t_ffloat_2d d_shield, d_s_hist, d_s_hist_X, d_s_hist_last;
  typename AT::t_ffloat_2d_randomread r_s_hist, r_s_hist_X, r_s_hist_last;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout> dup_X_diag;
  NonDupScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout> ndup_X_diag;

  DupScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout> dup_bb;
  NonDupScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout> ndup_bb;

  void init_shielding_k();
  void init_hist();
  void allocate_matrix() override;
  void allocate_array();
  void deallocate_array();
  int bicgstab_solve();
  void calculate_Q() override;

  int neighflag;
  int nlocal,nall,nmax,newton_pair;
  int count, isuccess;
  double alpha, beta, omega, cutsq;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um v_buf;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void get_chi_field() override;
  double memory_usage() override;

  void sparse_matvec_acks2(typename AT::t_ffloat_1d &, typename AT::t_ffloat_1d &);
};

template <class DeviceType>
struct FixACKS2ReaxFFKokkosNumNeighFunctor  {
  typedef DeviceType device_type;
  typedef int value_type;
  FixACKS2ReaxFFKokkos<DeviceType> c;
  FixACKS2ReaxFFKokkosNumNeighFunctor(FixACKS2ReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &maxneigh) const {
    c.num_neigh_item(ii, maxneigh);
  }
};

template <class DeviceType, int NEIGHFLAG>
struct FixACKS2ReaxFFKokkosComputeHFunctor {
  int atoms_per_team, vector_length;
  typedef int value_type;
  typedef Kokkos::ScratchMemorySpace<DeviceType> scratch_space;
  FixACKS2ReaxFFKokkos<DeviceType> c;

  FixACKS2ReaxFFKokkosComputeHFunctor(FixACKS2ReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };

  FixACKS2ReaxFFKokkosComputeHFunctor(FixACKS2ReaxFFKokkos<DeviceType> *c_ptr,
                                  int _atoms_per_team, int _vector_length)
      : atoms_per_team(_atoms_per_team), vector_length(_vector_length), c(*c_ptr) {
    c.cleanup_copy();
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &m_fill, const bool &final) const {
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
        Kokkos::View<F_FLOAT **, scratch_space,
                     Kokkos::MemoryUnmanaged>::shmem_size(atoms_per_team,
                                                          vector_length); // s_r
    return shmem_size;
  }
};

template <class DeviceType, int NEIGHFLAG>
struct FixACKS2ReaxFFKokkosComputeXFunctor {
  int atoms_per_team, vector_length;
  typedef int value_type;
  typedef Kokkos::ScratchMemorySpace<DeviceType> scratch_space;
  FixACKS2ReaxFFKokkos<DeviceType> c;

  FixACKS2ReaxFFKokkosComputeXFunctor(FixACKS2ReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };

  FixACKS2ReaxFFKokkosComputeXFunctor(FixACKS2ReaxFFKokkos<DeviceType> *c_ptr,
                                  int _atoms_per_team, int _vector_length)
    : atoms_per_team(_atoms_per_team), vector_length(_vector_length), c(*c_ptr) {
    c.cleanup_copy();
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &m_fill, const bool &final) const {
    c.template compute_x_item<NEIGHFLAG>(ii,m_fill,final);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(
      const typename Kokkos::TeamPolicy<DeviceType>::member_type &team) const {
    c.template compute_x_team<NEIGHFLAG>(team, atoms_per_team, vector_length);
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
        Kokkos::View<F_FLOAT **, scratch_space,
                     Kokkos::MemoryUnmanaged>::shmem_size(atoms_per_team,
                                                          vector_length); // s_r
    return shmem_size;
  }
};

}

#endif
#endif
