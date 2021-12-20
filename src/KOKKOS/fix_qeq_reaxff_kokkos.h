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

#ifdef HIP_OPT_CG_SOLVE_FUSED
struct TagSparseMatvec13 {};
struct TagSparseMatvec13Vector {};
struct TagSparseMatvec2Fused {};
struct TagSparseMatvec2FusedVector {};
#else
struct TagSparseMatvec1 {};
struct TagSparseMatvec1Vector {};
struct TagSparseMatvec2 {};
struct TagSparseMatvec2Vector {};
struct TagSparseMatvec3 {};
struct TagSparseMatvec3Vector {};
struct TagZeroQGhosts{};
#endif

struct TagFixQEqReaxFFPackForwardComm {};
struct TagFixQEqReaxFFUnpackForwardComm {};

template<class DeviceType>
class FixQEqReaxFFKokkos : public FixQEqReaxFF, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  FixQEqReaxFFKokkos(class LAMMPS *, int, char **);
  ~FixQEqReaxFFKokkos();

  void cleanup_copy();
  void init();
  void setup_pre_force(int);
  void pre_force(int);

  KOKKOS_INLINE_FUNCTION
  void num_neigh_item(int, int&) const;

  KOKKOS_INLINE_FUNCTION
  void zero_item(int) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_item(int, int &, const bool &) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void compute_h_team(const typename Kokkos::TeamPolicy <DeviceType> ::member_type &team, int, int) const;

  KOKKOS_INLINE_FUNCTION
  void matvec_item(int) const;

#ifdef HIP_OPT_CG_SOLVE_FUSED

  KOKKOS_INLINE_FUNCTION
  void sparse22_fused_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void sparse12_32_item(int) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2Fused> ::member_type membertype2fused;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec2Fused, const membertype2fused &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2FusedVector> ::member_type membertype2fusedvec;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec2FusedVector, const membertype2fusedvec &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec13> ::member_type membertype13;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec13, const membertype13 &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec13Vector> ::member_type membertype13vec;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec13Vector, const membertype13vec &team) const;

  KOKKOS_INLINE_FUNCTION
  void vecsum2_fused_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void norm12_item(int, F_FLOAT2&) const;

  KOKKOS_INLINE_FUNCTION
  void dot11_item(int, F_FLOAT2&) const;

  KOKKOS_INLINE_FUNCTION
  void dot22_item(int, F_FLOAT2&) const;

  KOKKOS_INLINE_FUNCTION
  void precon12_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void precon_fused_item(int, F_FLOAT2&) const;

#else

  KOKKOS_INLINE_FUNCTION
  void sparse12_item(int) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void sparse13_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void sparse22_item(int) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void sparse23_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void sparse32_item(int) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void sparse33_item(int) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec1> ::member_type membertype1;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec1, const membertype1 &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec1Vector> ::member_type membertype1vec;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec1Vector, const membertype1vec &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2> ::member_type membertype2;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec2, const membertype2 &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2Vector> ::member_type membertype2vec;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec2Vector, const membertype2vec &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec3> ::member_type membertype3;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec3, const membertype3 &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec3Vector> ::member_type membertype3vec;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec3Vector, const membertype3vec &team) const;

  KOKKOS_INLINE_FUNCTION
  void vecsum2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double norm1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double norm2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double dot1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double dot2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void precon1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void precon2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double precon_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagZeroQGhosts, const int&) const;

#endif

  

  KOKKOS_INLINE_FUNCTION
  double vecacc1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double vecacc2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void calculate_q_item(int) const;

  KOKKOS_INLINE_FUNCTION
  double calculate_H_k(const F_FLOAT &r, const F_FLOAT &shld) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixQEqReaxFFPackForwardComm, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixQEqReaxFFUnpackForwardComm, const int&) const;

  struct params_qeq{
    KOKKOS_INLINE_FUNCTION
    params_qeq() {chi=0;eta=0;gamma=0;};
    KOKKOS_INLINE_FUNCTION
    params_qeq(int /*i*/) {chi=0;eta=0;gamma=0;};
    F_FLOAT chi, eta, gamma;
  };

  int pack_forward_comm_fix_kokkos(int, DAT::tdual_int_2d, int, DAT::tdual_xfloat_1d&,
                       int, int *);
  void unpack_forward_comm_fix_kokkos(int, int, DAT::tdual_xfloat_1d&);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

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

  typename AT::t_int_scalar d_mfill_offset;

  typedef Kokkos::DualView<int***,DeviceType> tdual_int_1d;
  Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_qeq*, Kokkos::LayoutRight,DeviceType>::t_dev_const params;

  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_f_array_const f;
  //typename ArrayTypes<DeviceType>::t_float_1d_randomread mass, q;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread mass;
  typename ArrayTypes<DeviceType>::t_float_1d q;
  typename ArrayTypes<DeviceType>::t_int_1d type, mask;
  typename ArrayTypes<DeviceType>::t_tagint_1d tag;

  DAT::tdual_float_1d k_q;
  typename AT::t_float_1d d_q;
  HAT::t_float_1d h_q;

  typename ArrayTypes<DeviceType>::t_neighbors_2d d_neighbors;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread d_ilist, d_numneigh;

  DAT::tdual_ffloat_1d k_tap;
  typename AT::t_ffloat_1d d_tap;

  typename AT::t_int_1d d_firstnbr;
  typename AT::t_int_1d d_numnbrs;
  typename AT::t_int_1d d_jlist;
  typename AT::t_ffloat_1d d_val;

  DAT::tdual_ffloat_1d k_t, k_s, k_chi_field;
  typename AT::t_ffloat_1d d_Hdia_inv, d_b_s, d_b_t, d_t, d_s, d_chi_field;
  HAT::t_ffloat_1d h_t, h_s;
  typename AT::t_ffloat_1d_randomread r_b_s, r_b_t, r_t, r_s;

#ifdef HIP_OPT_CG_SOLVE_FUSED
  DAT::tdual_ffloat2_1d k_o_fused, k_d_fused;
  typename AT::t_ffloat2_1d d_p_fused, d_o_fused, d_r_fused, d_d_fused;
  HAT::t_ffloat2_1d h_o_fused, h_d_fused;
  int converged = 0;
#else
  DAT::tdual_ffloat_1d k_o, k_d;
  typename AT::t_ffloat_1d d_p, d_o, d_r, d_d;
  HAT::t_ffloat_1d h_o, h_d;

  Kokkos::Experimental::ScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated> dup_o;
  Kokkos::Experimental::ScatterView<F_FLOAT*, typename AT::t_ffloat_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated> ndup_o;
#endif

  typename AT::t_ffloat_1d_randomread r_p, r_o, r_r, r_d;

  DAT::tdual_ffloat_2d k_shield, k_s_hist, k_t_hist;
  typename AT::t_ffloat_2d d_shield, d_s_hist, d_t_hist;
  HAT::t_ffloat_2d h_s_hist, h_t_hist;
  typename AT::t_ffloat_2d_randomread r_s_hist, r_t_hist;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_xfloat_1d_um d_buf;

  void init_shielding_k();
  void init_hist();
  void allocate_matrix();
  void allocate_array();

#ifdef HIP_OPT_CG_SOLVE_FUSED
  int cg_solve_fused();
#else
  int cg_solve1();
  int cg_solve2();
#endif
  void calculate_q();

  int neighflag, pack_flag;
  int nlocal,nall,nmax,newton_pair;
  int count, isuccess;
#ifdef HIP_OPT_CG_SOLVE_FUSED
  F_FLOAT alpha[2];
  F_FLOAT beta[2];
#else
  double alpha, beta;
#endif

  double delta, cutsq;

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void get_chi_field();
};

template <class DeviceType>
struct FixQEqReaxFFKokkosNumNeighFunctor  {
  typedef DeviceType  device_type ;
  typedef int value_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosNumNeighFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &maxneigh) const {
    c.num_neigh_item(ii, maxneigh);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosMatVecFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosMatVecFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.matvec_item(ii);
  }
};

template <class DeviceType, int NEIGHFLAG>
struct FixQEqReaxFFKokkosComputeHFunctor {
  int atoms_per_team, vector_length;
  typedef int value_type;
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

template <class DeviceType>
struct FixQEqReaxFFKokkosZeroFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosZeroFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.zero_item(ii);
  }
};

#ifdef HIP_OPT_CG_SOLVE_FUSED
// Functors exclusive to fused CG

template <class DeviceType>
struct FixQEqReaxFFKokkosSparse22FusedFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse22FusedFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse22_fused_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosSparse12_32Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse12_32Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse12_32_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosVecSum2FusedFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosVecSum2FusedFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.vecsum2_fused_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosNorm12Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef F_FLOAT2 value_type;
  FixQEqReaxFFKokkosNorm12Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    c.norm12_item(ii, tmp);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosDot11Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef F_FLOAT2 value_type;
  FixQEqReaxFFKokkosDot11Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    c.dot11_item(ii, tmp);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosDot22Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef F_FLOAT2 value_type;
  FixQEqReaxFFKokkosDot22Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    c.dot22_item(ii, tmp);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosPrecon12Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosPrecon12Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.precon12_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosPreconFusedFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef F_FLOAT2 value_type;
  FixQEqReaxFFKokkosPreconFusedFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    c.precon_fused_item(ii, tmp);
  }
};

#else
// Functors exclusive to non-fused CG

template <class DeviceType>
struct FixQEqReaxFFKokkosSparse12Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse12Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse12_item(ii);
  }
};

template <class DeviceType,int NEIGHFLAG>
struct FixQEqReaxFFKokkosSparse13Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse13Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse13_item<NEIGHFLAG>(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosSparse22Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse22Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse22_item(ii);
  }
};

template <class DeviceType,int NEIGHFLAG>
struct FixQEqReaxFFKokkosSparse23Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse23Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse23_item<NEIGHFLAG>(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosSparse32Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse32Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse32_item(ii);
  }
};

template <class DeviceType,int NEIGHFLAG>
struct FixQEqReaxFFKokkosSparse33Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosSparse33Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse33_item<NEIGHFLAG>(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosVecSum2Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosVecSum2Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.vecsum2_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosNorm1Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosNorm1Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.norm1_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosNorm2Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosNorm2Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.norm2_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosDot1Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosDot1Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.dot1_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosDot2Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosDot2Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.dot2_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosPrecon1Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosPrecon1Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.precon1_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosPrecon2Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosPrecon2Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.precon2_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosPreconFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosPreconFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.precon_item(ii);
  }
};

#endif


template <class DeviceType>
struct FixQEqReaxFFKokkosVecAcc1Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosVecAcc1Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.vecacc1_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosVecAcc2Functor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  typedef double value_type;
  FixQEqReaxFFKokkosVecAcc2Functor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.vecacc2_item(ii);
  }
};

template <class DeviceType>
struct FixQEqReaxFFKokkosCalculateQFunctor  {
  typedef DeviceType  device_type ;
  FixQEqReaxFFKokkos<DeviceType> c;
  FixQEqReaxFFKokkosCalculateQFunctor(FixQEqReaxFFKokkos<DeviceType>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.calculate_q_item(ii);
  }
};

}

namespace Kokkos {
  // reduction identity must be defined in Kokkos namespace
  template<>
  struct reduction_identity< F_FLOAT2 > {
    KOKKOS_FORCEINLINE_FUNCTION static F_FLOAT2 sum() {
      return F_FLOAT2();
    }
  };
}

#endif
#endif
