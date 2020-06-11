/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qeq/reax/kk,FixQEqReaxKokkos<Device>)
FixStyle(qeq/reax/kk/device,FixQEqReaxKokkos<Device>)
FixStyle(qeq/reax/kk/host,FixQEqReaxKokkos<Host>)

#else

#ifndef LMP_FIX_QEQ_REAX_KOKKOS_H
#define LMP_FIX_QEQ_REAX_KOKKOS_H

#include "fix_qeq_reax.h"
#include "kokkos_type.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct TagSparseMatvec1 {};
struct TagSparseMatvec2 {};
struct TagSparseMatvec3 {};
struct TagZeroQGhosts{};

template<ExecutionSpace Space>
class FixQEqReaxKokkos : public FixQEqReax {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  FixQEqReaxKokkos(class LAMMPS *, int, char **);
  ~FixQEqReaxKokkos();

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

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2> ::member_type membertype2;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec2, const membertype2 &team) const;

  typedef typename Kokkos::TeamPolicy <DeviceType, TagSparseMatvec3> ::member_type membertype3;
  KOKKOS_INLINE_FUNCTION
  void operator() (TagSparseMatvec3, const membertype3 &team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagZeroQGhosts, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void vecsum2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT norm1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT norm2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT dot1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT dot2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void precon1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void precon2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT precon_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT vecacc1_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT vecacc2_item(int) const;

  KOKKOS_INLINE_FUNCTION
  void calculate_q_item(int) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT calculate_H_k(const KK_FLOAT &r, const KK_FLOAT &shld) const;

  struct params_qeq{
    KOKKOS_INLINE_FUNCTION
    params_qeq(){chi=0;eta=0;gamma=0;};
    KOKKOS_INLINE_FUNCTION
    params_qeq(int i){chi=0;eta=0;gamma=0;};
    KK_FLOAT chi, eta, gamma;
  };

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int inum;
  int allocated_flag;
  int need_dup, cuda_flag;

  typename AT::t_int_scalar d_mfill_offset;

  typedef Kokkos::DualView<int***,DeviceType> tdual_int_1d;
  Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_qeq*, Kokkos::LayoutRight,DeviceType>::t_dev_const params;

  typename AT::t_float_1d_3_lr x;
  typename AT::t_float_1d_3 v;
  typename AT::t_float_1d_3_const f;
  //typename AT::t_float_1d_randomread mass, q;
  typename AT::t_float_1d_randomread mass;
  typename AT::t_float_1d q;
  typename AT::t_int_1d type, mask;
  typename AT::t_tagint_1d tag;

  DAT::tdual_float_1d k_q;
  typename AT::t_float_1d d_q;
  HAT::t_float_1d h_q;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist, d_numneigh;

  DAT::tdual_float_1d k_tap;
  typename AT::t_float_1d d_tap;

  typename AT::t_int_1d d_firstnbr;
  typename AT::t_int_1d d_numnbrs;
  typename AT::t_int_1d d_jlist;
  typename AT::t_float_1d d_val;

  DAT::tdual_float_1d k_t, k_s;
  typename AT::t_float_1d d_Hdia_inv, d_b_s, d_b_t, d_t, d_s;
  HAT::t_float_1d h_t, h_s;
  typename AT::t_float_1d_randomread r_b_s, r_b_t, r_t, r_s;

  DAT::tdual_float_1d k_o, k_d;
  typename AT::t_float_1d d_p, d_o, d_r, d_d;
  HAT::t_float_1d h_o, h_d;
  typename AT::t_float_1d_randomread r_p, r_o, r_r, r_d;

  DAT::tdual_float_2d k_shield, k_s_hist, k_t_hist;
  typename AT::t_float_2d d_shield, d_s_hist, d_t_hist;
  HAT::t_float_2d h_s_hist, h_t_hist;
  typename AT::t_float_2d_randomread r_s_hist, r_t_hist;

  Kokkos::Experimental::ScatterView<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated> dup_o;
  Kokkos::Experimental::ScatterView<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout, typename KKDevice<DeviceType>::value, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated> ndup_o;

  void init_shielding_k();
  void init_hist();
  void allocate_matrix();
  void allocate_array();
  void cg_solve1();
  void cg_solve2();
  void calculate_q();

  int neighflag, pack_flag;
  int nlocal,nall,nmax,newton_pair;
  int count, isuccess;
  KK_FLOAT alpha, beta, delta, cutsq;

  int iswap;
  int first;
  typename AT::t_int_2d d_sendlist;
  typename AT::t_float_1d_um v_buf;

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosNumNeighFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  typedef int value_type ;
  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosNumNeighFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, int &maxneigh) const {
    c.num_neigh_item(ii, maxneigh);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosMatVecFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosMatVecFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.matvec_item(ii);
  }
};

template <ExecutionSpace Space, int NEIGHFLAG>
struct FixQEqReaxKokkosComputeHFunctor {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  int atoms_per_team, vector_length;
  typedef int value_type;
  typedef Kokkos::ScratchMemorySpace<DeviceType> scratch_space;
  FixQEqReaxKokkos<Space> c;

  FixQEqReaxKokkosComputeHFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };

  FixQEqReaxKokkosComputeHFunctor(FixQEqReaxKokkos<Space> *c_ptr,
                                  int _atoms_per_team, int _vector_length)
      : c(*c_ptr), atoms_per_team(_atoms_per_team),
        vector_length(_vector_length) {
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

  size_t team_shmem_size(int team_size) const {
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

template <ExecutionSpace Space>
struct FixQEqReaxKokkosZeroFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosZeroFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.zero_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosSparse12Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse12Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse12_item(ii);
  }
};

template <ExecutionSpace Space,int NEIGHFLAG>
struct FixQEqReaxKokkosSparse13Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse13Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse13_item<NEIGHFLAG>(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosSparse22Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse22Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse22_item(ii);
  }
};

template <ExecutionSpace Space,int NEIGHFLAG>
struct FixQEqReaxKokkosSparse23Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse23Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse23_item<NEIGHFLAG>(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosSparse32Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse32Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.sparse32_item(ii);
  }
};

template <ExecutionSpace Space,int NEIGHFLAG>
struct FixQEqReaxKokkosSparse33Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosSparse33Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.template sparse33_item<NEIGHFLAG>(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosVecSum2Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosVecSum2Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.vecsum2_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosNorm1Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosNorm1Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.norm1_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosNorm2Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosNorm2Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.norm2_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosDot1Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosDot1Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.dot1_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosDot2Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosDot2Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.dot2_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosPrecon1Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosPrecon1Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.precon1_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosPrecon2Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosPrecon2Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.precon2_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosPreconFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosPreconFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.precon_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosVecAcc1Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosVecAcc1Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.vecacc1_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosVecAcc2Functor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  typedef KK_FLOAT value_type;
  FixQEqReaxKokkosVecAcc2Functor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii, value_type &tmp) const {
    tmp += c.vecacc2_item(ii);
  }
};

template <ExecutionSpace Space>
struct FixQEqReaxKokkosCalculateQFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;

  FixQEqReaxKokkos<Space> c;
  FixQEqReaxKokkosCalculateQFunctor(FixQEqReaxKokkos<Space>* c_ptr):c(*c_ptr) {
    c.cleanup_copy();
  };
  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    c.calculate_q_item(ii);
  }
};

}

#endif
#endif
