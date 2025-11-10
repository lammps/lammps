// clang-format off
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

#ifndef LMP_ATOM_VEC_KOKKOS_H
#define LMP_ATOM_VEC_KOKKOS_H

#include "atom_vec.h"           //  IWYU pragma: export

#include "kokkos_type.h"
#include <type_traits>

#include <Kokkos_Sort.hpp>

namespace LAMMPS_NS {

class AtomVecKokkos : virtual public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  ~AtomVecKokkos() override;

  using KeyViewType = DAT::t_kkfloat_1d_3_lr;
  using BinOp = BinOp3DLAMMPS<KeyViewType>;
  virtual void
    sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) = 0;

  virtual void sync(ExecutionSpace space, unsigned int mask) = 0;
  virtual void modified(ExecutionSpace space, unsigned int mask) = 0;
  virtual void sync_pinned(ExecutionSpace space, unsigned int mask, int async_flag = 0) = 0;

  virtual int
    pack_comm_self(const int &n, const DAT::tdual_int_1d &list,
                   const int nfirst,
                   const int &pbc_flag, const int pbc[]);

  virtual int
    pack_comm_self_fused(const int &n, const DAT::tdual_int_2d_lr &list,
                         const DAT::tdual_int_1d &sendnum_scan,
                         const DAT::tdual_int_1d &firstrecv,
                         const DAT::tdual_int_1d &pbc_flag,
                         const DAT::tdual_int_2d &pbc,
                         const DAT::tdual_int_1d &g2l);

  virtual int
    pack_comm_kokkos(const int &n, const DAT::tdual_int_1d &list,
                     const DAT::tdual_double_2d_lr &buf,
                     const int &pbc_flag, const int pbc[]);

  virtual void
    unpack_comm_kokkos(const int &n, const int &nfirst,
                       const DAT::tdual_double_2d_lr &buf);

  virtual int
    pack_comm_vel_kokkos(const int &n, const DAT::tdual_int_1d &list,
                         const DAT::tdual_double_2d_lr &buf,
                         const int &pbc_flag, const int pbc[]);

  virtual void
    unpack_comm_vel_kokkos(const int &n, const int &nfirst,
                           const DAT::tdual_double_2d_lr &buf);

  virtual int
    pack_reverse_self(const int &n, const DAT::tdual_int_1d &list,
                      const int nfirst);

  virtual int
    pack_reverse_kokkos(const int &n, const int &nfirst,
                        const DAT::tdual_double_2d_lr &buf);

  virtual void
    unpack_reverse_kokkos(const int &n, const DAT::tdual_int_1d &list,
                          const DAT::tdual_double_2d_lr &buf);

  virtual int
    pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                       DAT::tdual_double_2d_lr buf,
                       int pbc_flag, int *pbc, ExecutionSpace space) = 0;

  virtual void
    unpack_border_kokkos(const int &n, const int &nfirst,
                         const DAT::tdual_double_2d_lr &buf,
                         ExecutionSpace space) = 0;

  virtual int
    pack_border_vel_kokkos(int /*n*/, DAT::tdual_int_1d /*k_sendlist*/,
                           DAT::tdual_double_2d_lr /*buf*/,
                           int /*pbc_flag*/, int * /*pbc*/, ExecutionSpace /*space*/) { return 0; }

  virtual void
    unpack_border_vel_kokkos(const int &/*n*/, const int & /*nfirst*/,
                             const DAT::tdual_double_2d_lr & /*buf*/,
                             ExecutionSpace /*space*/) {}

  virtual int
    pack_exchange_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space) = 0;

  virtual int
    unpack_exchange_kokkos(DAT::tdual_double_2d_lr &k_buf, int nrecv,
                           int nlocal, int dim, double lo, double hi,
                           ExecutionSpace space,
                           DAT::tdual_int_1d &k_indices) = 0;

  int no_comm_vel_flag,no_border_vel_flag;
  int unpack_exchange_indices_flag;
  int size_exchange;

 protected:
  HAT::t_kkfloat_1d_3_lr h_x;
  HAT::t_kkfloat_1d_3 h_v;
  HAT::t_kkacc_1d_3 h_f;

  size_t buffer_size;
  void* buffer;

  DAT::tdual_int_1d k_count;

 public:

  #ifdef LMP_KOKKOS_GPU
  template<class ViewType>
  void perform_pinned_copy(ViewType& src, unsigned int space, int async_flag = 0) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same_v<typename ViewType::execution_space,LMPDeviceType>,
                   LMPPinnedHostType,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer_size = src.span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_malloc<LMPPinnedHostType>(buffer_size);
    } else if (buffer_size < src.span()*sizeof(typename ViewType::value_type)) {
       buffer_size = src.span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_realloc<LMPPinnedHostType>(buffer,buffer_size);
    }

    mirror_type tmp_view((typename ViewType::value_type*)buffer, src.view_device().layout());

    if (src.view_device().data() && space == Device) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.view_host()),
      Kokkos::deep_copy(LMPHostType(),src.view_device(),tmp_view);
      src.clear_sync_state();
      if (!async_flag) Kokkos::fence(); // change to less agressive fence?
    } else if (src.view_host().data()) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.view_device()),
      Kokkos::deep_copy(LMPHostType(),src.view_host(),tmp_view);
      src.clear_sync_state();
      if (!async_flag) Kokkos::fence(); // change to less agressive fence?
    }
  }
  #else
  template<class ViewType>
  void perform_pinned_copy(ViewType& src, unsigned int space, int /*async_flag*/ = 0) {
    if (space == Device)
      src.sync_device();
    else
      src.sync_host();
  }
  #endif

  #ifdef LMP_KOKKOS_GPU
  template<class TransformViewType>
  void perform_pinned_copy_transform(TransformViewType& src, unsigned int space, int async_flag = 0) {
    typedef typename TransformViewType::kk_view ViewType;
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same_v<typename ViewType::execution_space,LMPDeviceType>,
                   LMPPinnedHostType,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer_size = src.view_device().span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_malloc<LMPPinnedHostType>(buffer_size);
    } else if (buffer_size < src.view_device().span()*sizeof(typename ViewType::value_type)) {
       buffer_size = src.view_device().span()*sizeof(typename ViewType::value_type);
       buffer = Kokkos::kokkos_realloc<LMPPinnedHostType>(buffer,buffer_size);
    }

    if (space == Device)
      src.sync_device(buffer,async_flag);
    else if (space == Host)
      src.sync_host(buffer,async_flag);
    else if (space == HostKK)
      src.sync_hostkk(buffer,async_flag);
  }
  #else
  template<class TransformViewType>
  void perform_pinned_copy_transform(TransformViewType& src, unsigned int space, int /*async_flag*/ = 0) {
    if (space == Device)
      src.sync_device();
    else if (space == Host)
      src.sync_host();
    else if (space == HostKK)
      src.sync_hostkk();
  }
  #endif
};

}

#endif
