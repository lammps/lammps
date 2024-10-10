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
// #include <type_traits>

#include <Kokkos_Sort.hpp>

namespace LAMMPS_NS {

union d_ubuf {
  double d;
  int64_t i;
  KOKKOS_INLINE_FUNCTION
  d_ubuf(double arg) : d(arg) {}
  KOKKOS_INLINE_FUNCTION
  d_ubuf(int64_t arg) : i(arg) {}
  KOKKOS_INLINE_FUNCTION
  d_ubuf(int arg) : i(arg) {}
};

class AtomVecKokkos : virtual public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  ~AtomVecKokkos() override;

  using KeyViewType = DAT::t_x_array;
  using BinOp = BinOp3DLAMMPS<KeyViewType>;
  virtual void
    sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) = 0;

  virtual void sync(ExecutionSpace space, unsigned int mask) = 0;
  virtual void modified(ExecutionSpace space, unsigned int mask) = 0;
  virtual void sync_overlapping_device(ExecutionSpace space, unsigned int mask) = 0;

  virtual int
    pack_comm_self(const int &n, const DAT::tdual_int_1d &list,
                   const int nfirst,
                   const int &pbc_flag, const int pbc[]);

  virtual int
    pack_comm_self_fused(const int &n, const DAT::tdual_int_2d &list,
                         const DAT::tdual_int_1d &sendnum_scan,
                         const DAT::tdual_int_1d &firstrecv,
                         const DAT::tdual_int_1d &pbc_flag,
                         const DAT::tdual_int_2d &pbc,
                         const DAT::tdual_int_1d &g2l);

  virtual int
    pack_comm_kokkos(const int &n, const DAT::tdual_int_1d &list,
                     const DAT::tdual_xfloat_2d &buf,
                     const int &pbc_flag, const int pbc[]);


  virtual int
    pack_comm_direct(const int &n, const DAT::tdual_int_2d &list,
                     const DAT::tdual_int_1d &sendnum_scan,
                     const DAT::tdual_int_1d &firstrecv,
                     const DAT::tdual_int_1d &pbc_flag,
                     const DAT::tdual_int_2d &pbc,
                     const DAT::tdual_int_1d &swap2llist,
                     const DAT::tdual_xfloat_1d &buf,
                     const DAT::tdual_int_1d &k_self_flag);

  virtual void
    unpack_comm_kokkos(const int &n, const int &nfirst,
                       const DAT::tdual_xfloat_2d &buf);

  virtual int
    pack_comm_vel_kokkos(const int &n, const DAT::tdual_int_1d &list,
                         const DAT::tdual_xfloat_2d &buf,
                         const int &pbc_flag, const int pbc[]);

  virtual void
    unpack_comm_vel_kokkos(const int &n, const int &nfirst,
                           const DAT::tdual_xfloat_2d &buf);

  virtual int
    pack_reverse_self(const int &n, const DAT::tdual_int_1d &list,
                      const int nfirst);

  virtual int
    pack_reverse_kokkos(const int &n, const int &nfirst,
                        const DAT::tdual_ffloat_2d &buf);

  virtual void
    unpack_reverse_kokkos(const int &n, const DAT::tdual_int_1d &list,
                          const DAT::tdual_ffloat_2d &buf);

  virtual int
    pack_border_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                       DAT::tdual_xfloat_2d buf,
                       int pbc_flag, int *pbc, ExecutionSpace space) = 0;

  virtual void
    unpack_border_kokkos(const int &n, const int &nfirst,
                         const DAT::tdual_xfloat_2d &buf,
                         ExecutionSpace space) = 0;

  virtual int
    pack_border_vel_kokkos(int /*n*/, DAT::tdual_int_1d /*k_sendlist*/,
                           DAT::tdual_xfloat_2d /*buf*/,
                           int /*pbc_flag*/, int * /*pbc*/, ExecutionSpace /*space*/) { return 0; }

  virtual void
    unpack_border_vel_kokkos(const int &/*n*/, const int & /*nfirst*/,
                             const DAT::tdual_xfloat_2d & /*buf*/,
                             ExecutionSpace /*space*/) {}

  virtual int
    pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &buf,
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space) = 0;

  virtual int
    unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                           int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                           ExecutionSpace space,
                           DAT::tdual_int_1d &k_indices) = 0;

  int no_comm_vel_flag,no_border_vel_flag;
  int unpack_exchange_indices_flag;
  int size_exchange;

 protected:
  HAT::t_x_array h_x;
  HAT::t_v_array h_v;
  HAT::t_f_array h_f;

  size_t buffer_size;
  void* buffer;

  DAT::tdual_int_1d k_count;

 public:

  #ifdef LMP_KOKKOS_GPU
  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
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
    mirror_type tmp_view((typename ViewType::value_type*)buffer, src.d_view.layout());

    if (space == Device) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.h_view),
      Kokkos::deep_copy(LMPHostType(),src.d_view,tmp_view);
      src.clear_sync_state();
    } else {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.d_view),
      Kokkos::deep_copy(LMPHostType(),src.h_view,tmp_view);
      src.clear_sync_state();
    }
  }
  #else
  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
    if (space == Device)
      src.template sync<LMPDeviceType>();
    else
      src.template sync<LMPHostType>();
  }
  #endif
};

}

#endif
