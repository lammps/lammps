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

#ifndef LMP_ATOM_VEC_KOKKOS_H
#define LMP_ATOM_VEC_KOKKOS_H

#include "atom_vec.h"
#include "kokkos_type.h"
#include <type_traits>

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

class AtomVecKokkos : public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  virtual ~AtomVecKokkos() {}
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);

  virtual void sync(ExecutionSpace space, unsigned int mask) = 0;
  virtual void modified(ExecutionSpace space, unsigned int mask) = 0;
  virtual void sync_overlapping_device(ExecutionSpace space, unsigned int mask) = 0;

  virtual int
    pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                   const int & iswap, const int nfirst,
                   const int &pbc_flag, const int pbc[]);

  virtual int
    pack_comm_self_fused(const int &n, const DAT::tdual_int_2d &list,
                         const DAT::tdual_int_1d &sendnum_scan,
                         const DAT::tdual_int_1d &firstrecv,
                         const DAT::tdual_int_1d &pbc_flag,
                         const DAT::tdual_int_2d &pbc,
                         const DAT::tdual_int_1d &g2l);

  virtual int
    pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &list,
                     const int & iswap, const DAT::tdual_xfloat_2d &buf,
                     const int &pbc_flag, const int pbc[]);

  virtual void
    unpack_comm_kokkos(const int &n, const int &nfirst,
                       const DAT::tdual_xfloat_2d &buf);

  virtual int
    pack_comm_vel_kokkos(const int &n, const DAT::tdual_int_2d &list,
                         const int & iswap, const DAT::tdual_xfloat_2d &buf,
                         const int &pbc_flag, const int pbc[]);

  virtual void
    unpack_comm_vel_kokkos(const int &n, const int &nfirst,
                           const DAT::tdual_xfloat_2d &buf);

  virtual int
    unpack_reverse_self(const int &n, const DAT::tdual_int_2d &list,
                      const int & iswap, const int nfirst);

  virtual int
    pack_reverse_kokkos(const int &n, const int &nfirst,
                        const DAT::tdual_ffloat_2d &buf);

  virtual void
    unpack_reverse_kokkos(const int &n, const DAT::tdual_int_2d &list,
                          const int & iswap, const DAT::tdual_ffloat_2d &buf);

  virtual int
    pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                       DAT::tdual_xfloat_2d buf,int iswap,
                       int pbc_flag, int *pbc, ExecutionSpace space) = 0;

  virtual void
    unpack_border_kokkos(const int &n, const int &nfirst,
                         const DAT::tdual_xfloat_2d &buf,
                         ExecutionSpace space) = 0;

  virtual int
    pack_border_vel_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                           DAT::tdual_xfloat_2d buf,int iswap,
                           int pbc_flag, int *pbc, ExecutionSpace space) { return 0; }

  virtual void
    unpack_border_vel_kokkos(const int &n, const int &nfirst,
                             const DAT::tdual_xfloat_2d &buf,
                             ExecutionSpace space) {}

  virtual int
    pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &buf,
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space, int dim, X_FLOAT lo, X_FLOAT hi) = 0;

  virtual int
    unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                           int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                           ExecutionSpace space) = 0;


  int no_comm_vel_flag,no_border_vel_flag;

 protected:

  HAT::t_x_array h_x;
  HAT::t_v_array h_v;
  HAT::t_f_array h_f;

  class CommKokkos *commKK;
  size_t buffer_size;
  void* buffer;

  #ifdef KOKKOS_ENABLE_CUDA
  template<class ViewType>
  Kokkos::View<typename ViewType::data_type,
               typename ViewType::array_layout,
               Kokkos::CudaHostPinnedSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> >
  create_async_copy(const ViewType& src) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.span());
       buffer_size = src.span();
    } else if (buffer_size < src.span()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.span());
       buffer_size = src.span();
    }
    return mirror_type(buffer, src.d_view.layout());
  }

  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.span()*sizeof(typename ViewType::value_type));
       buffer_size = src.span();
    } else if (buffer_size < src.span()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.span()*sizeof(typename ViewType::value_type));
       buffer_size = src.span();
    }
    mirror_type tmp_view((typename ViewType::value_type*)buffer, src.d_view.layout());

    if(space == Device) {
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
    if(space == Device)
      src.template sync<LMPDeviceType>();
    else
      src.template sync<LMPHostType>();
  }
  #endif
};

}

#endif

/* ERROR/WARNING messages:

*/
