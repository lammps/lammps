/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXPERIMENTAL_VIEWETIAVAIL_HPP
#define KOKKOS_EXPERIMENTAL_VIEWETIAVAIL_HPP

namespace Kokkos {
namespace Impl {

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          int Rank, typename iType>
struct ViewCopyETIAvail {
  enum { value = false };
};

#define KOKKOS_IMPL_VIEWCOPY_ETI_AVAIL(DATATYPE, LAYOUTA, LAYOUTB, EXECSPACE, \
                                       ITYPE)                                 \
  template <>                                                                 \
  struct ViewCopyETIAvail<                                                    \
      Kokkos::View<DATATYPE, LAYOUTA,                                         \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,         \
                   Kokkos::MemoryTraits<0>>,                                  \
      Kokkos::View<const DATATYPE, LAYOUTB,                                   \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,         \
                   Kokkos::MemoryTraits<0>>,                                  \
      Kokkos::LayoutLeft, EXECSPACE, Kokkos::View<DATATYPE>::rank, ITYPE> {   \
    enum { value = true };                                                    \
  };                                                                          \
  template <>                                                                 \
  struct ViewCopyETIAvail<                                                    \
      Kokkos::View<DATATYPE, LAYOUTA,                                         \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,         \
                   Kokkos::MemoryTraits<0>>,                                  \
      Kokkos::View<const DATATYPE, LAYOUTB,                                   \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,         \
                   Kokkos::MemoryTraits<0>>,                                  \
      Kokkos::LayoutRight, EXECSPACE, Kokkos::View<DATATYPE>::rank, ITYPE> {  \
    enum { value = true };                                                    \
  };

template <class ViewType, class Layout, class ExecSpace, int Rank,
          typename iType>
struct ViewFillETIAvail {
  enum { value = false };
};

#define KOKKOS_IMPL_VIEWFILL_ETI_AVAIL(DATATYPE, LAYOUT, EXECSPACE, ITYPE)   \
  template <>                                                                \
  struct ViewFillETIAvail<                                                   \
      Kokkos::View<DATATYPE, LAYOUT,                                         \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,        \
                   Kokkos::MemoryTraits<0>>,                                 \
      Kokkos::LayoutLeft, EXECSPACE, Kokkos::View<DATATYPE>::rank, ITYPE> {  \
    enum { value = true };                                                   \
  };                                                                         \
  template <>                                                                \
  struct ViewFillETIAvail<                                                   \
      Kokkos::View<DATATYPE, LAYOUT,                                         \
                   Kokkos::Device<EXECSPACE, Kokkos::AnonymousSpace>,        \
                   Kokkos::MemoryTraits<0>>,                                 \
      Kokkos::LayoutRight, EXECSPACE, Kokkos::View<DATATYPE>::rank, ITYPE> { \
    enum { value = true };                                                   \
  };

}  // namespace Impl
}  // namespace Kokkos
#ifdef KOKKOS_ENABLE_ETI
#ifdef KOKKOS_ENABLE_Serial
#include <Serial/Kokkos_Serial_ViewCopyETIAvail.hpp>
#endif
#ifdef KOKKOS_ENABLE_OPENMP
#include <OpenMP/Kokkos_OpenMP_ViewCopyETIAvail.hpp>
#endif
#ifdef KOKKOS_ENABLE_THREADS
#include <Threads/Kokkos_Threads_ViewCopyETIAvail.hpp>
#endif
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_ViewCopyETIAvail.hpp>
#endif
#ifdef KOKKOS_ENABLE_ROCM
#include <ROCm/Kokkos_ROCm_ViewCopyETIAvail.hpp>
#endif
#endif

#endif
