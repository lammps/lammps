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

#ifndef KOKKOS_FUNCTORADAPTER_HPP
#define KOKKOS_FUNCTORADAPTER_HPP

#include <cstddef>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class Enable = void>
struct ReduceFunctorHasInit {
  enum : bool { value = false };
};

// The else clause idiom failed with NVCC+MSVC, causing some symbols not being
// compiled for the device. The code in there is anyway sketchy, and likely not
// standard compliant (just happens to work on all compilers we ever used)
// We intend to replace all of this long term with proper detection idiom.
#if defined(KOKKOS_COMPILER_MSVC) || defined(KOKKOS_IMPL_WINDOWS_CUDA)
template <class>
using impl_void_t_workaround = void;

template <class F>
using init_archetype = decltype(&F::init);

template <class FunctorType>
struct ReduceFunctorHasInit<
    FunctorType, impl_void_t_workaround<init_archetype<FunctorType>>> {
  enum : bool { value = true };
};
#else
template <class FunctorType>
struct ReduceFunctorHasInit<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::init)>::type> {
  enum : bool { value = true };
};
// FIXME_SYCL not all compilers distinguish between the FunctorType::init and
// the FunctorType::template init<> specialization
#ifdef KOKKOS_ENABLE_SYCL
template <class FunctorType>
struct ReduceFunctorHasInit<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::template init<>)>::type> {
  enum : bool { value = true };
};
#endif
#endif

template <class FunctorType, class Enable = void>
struct ReduceFunctorHasJoin {
  enum : bool { value = false };
};

#if defined(KOKKOS_COMPILER_MSVC) || defined(KOKKOS_IMPL_WINDOWS_CUDA)
template <class F>
using join_archetype = decltype(&F::join);

template <class FunctorType>
struct ReduceFunctorHasJoin<
    FunctorType, impl_void_t_workaround<join_archetype<FunctorType>>> {
  enum : bool { value = true };
};
#else
template <class FunctorType>
struct ReduceFunctorHasJoin<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::join)>::type> {
  enum : bool { value = true };
};
// FIXME_SYCL not all compilers distinguish between the FunctorType::join and
// the FunctorType::template join<> specialization
#ifdef KOKKOS_ENABLE_SYCL
template <class FunctorType>
struct ReduceFunctorHasJoin<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::template join<>)>::type> {
  enum : bool { value = true };
};
#endif
#endif

template <class FunctorType, class Enable = void>
struct ReduceFunctorHasFinal {
  enum : bool { value = false };
};

#if defined(KOKKOS_COMPILER_MSVC) || defined(KOKKOS_IMPL_WINDOWS_CUDA)
template <class F>
using final_archetype = decltype(&F::final);

template <class FunctorType>
struct ReduceFunctorHasFinal<
    FunctorType, impl_void_t_workaround<final_archetype<FunctorType>>> {
  enum : bool { value = true };
};
#else
template <class FunctorType>
struct ReduceFunctorHasFinal<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::final)>::type> {
  enum : bool { value = true };
};
// FIXME_SYCL not all compilers distinguish between the FunctorType::final and
// the FunctorType::template final<> specialization
#ifdef KOKKOS_ENABLE_SYCL
template <class FunctorType>
struct ReduceFunctorHasFinal<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::template final<>)>::type> {
  enum : bool { value = true };
};
#endif
#endif

template <class FunctorType, class Enable = void>
struct ReduceFunctorHasShmemSize {
  enum : bool { value = false };
};

#if defined(KOKKOS_COMPILER_MSVC) || defined(KOKKOS_IMPL_WINDOWS_CUDA)
template <class F>
using shmemsize_archetype = decltype(&F::team_shmem_size);

template <class FunctorType>
struct ReduceFunctorHasShmemSize<
    FunctorType, impl_void_t_workaround<shmemsize_archetype<FunctorType>>> {
  enum : bool { value = true };
};
#else
template <class FunctorType>
struct ReduceFunctorHasShmemSize<
    FunctorType,
    typename std::enable_if<0 < sizeof(&FunctorType::team_shmem_size)>::type> {
  enum : bool { value = true };
};
// FIXME_SYCL not all compilers distinguish between the
// FunctorType::team_shmem_size and the FunctorType::template team_shmem_size<>
// specialization
#ifdef KOKKOS_ENABLE_SYCL
template <class FunctorType>
struct ReduceFunctorHasShmemSize<
    FunctorType,
    typename std::enable_if<
        0 < sizeof(&FunctorType::template team_shmem_size<>)>::type> {
  enum : bool { value = true };
};
#endif
#endif

template <class FunctorType, class ArgTag, class Enable = void>
struct FunctorDeclaresValueType : public std::false_type {};

template <class FunctorType, class ArgTag>
struct FunctorDeclaresValueType<
    FunctorType, ArgTag,
    typename Impl::enable_if_type<typename FunctorType::value_type>::type>
    : public std::true_type {};

template <class FunctorType,
          bool Enable = (FunctorDeclaresValueType<FunctorType, void>::value) ||
                        (ReduceFunctorHasInit<FunctorType>::value) ||
                        (ReduceFunctorHasJoin<FunctorType>::value) ||
                        (ReduceFunctorHasFinal<FunctorType>::value) ||
                        (ReduceFunctorHasShmemSize<FunctorType>::value)>
struct IsNonTrivialReduceFunctor {
  enum : bool { value = false };
};

template <class FunctorType>
struct IsNonTrivialReduceFunctor<FunctorType, true> {
  enum : bool { value = true };
};

/** \brief  Query Functor and execution policy argument tag for value type.
 *
 *  If C++11 enabled and 'value_type' is not explicitly declared then attempt
 *  to deduce the type from FunctorType::operator().
 */
template <class FunctorType, class ArgTag,
          bool Dec = FunctorDeclaresValueType<FunctorType, ArgTag>::value>
struct FunctorValueTraits {
  using value_type     = void;
  using pointer_type   = void;
  using reference_type = void;
  using functor_type   = void;

  enum { StaticValueSize = 0 };

  KOKKOS_FORCEINLINE_FUNCTION static unsigned value_count(const FunctorType&) {
    return 0;
  }

  KOKKOS_FORCEINLINE_FUNCTION static unsigned value_size(const FunctorType&) {
    return 0;
  }
};

template <class ArgTag>
struct FunctorValueTraits<void, ArgTag, false> {
  using value_type     = void;
  using pointer_type   = void;
  using reference_type = void;
  using functor_type   = void;
};

/** \brief  FunctorType::value_type is explicitly declared so use it.
 *
 * Two options for declaration
 *
 *   1) A plain-old-data (POD) type
 *        using value_type = {pod_type};
 *
 *   2) An array of POD of a runtime specified count.
 *        using value_type = {pod_type}[];
 *        const unsigned     value_count ;
 */
template <class FunctorType, class ArgTag>
struct FunctorValueTraits<FunctorType, ArgTag,
                          true /* == exists FunctorType::value_type */> {
  using value_type =
      typename std::remove_extent<typename FunctorType::value_type>::type;
  using functor_type = FunctorType;

  static_assert((sizeof(value_type) < sizeof(int)) ||
                    0 == (sizeof(value_type) % sizeof(int)),
                "Reduction functor's declared value_type requires: 0 == "
                "sizeof(value_type) % sizeof(int)");

  /* this cast to bool is needed for correctness by NVCC */
  enum : bool {
    IsArray = static_cast<bool>(
        std::is_array<typename FunctorType::value_type>::value)
  };

  // If not an array then what is the sizeof(value_type)
  enum { StaticValueSize = IsArray ? 0 : sizeof(value_type) };

  using pointer_type = value_type*;

  // The reference_type for an array is 'value_type *'
  // The reference_type for a single value is 'value_type &'

  using reference_type =
      typename Impl::if_c<IsArray, value_type*, value_type&>::type;

  // Number of values if single value
  template <class F>
  KOKKOS_FORCEINLINE_FUNCTION static
      typename std::enable_if<std::is_same<F, FunctorType>::value && !IsArray,
                              unsigned>::type
      value_count(const F&) {
    return 1;
  }

  // Number of values if an array, protect via templating because
  // 'f.value_count' will only exist when the functor declares the value_type to
  // be an array.
  template <class F>
  KOKKOS_FORCEINLINE_FUNCTION static
      typename std::enable_if<std::is_same<F, FunctorType>::value && IsArray,
                              unsigned>::type
      value_count(const F& f) {
    return f.value_count;
  }

  // Total size of the value
  KOKKOS_INLINE_FUNCTION static unsigned value_size(const FunctorType& f) {
    return value_count(f) * sizeof(value_type);
  }
};

template <class FunctorType, class ArgTag>
struct FunctorValueTraits<FunctorType, ArgTag,
                          false /* == exists FunctorType::value_type */
                          > {
 private:
  struct VOIDTAG {
  };  // Allow declaration of non-matching operator() with void argument tag.
  struct REJECTTAG {
  };  // Reject tagged operator() when using non-tagged execution policy.

  using tag_type = typename Impl::if_c<std::is_same<ArgTag, void>::value,
                                       VOIDTAG, ArgTag>::type;

  //----------------------------------------
  // parallel_for operator without a tag:

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember)
                   const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember)
                   const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, const ArgMember&,
                                     const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, const ArgMember&,
                                     const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class TagType, class ArgMember>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  //----------------------------------------
  // parallel_for operator with a tag:

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, const ArgMember&,
                                      const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, const ArgMember&,
                                      const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&) const) {}

  template <class ArgMember>
  KOKKOS_INLINE_FUNCTION static VOIDTAG deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&) const) {}

  //----------------------------------------
  // parallel_reduce operator without a tag:
  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&, const ArgMember&,
                                     const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, T&) const) {
  }

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, const ArgMember&,
                                     const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, ArgMember,
                                     ArgMember, ArgMember, ArgMember, ArgMember,
                                     ArgMember, ArgMember, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, const ArgMember&,
                                     const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  //----------------------------------------
  // parallel_reduce operator with a tag:

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, ArgMember,
                                      ArgMember, ArgMember, ArgMember,
                                      ArgMember, ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, ArgMember, ArgMember,
                            ArgMember, ArgMember, ArgMember, ArgMember,
                            ArgMember, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, const ArgMember&,
                                      const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, const ArgMember&,
                                      const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, T&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&,
                            const ArgMember&, const ArgMember&, T&) const) {}

  //----------------------------------------
  // parallel_scan operator without a tag:

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, T&, bool) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, T&, bool) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, ArgMember, T&, bool) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, const ArgMember&, T&, bool) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(const TagType&, ArgMember, T&, bool) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, const ArgMember&, T&, bool)
                   const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(ArgMember, T&, const bool&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const ArgMember&, T&, const bool&) const) {
  }

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG,
      void (FunctorType::*)(TagType, ArgMember, T&, const bool&) const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(TagType, const ArgMember&, T&, const bool&)
                   const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, ArgMember, T&, const bool&)
                   const) {}

  template <class TagType, class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static REJECTTAG deduce_reduce_type(
      VOIDTAG, void (FunctorType::*)(const TagType&, const ArgMember&, T&,
                                     const bool&) const) {}
  //----------------------------------------
  // parallel_scan operator with a tag:

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, ArgMember, T&, bool) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(const tag_type&, ArgMember, T&, bool) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, const ArgMember&, T&, bool) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, const ArgMember&, T&,
                                      bool) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type,
      void (FunctorType::*)(tag_type, ArgMember, T&, const bool&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, ArgMember, T&,
                                      const bool&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(tag_type, const ArgMember&, T&,
                                      const bool&) const) {}

  template <class ArgMember, class T>
  KOKKOS_INLINE_FUNCTION static T deduce_reduce_type(
      tag_type, void (FunctorType::*)(const tag_type&, const ArgMember&, T&,
                                      const bool&) const) {}
  //----------------------------------------

  using ValueType =
      decltype(deduce_reduce_type(tag_type(), &FunctorType::operator()));

  enum { IS_VOID = std::is_same<VOIDTAG, ValueType>::value };
  enum { IS_REJECT = std::is_same<REJECTTAG, ValueType>::value };

 public:
  using value_type =
      typename Impl::if_c<IS_VOID || IS_REJECT, void, ValueType>::type;
  using pointer_type =
      typename Impl::if_c<IS_VOID || IS_REJECT, void, ValueType*>::type;
  using reference_type =
      typename Impl::if_c<IS_VOID || IS_REJECT, void, ValueType&>::type;
  using functor_type = FunctorType;

  static_assert(
      IS_VOID || IS_REJECT || 0 == (sizeof(ValueType) % sizeof(int)),
      "Reduction functor's value_type deduced from functor::operator() "
      "requires: 0 == sizeof(value_type) % sizeof(int)");

  enum { StaticValueSize = IS_VOID || IS_REJECT ? 0 : sizeof(ValueType) };

  KOKKOS_FORCEINLINE_FUNCTION static unsigned value_size(const FunctorType&) {
    return StaticValueSize;
  }

  KOKKOS_FORCEINLINE_FUNCTION static unsigned value_count(const FunctorType&) {
    return IS_VOID || IS_REJECT ? 0 : 1;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** Function signatures for FunctorType::init function with a tag.
 *  reference_type is 'value_type &' for scalar and 'value_type *' for array.
 */
template <class FunctorType, class ArgTag>
struct FunctorValueInitFunction {
  using reference_type =
      typename FunctorValueTraits<FunctorType, ArgTag>::reference_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, reference_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, reference_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag,
                                                        reference_type));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        reference_type));
};

/** Function signatures for FunctorType::init function without a tag.
 *  reference_type is 'value_type &' for scalar and 'value_type *' for array.
 */
template <class FunctorType>
struct FunctorValueInitFunction<FunctorType, void> {
  using reference_type =
      typename FunctorValueTraits<FunctorType, void>::reference_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(reference_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(reference_type));
};

// Adapter for value initialization function.
// If a proper FunctorType::init is declared then use it,
// otherwise use default constructor.
template <class FunctorType, class ArgTag,
          class T = typename FunctorValueTraits<FunctorType, ArgTag>::
              reference_type  // FIXME Fix FunctorValueTraits for multi-dim
                              // operator
          ,
          class Enable = void>
struct FunctorValueInit;

/* No 'init' function provided for single value */
template <class FunctorType, class ArgTag, class T, class Enable>
struct FunctorValueInit<FunctorType, ArgTag, T&, Enable> {
  KOKKOS_FORCEINLINE_FUNCTION static T& init(const FunctorType&, void* p) {
    return *(new (p) T());
  };
};

/* No 'init' function provided for array value */
template <class FunctorType, class ArgTag, class T, class Enable>
struct FunctorValueInit<FunctorType, ArgTag, T*, Enable> {
  KOKKOS_FORCEINLINE_FUNCTION static T* init(const FunctorType& f, void* p) {
    const int n = FunctorValueTraits<FunctorType, ArgTag>::value_count(f);
    for (int i = 0; i < n; ++i) {
      new (((T*)p) + i) T();
    }
    return (T*)p;
  }
};

/* 'init' function provided for single value */
template <class FunctorType, class T>
struct FunctorValueInit<
    FunctorType, void,
    T&
    // First  substitution failure when FunctorType::init does not exist.
    // Second substitution failure when FunctorType::init is not compatible.
    ,
    decltype(FunctorValueInitFunction<FunctorType, void>::enable_if(
        &FunctorType::init))> {
  KOKKOS_FORCEINLINE_FUNCTION static T& init(const FunctorType& f, void* p) {
    f.init(*((T*)p));
    return *((T*)p);
  }
};

/* 'init' function provided for array value */
template <class FunctorType, class T>
struct FunctorValueInit<
    FunctorType, void,
    T*
    // First  substitution failure when FunctorType::init does not exist.
    // Second substitution failure when FunctorType::init is not compatible
    ,
    decltype(FunctorValueInitFunction<FunctorType, void>::enable_if(
        &FunctorType::init))> {
  KOKKOS_FORCEINLINE_FUNCTION static T* init(const FunctorType& f, void* p) {
    f.init((T*)p);
    return (T*)p;
  }
};

/* 'init' function provided for single value */
template <class FunctorType, class ArgTag, class T>
struct FunctorValueInit<
    FunctorType, ArgTag,
    T&
    // First  substitution failure when FunctorType::init does not exist.
    // Second substitution failure when FunctorType::init is not compatible.
    ,
    typename std::enable_if<
        !std::is_same<ArgTag, void>::value,
        decltype(FunctorValueInitFunction<FunctorType, ArgTag>::enable_if(
            &FunctorType::init))>::type> {
  KOKKOS_FORCEINLINE_FUNCTION static T& init(const FunctorType& f, void* p) {
    f.init(ArgTag(), *((T*)p));
    return *((T*)p);
  }
};

/* 'init' function provided for array value */
template <class FunctorType, class ArgTag, class T>
struct FunctorValueInit<
    FunctorType, ArgTag,
    T*
    // First  substitution failure when FunctorType::init does not exist.
    // Second substitution failure when FunctorType::init is not compatible
    ,
    typename std::enable_if<
        !std::is_same<ArgTag, void>::value,
        decltype(FunctorValueInitFunction<FunctorType, ArgTag>::enable_if(
            &FunctorType::init))>::type> {
  KOKKOS_FORCEINLINE_FUNCTION static T* init(const FunctorType& f, void* p) {
    f.init(ArgTag(), (T*)p);
    return (T*)p;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Signatures for compatible FunctorType::join with tag and not an array
template <class FunctorType, class ArgTag,
          bool IsArray =
              0 == FunctorValueTraits<FunctorType, ArgTag>::StaticValueSize>
struct FunctorValueJoinFunction {
  using value_type =
      typename FunctorValueTraits<FunctorType, ArgTag>::value_type;

  using vref_type  = volatile value_type&;
  using cvref_type = const volatile value_type&;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, vref_type, cvref_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, vref_type, cvref_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag, vref_type,
                                                        cvref_type));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        vref_type, cvref_type));
};

// Signatures for compatible FunctorType::join with tag and is an array
template <class FunctorType, class ArgTag>
struct FunctorValueJoinFunction<FunctorType, ArgTag, true> {
  using value_type =
      typename FunctorValueTraits<FunctorType, ArgTag>::value_type;

  using vptr_type  = volatile value_type*;
  using cvptr_type = const volatile value_type*;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, vptr_type, cvptr_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, vptr_type, cvptr_type) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag, vptr_type,
                                                        cvptr_type));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        vptr_type, cvptr_type));
};

// Signatures for compatible FunctorType::join without tag and not an array
template <class FunctorType>
struct FunctorValueJoinFunction<FunctorType, void, false> {
  using value_type = typename FunctorValueTraits<FunctorType, void>::value_type;

  using vref_type  = volatile value_type&;
  using cvref_type = const volatile value_type&;

  KOKKOS_INLINE_FUNCTION static void enable_if(void (FunctorType::*)(vref_type,
                                                                     cvref_type)
                                                   const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(vref_type, cvref_type));
};

// Signatures for compatible FunctorType::join without tag and is an array
template <class FunctorType>
struct FunctorValueJoinFunction<FunctorType, void, true> {
  using value_type = typename FunctorValueTraits<FunctorType, void>::value_type;

  using vptr_type  = volatile value_type*;
  using cvptr_type = const volatile value_type*;

  KOKKOS_INLINE_FUNCTION static void enable_if(void (FunctorType::*)(vptr_type,
                                                                     cvptr_type)
                                                   const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(vptr_type, cvptr_type));
};

template <class FunctorType, class ArgTag,
          class T =
              typename FunctorValueTraits<FunctorType, ArgTag>::reference_type,
          class Enable = void>
struct FunctorValueJoin;

/* No 'join' function provided, single value */
template <class FunctorType, class ArgTag, class T, class Enable>
struct FunctorValueJoin<FunctorType, ArgTag, T&, Enable> {
  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType&) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& /*f*/,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    *((volatile T*)lhs) += *((const volatile T*)rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T& lhs, const volatile T& rhs) const { lhs += rhs; }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T& lhs, const T& rhs) const { lhs += rhs; }
};

/* No 'join' function provided, array of values */
template <class FunctorType, class ArgTag, class T, class Enable>
struct FunctorValueJoin<FunctorType, ArgTag, T*, Enable> {
  const FunctorType& f;

  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType& f_) : f(f_) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& f_,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    const int n = FunctorValueTraits<FunctorType, ArgTag>::value_count(f_);

    for (int i = 0; i < n; ++i) {
      ((volatile T*)lhs)[i] += ((const volatile T*)rhs)[i];
    }
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T* const lhs, const volatile T* const rhs) const {
    const int n = FunctorValueTraits<FunctorType, ArgTag>::value_count(f);

    for (int i = 0; i < n; ++i) {
      lhs[i] += rhs[i];
    }
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T* lhs, const T* rhs) const {
    const int n = FunctorValueTraits<FunctorType, ArgTag>::value_count(f);

    for (int i = 0; i < n; ++i) {
      lhs[i] += rhs[i];
    }
  }
};

/* 'join' function provided, single value */
template <class FunctorType, class ArgTag, class T>
struct FunctorValueJoin<
    FunctorType, ArgTag,
    T&
    // First  substitution failure when FunctorType::join does not exist.
    // Second substitution failure when enable_if( & Functor::join ) does not
    // exist
    ,
    decltype(FunctorValueJoinFunction<FunctorType, ArgTag>::enable_if(
        &FunctorType::join))> {
  const FunctorType& f;

  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType& f_) : f(f_) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& f_,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    f_.join(ArgTag(), *((volatile T*)lhs), *((const volatile T*)rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T& lhs, const volatile T& rhs) const {
    f.join(ArgTag(), lhs, rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T& lhs, const T& rhs) const { f.join(ArgTag(), lhs, rhs); }
};

/* 'join' function provided, no tag, single value */
template <class FunctorType, class T>
struct FunctorValueJoin<
    FunctorType, void,
    T&
    // First  substitution failure when FunctorType::join does not exist.
    // Second substitution failure when enable_if( & Functor::join ) does not
    // exist
    ,
    decltype(FunctorValueJoinFunction<FunctorType, void>::enable_if(
        &FunctorType::join))> {
  const FunctorType& f;

  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType& f_) : f(f_) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& f_,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    f_.join(*((volatile T*)lhs), *((const volatile T*)rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T& lhs, const volatile T& rhs) const {
    f.join(lhs, rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T& lhs, const T& rhs) const { f.join(lhs, rhs); }
};

/* 'join' function provided for array value */
template <class FunctorType, class ArgTag, class T>
struct FunctorValueJoin<
    FunctorType, ArgTag,
    T*
    // First  substitution failure when FunctorType::join does not exist.
    // Second substitution failure when enable_if( & Functor::join ) does not
    // exist
    ,
    decltype(FunctorValueJoinFunction<FunctorType, ArgTag>::enable_if(
        &FunctorType::join))> {
  const FunctorType& f;

  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType& f_) : f(f_) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& f_,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    f_.join(ArgTag(), (volatile T*)lhs, (const volatile T*)rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T* const lhs, const volatile T* const rhs) const {
    f.join(ArgTag(), lhs, rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T* lhs, const T* rhs) const { f.join(ArgTag(), lhs, rhs); }
};

/* 'join' function provided, no tag, array value */
template <class FunctorType, class T>
struct FunctorValueJoin<
    FunctorType, void,
    T*
    // First  substitution failure when FunctorType::join does not exist.
    // Second substitution failure when enable_if( & Functor::join ) does not
    // exist
    ,
    decltype(FunctorValueJoinFunction<FunctorType, void>::enable_if(
        &FunctorType::join))> {
  const FunctorType& f;

  KOKKOS_FORCEINLINE_FUNCTION
  FunctorValueJoin(const FunctorType& f_) : f(f_) {}

  KOKKOS_FORCEINLINE_FUNCTION static void join(const FunctorType& f_,
                                               volatile void* const lhs,
                                               const volatile void* const rhs) {
    f_.join((volatile T*)lhs, (const volatile T*)rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(volatile T* const lhs, const volatile T* const rhs) const {
    f.join(lhs, rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  void operator()(T* lhs, const T* rhs) const { f.join(lhs, rhs); }
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

namespace Impl {

template <typename ValueType, class JoinOp, class Enable = void>
struct JoinLambdaAdapter {
  using value_type = ValueType;
  const JoinOp& lambda;
  KOKKOS_INLINE_FUNCTION
  JoinLambdaAdapter(const JoinOp& lambda_) : lambda(lambda_) {}

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dst, const volatile value_type& src) const {
    lambda(dst, src);
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const { lambda(dst, src); }

  KOKKOS_INLINE_FUNCTION
  void operator()(volatile value_type& dst,
                  const volatile value_type& src) const {
    lambda(dst, src);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(value_type& dst, const value_type& src) const {
    lambda(dst, src);
  }
};

template <typename ValueType, class JoinOp>
struct JoinLambdaAdapter<ValueType, JoinOp,
                         decltype(FunctorValueJoinFunction<
                                  JoinOp, void>::enable_if(&JoinOp::join))> {
  using value_type = ValueType;
  static_assert(
      std::is_same<ValueType, typename JoinOp::value_type>::value,
      "JoinLambdaAdapter static_assert Fail: ValueType != JoinOp::value_type");

  const JoinOp& lambda;
  KOKKOS_INLINE_FUNCTION
  JoinLambdaAdapter(const JoinOp& lambda_) : lambda(lambda_) {}

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dst, const volatile value_type& src) const {
    lambda.join(dst, src);
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const {
    lambda.join(dst, src);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(volatile value_type& dst,
                  const volatile value_type& src) const {
    lambda.join(dst, src);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(value_type& dst, const value_type& src) const {
    lambda.join(dst, src);
  }
};

template <typename ValueType>
struct JoinAdd {
  using value_type = ValueType;

  KOKKOS_DEFAULTED_FUNCTION
  JoinAdd() = default;

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dst, const volatile value_type& src) const {
    dst += src;
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(value_type& dst, const value_type& src) const { dst += src; }
  KOKKOS_INLINE_FUNCTION
  void operator()(volatile value_type& dst,
                  const volatile value_type& src) const {
    dst += src;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class ArgTag,
          class T =
              typename FunctorValueTraits<FunctorType, ArgTag>::reference_type>
struct FunctorValueOps;

template <class FunctorType, class ArgTag, class T>
struct FunctorValueOps<FunctorType, ArgTag, T&> {
  KOKKOS_FORCEINLINE_FUNCTION static T* pointer(T& r) { return &r; }

  KOKKOS_FORCEINLINE_FUNCTION static T& reference(void* p) { return *((T*)p); }

  KOKKOS_FORCEINLINE_FUNCTION static void copy(const FunctorType&,
                                               void* const lhs,
                                               const void* const rhs) {
    *((T*)lhs) = *((const T*)rhs);
  }
};

/* No 'join' function provided, array of values */
template <class FunctorType, class ArgTag, class T>
struct FunctorValueOps<FunctorType, ArgTag, T*> {
  KOKKOS_FORCEINLINE_FUNCTION static T* pointer(T* p) { return p; }

  KOKKOS_FORCEINLINE_FUNCTION static T* reference(void* p) { return ((T*)p); }

  KOKKOS_FORCEINLINE_FUNCTION static void copy(const FunctorType& f,
                                               void* const lhs,
                                               const void* const rhs) {
    const int n = FunctorValueTraits<FunctorType, ArgTag>::value_count(f);
    for (int i = 0; i < n; ++i) {
      ((T*)lhs)[i] = ((const T*)rhs)[i];
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Compatible functions for 'final' function and value_type not an array
template <class FunctorType, class ArgTag,
          bool IsArray =
              0 == FunctorValueTraits<FunctorType, ArgTag>::StaticValueSize>
struct FunctorFinalFunction {
  using value_type =
      typename FunctorValueTraits<FunctorType, ArgTag>::value_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type&));
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag, value_type&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        value_type&));

  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // , value_type volatile & ) const ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (FunctorType::*)( ArgTag const & , value_type volatile & )
  // const ); KOKKOS_INLINE_FUNCTION static void enable_if( void
  // (FunctorType::*)( ArgTag         , value_type volatile & ) );
  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // const & , value_type volatile & ) ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (             *)( ArgTag         , value_type volatile & )
  // ); KOKKOS_INLINE_FUNCTION static void enable_if( void (             *)(
  // ArgTag const & , value_type volatile & ) );

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type const&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type const&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type const&));
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type const&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag,
                                                        value_type const&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        value_type const&));

  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // , value_type const volatile & ) const ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (FunctorType::*)( ArgTag const & , value_type const
  // volatile & ) const ); KOKKOS_INLINE_FUNCTION static void enable_if( void
  // (FunctorType::*)( ArgTag         , value_type const volatile & ) );
  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // const & , value_type const volatile & ) ); KOKKOS_INLINE_FUNCTION static
  // void enable_if( void (             *)( ArgTag         , value_type const
  // volatile & ) ); KOKKOS_INLINE_FUNCTION static void enable_if( void ( *)(
  // ArgTag const & , value_type const volatile & ) );
};

// Compatible functions for 'final' function and value_type is an array
template <class FunctorType, class ArgTag>
struct FunctorFinalFunction<FunctorType, ArgTag, true> {
  using value_type =
      typename FunctorValueTraits<FunctorType, ArgTag>::value_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type*));
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag, value_type*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        value_type*));

  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // , value_type volatile * ) const ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (FunctorType::*)( ArgTag const & , value_type volatile * )
  // const ); KOKKOS_INLINE_FUNCTION static void enable_if( void
  // (FunctorType::*)( ArgTag         , value_type volatile * ) );
  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // const & , value_type volatile * ) ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (             *)( ArgTag         , value_type volatile * )
  // ); KOKKOS_INLINE_FUNCTION static void enable_if( void (             *)(
  // ArgTag const & , value_type volatile * ) );

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type const*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type const*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, value_type const*));
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, value_type const*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag,
                                                        value_type const*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        value_type const*));

  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // , value_type const volatile * ) const ); KOKKOS_INLINE_FUNCTION static void
  // enable_if( void (FunctorType::*)( ArgTag const & , value_type const
  // volatile * ) const ); KOKKOS_INLINE_FUNCTION static void enable_if( void
  // (FunctorType::*)( ArgTag         , value_type const volatile * ) );
  // KOKKOS_INLINE_FUNCTION static void enable_if( void (FunctorType::*)( ArgTag
  // const & , value_type const volatile * ) ); KOKKOS_INLINE_FUNCTION static
  // void enable_if( void (             *)( ArgTag         , value_type const
  // volatile * ) ); KOKKOS_INLINE_FUNCTION static void enable_if( void ( *)(
  // ArgTag const & , value_type const volatile * ) );
};

template <class FunctorType>
struct FunctorFinalFunction<FunctorType, void, false> {
  using value_type = typename FunctorValueTraits<FunctorType, void>::value_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(value_type&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(value_type&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(value_type&));

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(const value_type&) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(const value_type&));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(const value_type&));
};

template <class FunctorType>
struct FunctorFinalFunction<FunctorType, void, true> {
  using value_type = typename FunctorValueTraits<FunctorType, void>::value_type;

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(value_type*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(value_type*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(value_type*));

  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(const value_type*) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(const value_type*));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(const value_type*));
};

/* No 'final' function provided */
template <class FunctorType, class ArgTag,
          class ResultType =
              typename FunctorValueTraits<FunctorType, ArgTag>::reference_type,
          class Enable = void>
struct FunctorFinal {
  KOKKOS_FORCEINLINE_FUNCTION static void final(const FunctorType&, void*) {}
};

/* 'final' function provided */
template <class FunctorType, class ArgTag, class T>
struct FunctorFinal<FunctorType, ArgTag,
                    T&
                    // First  substitution failure when FunctorType::final does
                    // not exist. Second substitution failure when enable_if( &
                    // Functor::final ) does not exist
                    ,
                    decltype(
                        FunctorFinalFunction<FunctorType, ArgTag>::enable_if(
                            &FunctorType::final))> {
  KOKKOS_FORCEINLINE_FUNCTION static void final(const FunctorType& f, void* p) {
    f.final(*((T*)p));
  }

  KOKKOS_FORCEINLINE_FUNCTION static void final(FunctorType& f, void* p) {
    f.final(*((T*)p));
  }
};

/* 'final' function provided for array value */
template <class FunctorType, class ArgTag, class T>
struct FunctorFinal<FunctorType, ArgTag,
                    T*
                    // First  substitution failure when FunctorType::final does
                    // not exist. Second substitution failure when enable_if( &
                    // Functor::final ) does not exist
                    ,
                    decltype(
                        FunctorFinalFunction<FunctorType, ArgTag>::enable_if(
                            &FunctorType::final))> {
  KOKKOS_FORCEINLINE_FUNCTION static void final(const FunctorType& f, void* p) {
    f.final((T*)p);
  }

  KOKKOS_FORCEINLINE_FUNCTION static void final(FunctorType& f, void* p) {
    f.final((T*)p);
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class ArgTag,
          class ReferenceType =
              typename FunctorValueTraits<FunctorType, ArgTag>::reference_type>
struct FunctorApplyFunction {
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, ReferenceType) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, ReferenceType) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag, ReferenceType));
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ArgTag const&, ReferenceType));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag, ReferenceType));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ArgTag const&,
                                                        ReferenceType));
};

template <class FunctorType, class ReferenceType>
struct FunctorApplyFunction<FunctorType, void, ReferenceType> {
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ReferenceType) const);
  KOKKOS_INLINE_FUNCTION static void enable_if(
      void (FunctorType::*)(ReferenceType));
  KOKKOS_INLINE_FUNCTION static void enable_if(void (*)(ReferenceType));
};

template <class FunctorType>
struct FunctorApplyFunction<FunctorType, void, void> {
  KOKKOS_INLINE_FUNCTION static void enable_if(void (FunctorType::*)() const);
  KOKKOS_INLINE_FUNCTION static void enable_if(void (FunctorType::*)());
};

template <class FunctorType, class ArgTag, class ReferenceType,
          class Enable = void>
struct FunctorApply {
  KOKKOS_FORCEINLINE_FUNCTION static void apply(const FunctorType&, void*) {}
};

/* 'apply' function provided for void value */
template <class FunctorType, class ArgTag>
struct FunctorApply<
    FunctorType, ArgTag,
    void
    // First  substitution failure when FunctorType::apply does not exist.
    // Second substitution failure when enable_if( & Functor::apply ) does not
    // exist
    ,
    decltype(FunctorApplyFunction<FunctorType, ArgTag, void>::enable_if(
        &FunctorType::apply))> {
  KOKKOS_FORCEINLINE_FUNCTION static void apply(FunctorType& f) { f.apply(); }

  KOKKOS_FORCEINLINE_FUNCTION static void apply(const FunctorType& f) {
    f.apply();
  }
};

/* 'apply' function provided for single value */
template <class FunctorType, class ArgTag, class T>
struct FunctorApply<FunctorType, ArgTag,
                    T&
                    // First  substitution failure when FunctorType::apply does
                    // not exist. Second substitution failure when enable_if( &
                    // Functor::apply ) does not exist
                    ,
                    decltype(
                        FunctorApplyFunction<FunctorType, ArgTag>::enable_if(
                            &FunctorType::apply))> {
  KOKKOS_FORCEINLINE_FUNCTION static void apply(const FunctorType& f, void* p) {
    f.apply(*((T*)p));
  }

  KOKKOS_FORCEINLINE_FUNCTION static void apply(FunctorType& f, void* p) {
    f.apply(*((T*)p));
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_FUNCTORADAPTER_HPP */
