/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_TAGS_HPP
#define KOKKOS_TAGS_HPP

#include <impl/Kokkos_Traits.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <type_traits>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** KOKKOS_HAVE_TYPE( Type )
 *
 * defines a meta-function that check if a type expose an internal typedef or
 * type alias which matches Type
 *
 * e.g.
 *   KOKKOS_HAVE_TYPE( array_layout );
 *   struct Foo { using array_layout = void; };
 *   have_array_layout<Foo>::value == 1;
 */
#define KOKKOS_HAVE_TYPE( Type )                                                \
template <typename T>                                                           \
struct have_##Type {                                                            \
  template <typename U> static std::false_type have_type(...);                  \
  template <typename U> static std::true_type  have_type( typename U::Type* );  \
  using type = decltype(have_type<T>(nullptr));                                 \
  static constexpr bool value = type::value;                                    \
}

/** KOKKOS_IS_CONCEPT( Concept )
 *
 * defines a meta-function that check if a type match the given Kokkos concept
 * type alias which matches Type
 *
 * e.g.
 *   KOKKOS_IS_CONCEPT( array_layout );
 *   struct Foo { using array_layout = Foo; };
 *   is_array_layout<Foo>::value == 1;
 */
#define KOKKOS_IS_CONCEPT( Concept )                                            \
template <typename T>                                                           \
struct is_##Concept {                                                           \
  template <typename U> static std::false_type have_concept(...);               \
  template <typename U> static auto have_concept( typename U::Concept* )        \
                          ->typename std::is_same<T, typename U::Concept>::type;\
  using type = decltype(have_concept<T>(nullptr));                              \
  static constexpr bool value = type::value;                                    \
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos { namespace Impl {

template <typename T>
using is_void = std::is_same<void,T>;

// is_memory_space<T>::value
KOKKOS_IS_CONCEPT( memory_space );

// is_memory_traits<T>::value
KOKKOS_IS_CONCEPT( memory_traits );

// is_execution_space<T>::value
KOKKOS_IS_CONCEPT( execution_space );

// is_execution_policy<T>::value
KOKKOS_IS_CONCEPT( execution_policy );

// is_array_layout<T>::value
KOKKOS_IS_CONCEPT( array_layout );

// is_iteration_pattern<T>::value
KOKKOS_IS_CONCEPT( iteration_pattern );

// is_schedule_type<T>::value
KOKKOS_IS_CONCEPT( schedule_type );

// is_index_type<T>::value
KOKKOS_IS_CONCEPT( index_type );

}} // namespace Kokkos::Impl


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class ExecutionSpace , class MemorySpace >
struct Device {
  static_assert( Impl::is_execution_space<ExecutionSpace>::value
               , "Execution space is not valid" );
  static_assert( Impl::is_memory_space<MemorySpace>::value
               , "Memory space is not valid" );
  typedef ExecutionSpace execution_space;
  typedef MemorySpace memory_space;
  typedef Device<execution_space,memory_space> device_type;
};
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class C , class Enable = void >
struct is_space : public Impl::false_type {};

template< class C >
struct is_space< C
                 , typename Impl::enable_if<(
                     Impl::is_same< C , typename C::execution_space >::value ||
                     Impl::is_same< C , typename C::memory_space    >::value ||
                     Impl::is_same< C , Device<
                                             typename C::execution_space,
                                             typename C::memory_space> >::value
                   )>::type
                 >
  : public Impl::true_type
{
  typedef typename C::execution_space  execution_space ;
  typedef typename C::memory_space     memory_space ;

  // The host_memory_space defines a space with host-resident memory.
  // If the execution space's memory space is host accessible then use that execution space.
  // else use the HostSpace.
  typedef
      typename Impl::if_c< Impl::is_same< memory_space , HostSpace >::value
#ifdef KOKKOS_HAVE_CUDA
                        || Impl::is_same< memory_space , CudaUVMSpace>::value
                        || Impl::is_same< memory_space , CudaHostPinnedSpace>::value
#endif
                          , memory_space , HostSpace >::type
      host_memory_space ;

  // The host_execution_space defines a space which has access to HostSpace.
  // If the execution space can access HostSpace then use that execution space.
  // else use the DefaultHostExecutionSpace.
#ifdef KOKKOS_HAVE_CUDA
  typedef
      typename Impl::if_c< Impl::is_same< execution_space , Cuda >::value
                          , DefaultHostExecutionSpace , execution_space >::type
      host_execution_space ;
#else
  typedef execution_space host_execution_space;
#endif

  typedef Device<host_execution_space,host_memory_space> host_mirror_space;
};
}
}

#endif
