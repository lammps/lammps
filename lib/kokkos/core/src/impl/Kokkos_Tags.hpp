/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_TAGS_HPP
#define KOKKOS_TAGS_HPP

#include <impl/Kokkos_Traits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct LayoutTag {};

struct MemorySpaceTag {};
struct MemoryTraitsTag {};

struct ExecutionPolicyTag {};
struct ExecutionSpaceTag {};


template< class C , class Enable = void >
struct is_memory_space : public bool_< false > {};

template< class C >
struct is_memory_space< C , typename Impl::enable_if_type< typename C::kokkos_tag >::type >
  : public bool_< Impl::is_same< typename C::kokkos_tag , Impl::MemorySpaceTag >::value > {};


template< class C , class Enable = void >
struct is_execution_space : public bool_< false > {};

template< class C >
struct is_execution_space< C , typename Impl::enable_if_type< typename C::kokkos_tag >::type >
  : public bool_< Impl::is_same< typename C::kokkos_tag , Impl::ExecutionSpaceTag >::value > {};


template< class C , class Enable = void >
struct is_execution_policy : public bool_< false > {};

template< class C >
struct is_execution_policy< C , typename Impl::enable_if_type< typename C::kokkos_tag >::type >
  : public bool_< Impl::is_same< typename C::kokkos_tag , Impl::ExecutionPolicyTag >::value > {};


template< class C , class Enable = void >
struct is_layout : public Impl::false_type {};

template<class C>
struct is_layout<C,typename Impl::enable_if_type< typename C::kokkos_tag >::type > {
  enum {value=bool(Impl::is_same<Impl::LayoutTag,typename C::kokkos_tag>::value)};
};


template< class C , class Enable = void >
struct is_memorytraits : public Impl::false_type {};

template<class C>
struct is_memorytraits<C,typename Impl::enable_if_type< typename C::kokkos_tag >::type > {
  enum {value=bool(Impl::is_same<Impl::MemoryTraitsTag,typename C::kokkos_tag>::value)};
};

}
}

#endif
