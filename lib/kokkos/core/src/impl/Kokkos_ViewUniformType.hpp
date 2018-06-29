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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXPERIMENTAL_VIEWUNIFORMTYPE_HPP
#define KOKKOS_EXPERIMENTAL_VIEWUNIFORMTYPE_HPP

namespace Kokkos {
namespace Impl {
  template< class ScalarType, int Rank>
  struct ViewScalarToDataType {
    typedef typename ViewScalarToDataType<ScalarType,Rank-1>::type* type;
  };

  template< class ScalarType>
  struct ViewScalarToDataType<ScalarType,0> {
    typedef ScalarType type;
  };

  template< class LayoutType, int Rank>
  struct ViewUniformLayout {
    typedef LayoutType array_layout;
  };

  template< class LayoutType>
  struct ViewUniformLayout<LayoutType, 0> {
    typedef Kokkos::LayoutLeft array_layout;
  };

  template<>
  struct ViewUniformLayout<Kokkos::LayoutRight, 1> {
    typedef Kokkos::LayoutLeft array_layout;
  };

  template< class ViewType , int Traits>
  struct ViewUniformType {
    typedef typename ViewType::data_type data_type;
    typedef typename std::add_const<typename ViewType::data_type>::type const_data_type;
    typedef typename ViewScalarToDataType<typename ViewType::value_type,ViewType::rank>::type runtime_data_type;
    typedef typename ViewScalarToDataType<typename std::add_const<typename ViewType::value_type>::type,ViewType::rank>::type runtime_const_data_type;

    typedef typename ViewUniformLayout<typename ViewType::array_layout, ViewType::rank>::array_layout array_layout;

    typedef typename ViewType::device_type device_type;
    typedef typename Kokkos::Device<typename device_type::execution_space,Kokkos::AnonymousSpace> anonymous_device_type;

    typedef typename Kokkos::MemoryTraits<Traits> memory_traits;
    typedef Kokkos::View<data_type,array_layout,device_type,memory_traits> type;
    typedef Kokkos::View<const_data_type,array_layout,device_type,memory_traits> const_type;
    typedef Kokkos::View<runtime_data_type,array_layout,device_type,memory_traits> runtime_type;
    typedef Kokkos::View<runtime_const_data_type,array_layout,device_type,memory_traits> runtime_const_type;

    typedef Kokkos::View<data_type,array_layout,anonymous_device_type,memory_traits> nomemspace_type;
    typedef Kokkos::View<const_data_type,array_layout,anonymous_device_type,memory_traits> const_nomemspace_type;
    typedef Kokkos::View<runtime_data_type,array_layout,anonymous_device_type,memory_traits> runtime_nomemspace_type;
    typedef Kokkos::View<runtime_const_data_type,array_layout,anonymous_device_type,memory_traits> runtime_const_nomemspace_type;
  };
}
}


#endif
