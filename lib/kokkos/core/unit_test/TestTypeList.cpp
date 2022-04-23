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

#include <impl/Kokkos_Utilities.hpp>

using TypeList2 = Kokkos::Impl::type_list<void, bool>;
using TypeList3 = Kokkos::Impl::type_list<char, short, int>;
using TypeList223 =
    Kokkos::Impl::type_list<void, bool, void, bool, char, short, int>;
using TypeList223Void   = Kokkos::Impl::type_list<void, void>;
using TypeList223NoVoid = Kokkos::Impl::type_list<bool, bool, char, short, int>;

// concat_type_list
using ConcatTypeList2 = Kokkos::Impl::concat_type_list_t<TypeList2>;
static_assert(std::is_same<TypeList2, ConcatTypeList2>::value,
              "concat_type_list of a single type_list failed");

using ConcatTypeList223 =
    Kokkos::Impl::concat_type_list_t<TypeList2, TypeList2, TypeList3>;
static_assert(std::is_same<TypeList223, ConcatTypeList223>::value,
              "concat_type_list of three type_lists failed");

// filter_type_list
using FilterTypeList223Void =
    Kokkos::Impl::filter_type_list_t<std::is_void, TypeList223>;
static_assert(std::is_same<TypeList223Void, FilterTypeList223Void>::value,
              "filter_type_list with predicate value==true failed");

using FilterTypeList223NoVoid =
    Kokkos::Impl::filter_type_list_t<std::is_void, TypeList223, false>;
static_assert(std::is_same<TypeList223NoVoid, FilterTypeList223NoVoid>::value,
              "filter_type_list with predicate value==false failed");
