//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
