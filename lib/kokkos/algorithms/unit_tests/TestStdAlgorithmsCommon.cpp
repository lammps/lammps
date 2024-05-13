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

#include <TestStdAlgorithmsCommon.hpp>

namespace Test {
namespace stdalgos {

std::string view_tag_to_string(DynamicTag) { return "dynamic_view"; }

std::string view_tag_to_string(DynamicLayoutLeftTag) {
  return "dynamic_layout_left_view";
}

std::string view_tag_to_string(DynamicLayoutRightTag) {
  return "dynamic_layout_right_view";
}

std::string view_tag_to_string(StridedTwoTag) { return "stride2_view"; }

std::string view_tag_to_string(StridedThreeTag) { return "stride3_view"; }

std::string view_tag_to_string(StridedTwoRowsTag) { return "stride2rows_view"; }

std::string view_tag_to_string(StridedThreeRowsTag) {
  return "stride3rows_view";
}

}  // namespace stdalgos
}  // namespace Test
