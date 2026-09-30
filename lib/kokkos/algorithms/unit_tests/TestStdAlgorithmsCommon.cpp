// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
