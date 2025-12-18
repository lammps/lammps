// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// This file tests the categorical tuner

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <unistd.h>
struct point {
  float x;
  float y;
  float z;
};
void do_computation(const point& test_point) {
  usleep(((unsigned int)test_point.x) * 100);
}
using namespace Kokkos::Tools::Experimental;
int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    VariableInfo info;
    info.category              = StatisticalCategory::kokkos_value_categorical;
    info.valueQuantity         = CandidateValueType::kokkos_value_unbounded;
    info.type                  = ValueType::kokkos_value_string;
    size_t input               = declare_input_type("kernel", info);
    VariableValue kernel_value = make_variable_value(input, "abs");
    size_t kernel_context      = get_new_context_id();
    begin_context(kernel_context);
    set_input_values(kernel_context, 1, &kernel_value);

    std::vector<point> points;
    points.push_back({1.0, 1.0, 1.0});
    points.push_back({10.0, 10.0, 10.0});
    points.push_back({0.0, 0.0, 0.0});
    auto tuner =
        Kokkos::Tools::Experimental::make_categorical_tuner("points", points);
    for (decltype(points)::size_type x = 0; x < 3000; ++x) {
      point test_point = tuner.begin();
      do_computation(test_point);
      tuner.end();
    }

    end_context(kernel_context);
  }
  Kokkos::finalize();
}
