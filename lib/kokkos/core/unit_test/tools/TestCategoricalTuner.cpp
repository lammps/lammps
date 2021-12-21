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

// This file tests the categorical tuner

#include <Kokkos_Core.hpp>
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
