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
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Tools_Generic.hpp>
using ExecSpace  = Kokkos::DefaultHostExecutionSpace;
using TeamMember = Kokkos::TeamPolicy<ExecSpace>::member_type;
struct TestTeamFunctor {
  KOKKOS_FUNCTION void operator()(TeamMember) const {}
};
struct TestMDFunctor {
  KOKKOS_FUNCTION void operator()(const int, const int) const {}
};
int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    Kokkos::TeamPolicy<ExecSpace> teamp(1, Kokkos::AUTO, Kokkos::AUTO);
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdp({0, 0}, {1, 1});
    Kokkos::Tools::Experimental::TeamSizeTuner team_tune_this(
        "team_tuner", teamp, TestTeamFunctor{}, Kokkos::ParallelForTag{},
        Kokkos::Tools::Experimental::Impl::Impl::SimpleTeamSizeCalculator{});

    Kokkos::Tools::Experimental::MDRangeTuner<2> md_tune_this(
        "md_tuner", mdp, TestMDFunctor{}, Kokkos::ParallelForTag{},
        Kokkos::Tools::Experimental::Impl::Impl::SimpleTeamSizeCalculator{});

    std::vector<int> options{1, 2, 3, 4, 5};

    auto new_team_tuner = team_tune_this.combine("options", options);
    auto new_md_tuner   = md_tune_this.combine("options", options);
    using namespace Kokkos::Tools::Experimental;
    VariableInfo info;
    info.category      = StatisticalCategory::kokkos_value_categorical;
    info.valueQuantity = CandidateValueType::kokkos_value_unbounded;
    info.type          = ValueType::kokkos_value_string;
    size_t input       = declare_input_type("kernel", info);
    VariableValue team_kernel_value = make_variable_value(input, "abs");
    VariableValue md_kernel_value   = make_variable_value(input, "abs");
    size_t kernel_context           = get_new_context_id();
    begin_context(kernel_context);
    set_input_values(kernel_context, 1, &team_kernel_value);
    for (int x = 0; x < 10000; ++x) {
      auto config = new_md_tuner.begin();
      int option  = std::get<0>(config);
      (void)option;
      int tile_x = std::get<1>(config);
      int tile_y = std::get<2>(config);
      Kokkos::parallel_for("mdrange",
                           Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
                               {0, 0}, {1, 1}, {tile_x, tile_y}),
                           TestMDFunctor{});
      new_md_tuner.end();
    }
    end_context(kernel_context);
    begin_context(kernel_context);
    set_input_values(kernel_context, 1, &md_kernel_value);

    /**
     * Note that 0.0 is basically a floating point index into
     * the outermost index in this, which is the options vector
     * above. The At 0.0, this will be the first element (1).
     * At 0.9 this will be the last element (5)
     */
    auto begin_point = new_team_tuner.get_point(0.0, 0.0, 0.0);
    assert(std::get<0>(begin_point) == 1);
    (void)begin_point;  // to avoid warnings in some compilers
    auto end_point = new_team_tuner.get_point(0.9, 0.0, 0.0);
    (void)end_point;  // to avoid warnings in some compilers
    assert(std::get<0>(end_point) == 5);
    for (int x = 0; x < 10000; ++x) {
      auto config = new_team_tuner.begin();
      int option  = std::get<0>(config);
      (void)option;
      int team   = std::get<1>(config);
      int vector = std::get<2>(config);
      Kokkos::parallel_for("mdrange",
                           Kokkos::TeamPolicy<ExecSpace>(1, team, vector),
                           TestTeamFunctor{});
      new_team_tuner.end();
    }
    end_context(kernel_context);
  }
  Kokkos::finalize();
}
