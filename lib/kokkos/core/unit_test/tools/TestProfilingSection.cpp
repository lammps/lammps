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

#include <Kokkos_Profiling_ProfileSection.hpp>

#include <gtest/gtest.h>

namespace {
struct Section {
  std::string name;
  int start_call_cnt;
  int stop_call_cnt;
  int is_destroyed;
  friend std::ostream& operator<<(std::ostream& os, Section const& s) {
    os << "( " << s.name << ", " << s.start_call_cnt << ", " << s.stop_call_cnt
       << ", " << s.is_destroyed << " )";
    return os;
  }
  friend bool operator==(Section const& l, Section const& r) {
    return (l.name == r.name) && (l.start_call_cnt == r.start_call_cnt) &&
           (l.stop_call_cnt == r.stop_call_cnt) &&
           (l.is_destroyed == r.is_destroyed);
  }
};

std::vector<Section> kokkosp_test_section_vector;

void kokkosp_test_create_section(char const* label, std::uint32_t* id) {
  *id = kokkosp_test_section_vector.size();
  kokkosp_test_section_vector.emplace_back(Section{label, 0, 0, 0});
}

void kokkosp_test_start_section(std::uint32_t id) {
  ++kokkosp_test_section_vector[id].start_call_cnt;
}

void kokkosp_test_stop_section(std::uint32_t id) {
  ++kokkosp_test_section_vector[id].stop_call_cnt;
}

void kokkosp_test_destroy_section(std::uint32_t id) {
  ++kokkosp_test_section_vector[id].is_destroyed;
}

}  // namespace

TEST(defaultdevicetype, profiling_section) {
  Kokkos::Profiling::Experimental::set_create_profile_section_callback(
      kokkosp_test_create_section);
  Kokkos::Profiling::Experimental::set_destroy_profile_section_callback(
      kokkosp_test_destroy_section);
  Kokkos::Profiling::Experimental::set_start_profile_section_callback(
      kokkosp_test_start_section);
  Kokkos::Profiling::Experimental::set_stop_profile_section_callback(
      kokkosp_test_stop_section);

  ASSERT_TRUE(kokkosp_test_section_vector.empty());

  {
    Kokkos::Profiling::ProfilingSection profile_1("one");
    ASSERT_EQ(kokkosp_test_section_vector.size(), 1u);
    ASSERT_EQ(kokkosp_test_section_vector[0], (Section{"one", 0, 0, 0}));

    // NOTE: ProfilingSection is a wrapper that manages the lifetime of the
    // underlying section but does not care whether the start and stop call
    // sequence makes any sense.
    profile_1.stop();
    profile_1.stop();
    profile_1.start();
    profile_1.start();
    profile_1.start();
    ASSERT_EQ(kokkosp_test_section_vector[0], (Section{"one", 3, 2, 0}));

    {
      Kokkos::Profiling::ProfilingSection profile_2("two");
      profile_2.start();
    }
    ASSERT_EQ(kokkosp_test_section_vector.size(), 2u);
    ASSERT_EQ(kokkosp_test_section_vector[1], (Section{"two", 1, 0, 1}));

    profile_1.start();
    profile_1.start();
  }

  ASSERT_EQ(kokkosp_test_section_vector.size(), 2u);
  ASSERT_EQ(kokkosp_test_section_vector[0], (Section{"one", 5, 2, 1}));
  ASSERT_EQ(kokkosp_test_section_vector[1], (Section{"two", 1, 0, 1}));

  // Cleanup
  kokkosp_test_section_vector.clear();
  Kokkos::Tools::Experimental::set_create_profile_section_callback(nullptr);
  Kokkos::Tools::Experimental::set_destroy_profile_section_callback(nullptr);
  Kokkos::Tools::Experimental::set_start_profile_section_callback(nullptr);
  Kokkos::Tools::Experimental::set_stop_profile_section_callback(nullptr);
}

using Kokkos::Profiling::ProfilingSection;
static_assert(!std::is_default_constructible<ProfilingSection>::value, "");
static_assert(!std::is_copy_constructible<ProfilingSection>::value, "");
static_assert(!std::is_move_constructible<ProfilingSection>::value, "");
static_assert(!std::is_copy_assignable<ProfilingSection>::value, "");
static_assert(!std::is_move_assignable<ProfilingSection>::value, "");
