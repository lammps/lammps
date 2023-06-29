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
#include <utility>

namespace Test {
namespace stdalgos {
namespace Search_n {

namespace KE = Kokkos::Experimental;

// search_n is only available from c++20, so I have to put it here
template <class ForwardIt, class Size, class T, class BinaryPredicate>
ForwardIt my_std_search_n(ForwardIt first, ForwardIt last, Size count,
                          const T& value, BinaryPredicate p) {
  if (count <= 0) {
    return first;
  }
  for (; first != last; ++first) {
    if (!p(*first, value)) {
      continue;
    }

    ForwardIt candidate = first;
    Size cur_count      = 0;

    while (true) {
      ++cur_count;
      if (cur_count >= count) {
        // success
        return candidate;
      }
      ++first;
      if (first == last) {
        // exhausted the list
        return last;
      }
      if (!p(*first, value)) {
        // too few in a row
        break;
      }
    }
  }

  return last;
}

template <class ForwardIt, class Size, class T>
ForwardIt my_std_search_n(ForwardIt first, ForwardIt last, Size count,
                          const T& value) {
  using iter_value_type = typename ForwardIt::value_type;
  using p_type          = IsEqualFunctor<iter_value_type, T>;
  return my_std_search_n(first, last, count, value, p_type());
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(0, 20) { m_gen.seed(1034343); }
  UnifDist(int a, int b) : m_dist(a, b) { m_gen.seed(234343); }

  int operator()() { return m_dist(m_gen); }
};

template <class ViewType, class ValueType>
void fill_view(ViewType dest_view, ValueType value, std::size_t count,
               const std::string& where_to_place_count_values) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  // fill with something
  for (std::size_t i = 0; i < ext; ++i) {
    v_h(i) = 15;
  }

  if (where_to_place_count_values == "none") {
    // do nothing
  }

  else if (where_to_place_count_values == "left") {
    for (std::size_t i = 0; i < count; ++i) {
      v_h(i) = value;
    }
  }

  else if (where_to_place_count_values == "left_and_1567") {
    for (std::size_t i = 0; i < count; ++i) {
      v_h(i) = value;
    }

    for (std::size_t i = 0; i < count; ++i) {
      v_h(1567 + i) = value;
    }
  }

  else if (where_to_place_count_values == "random") {
    // find a random location to start filling view
    using dist_type = std::uniform_int_distribution<int>;
    std::random_device r;
    // from this:
    // https://stackoverflow.com/questions/34490599/c11-how-to-set-seed-using-random
    std::seed_seq seed{r(), r(), r(), r(), r(), r()};
    std::mt19937 gen(seed);
    dist_type dist(0, ext - count);
    const auto start_at = dist(gen);
    // std::cout << "start_at " << start_at << std::endl;

    for (std::size_t i = 0; i < count; ++i) {
      v_h(start_at + i) = value;
    }
  }

  else if (where_to_place_count_values == "11133_and_right") {
    for (std::size_t i = 0; i < count; ++i) {
      v_h(11133 + i) = value;
    }

    for (std::size_t i = 0; i < count; ++i) {
      v_h(ext - count + i) = value;
    }
  }

  else if (where_to_place_count_values == "right") {
    for (std::size_t i = 0; i < count; ++i) {
      v_h(ext - count + i) = value;
    }
  }

  else {
    throw std::runtime_error("Kokkos: test: search_n: this should not happen");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class Tag, class ValueType>
void print_scenario_details(const std::string& name, std::size_t count,
                            const std::string& where_to_place_count_values) {
  std::cout << "search_n: default predicate: " << name << ", "
            << "count = " << count << ", " << where_to_place_count_values
            << ", " << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ValueType, class Predicate>
void print_scenario_details(const std::string& name, std::size_t count,
                            const std::string& where_to_place_count_values,
                            Predicate pred) {
  (void)pred;
  std::cout << "search_n: custom predicate: " << name << ", "
            << "count = " << count << ", " << where_to_place_count_values
            << ", " << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ValueType, class InfoType, class... Args>
void run_single_scenario(const InfoType& scenario_info, std::size_t count,
                         ValueType value, Args... args) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  const auto count_place     = std::get<2>(scenario_info);
  // print_scenario_details<Tag, ValueType>(name, count, count_place, args...);

  auto view = create_view<ValueType>(Tag{}, view_ext, "search_n_test_view");
  fill_view(view, value, count, count_place);

  // run std
  auto view_h = create_host_space_copy(view);
  auto stdrit = my_std_search_n(KE::cbegin(view_h), KE::cend(view_h), count,
                                value, args...);
  const auto stddiff = stdrit - KE::cbegin(view_h);

  {
    auto myrit = KE::search_n(exespace(), KE::cbegin(view), KE::cend(view),
                              count, value, args...);
    const auto mydiff = myrit - KE::cbegin(view);
    ASSERT_EQ(mydiff, stddiff);
  }

  {
    auto myrit        = KE::search_n("label", exespace(), KE::cbegin(view),
                              KE::cend(view), count, value, args...);
    const auto mydiff = myrit - KE::cbegin(view);
    ASSERT_EQ(mydiff, stddiff);
  }

  {
    auto myrit = KE::search_n("label", exespace(), view, count, value, args...);
    const auto mydiff = myrit - KE::begin(view);
    ASSERT_EQ(mydiff, stddiff);
  }

  {
    auto myrit        = KE::search_n(exespace(), view, count, value, args...);
    const auto mydiff = myrit - KE::begin(view);
    ASSERT_EQ(mydiff, stddiff);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  using scenario_t = std::tuple<std::string, std::size_t, std::string>;
  std::vector<scenario_t> scenarios(31);
  scenarios[0] = scenario_t("empty", 0u, "none");
  scenarios[1] = scenario_t("one-element-a", 1u, "none");
  scenarios[2] = scenario_t("one-element-b", 1u, "left");

  scenarios[3] = scenario_t("two-elements-a", 2u, "none");
  scenarios[4] = scenario_t("two-elements-b", 2u, "left");
  scenarios[5] = scenario_t("two-elements-c", 2u, "right");

  scenarios[6] = scenario_t("three-elements-a", 3u, "none");
  scenarios[7] = scenario_t("three-elements-b", 3u, "left");
  scenarios[8] = scenario_t("three-elements-c", 3u, "random");
  scenarios[9] = scenario_t("three-elements-d", 3u, "right");

  scenarios[10] = scenario_t("four-elements-a", 4u, "none");
  scenarios[11] = scenario_t("four-elements-b", 4u, "left");
  scenarios[12] = scenario_t("four-elements-c", 4u, "random");
  scenarios[13] = scenario_t("four-elements-d", 4u, "right");

  scenarios[14] = scenario_t("small-a", 13u, "none");
  scenarios[15] = scenario_t("small-b", 13u, "left");
  scenarios[16] = scenario_t("small-c", 13u, "random");
  scenarios[17] = scenario_t("small-d", 13u, "right");
  scenarios[18] = scenario_t("small-e", 131u, "none");
  scenarios[19] = scenario_t("small-f", 131u, "left");
  scenarios[20] = scenario_t("small-g", 131u, "random");
  scenarios[21] = scenario_t("small-h", 131u, "right");

  scenarios[22] = scenario_t("medium-a", 21103u, "none");
  scenarios[22] = scenario_t("medium-b", 21103u, "left");
  scenarios[23] = scenario_t("medium-c", 21103u, "random");
  scenarios[24] = scenario_t("medium-d", 21103u, "right");
  scenarios[25] = scenario_t("medium-e", 21103u, "left_and_1567");
  scenarios[26] = scenario_t("medium-f", 21103u, "11133_and_right");

  scenarios[27] = scenario_t("large-a", 101333u, "none");
  scenarios[28] = scenario_t("large-b", 101333u, "left");
  scenarios[29] = scenario_t("large-c", 101333u, "random");
  scenarios[30] = scenario_t("large-d", 101333u, "right");

  const std::vector<std::size_t> counts = {1,  2,  3,  4,   5,  8,
                                           11, 13, 31, 131, 523};

  const ValueType target_value = 3;

  // for each view scenario, run "search_n" for multiple counts
  for (const auto& it : scenarios) {
    const std::size_t view_ext = std::get<1>(it);

    if (view_ext == 0) {
      run_single_scenario<Tag, ValueType>(it, 0, target_value);
    } else {
      for (const auto& it2 : counts) {
        // only run if view is larger or equal than count
        if (view_ext >= it2) {
          run_single_scenario<Tag, ValueType>(it, it2, target_value);

          using func_t = IsEqualFunctor<ValueType>;
          run_single_scenario<Tag, ValueType>(it, it2, target_value, func_t());
        }
      }
    }
  }
}

TEST(std_algorithms_non_mod_seq_ops, search_n) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace Search_n
}  // namespace stdalgos
}  // namespace Test
