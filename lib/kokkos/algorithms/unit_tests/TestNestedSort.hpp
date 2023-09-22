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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_NESTED_SORT_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_NESTED_SORT_HPP

#include <gtest/gtest.h>
#include <unordered_set>
#include <random>
#include <Kokkos_Random.hpp>
#include <Kokkos_NestedSort.hpp>

namespace Test {
namespace NestedSortImpl {

// Comparator for sorting in descending order
template <typename Key>
struct GreaterThan {
  KOKKOS_FUNCTION constexpr bool operator()(const Key& lhs,
                                            const Key& rhs) const {
    return lhs > rhs;
  }
};

// Functor to test sort_team: each team responsible for sorting one array
template <typename ExecSpace, typename KeyViewType, typename OffsetViewType>
struct TeamSortFunctor {
  using TeamMem  = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using SizeType = typename KeyViewType::size_type;
  using KeyType  = typename KeyViewType::non_const_value_type;
  TeamSortFunctor(const KeyViewType& keys_, const OffsetViewType& offsets_,
                  bool sortDescending_)
      : keys(keys_), offsets(offsets_), sortDescending(sortDescending_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem& t) const {
    int i          = t.league_rank();
    SizeType begin = offsets(i);
    SizeType end   = offsets(i + 1);
    if (sortDescending)
      Kokkos::Experimental::sort_team(
          t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
          GreaterThan<KeyType>());
    else
      Kokkos::Experimental::sort_team(
          t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)));
  }
  KeyViewType keys;
  OffsetViewType offsets;
  bool sortDescending;
};

// Functor to test sort_by_key_team: each team responsible for sorting one array
template <typename ExecSpace, typename KeyViewType, typename ValueViewType,
          typename OffsetViewType>
struct TeamSortByKeyFunctor {
  using TeamMem  = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using SizeType = typename KeyViewType::size_type;
  using KeyType  = typename KeyViewType::non_const_value_type;
  TeamSortByKeyFunctor(const KeyViewType& keys_, const ValueViewType& values_,
                       const OffsetViewType& offsets_, bool sortDescending_)
      : keys(keys_),
        values(values_),
        offsets(offsets_),
        sortDescending(sortDescending_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem& t) const {
    int i          = t.league_rank();
    SizeType begin = offsets(i);
    SizeType end   = offsets(i + 1);
    if (sortDescending) {
      Kokkos::Experimental::sort_by_key_team(
          t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
          Kokkos::subview(values, Kokkos::make_pair(begin, end)),
          GreaterThan<KeyType>());
    } else {
      Kokkos::Experimental::sort_by_key_team(
          t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
          Kokkos::subview(values, Kokkos::make_pair(begin, end)));
    }
  }
  KeyViewType keys;
  ValueViewType values;
  OffsetViewType offsets;
  bool sortDescending;
};

// Functor to test sort_thread: each thread (multiple vector lanes) responsible
// for sorting one array
template <typename ExecSpace, typename KeyViewType, typename OffsetViewType>
struct ThreadSortFunctor {
  using TeamMem  = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using SizeType = typename KeyViewType::size_type;
  using KeyType  = typename KeyViewType::non_const_value_type;
  ThreadSortFunctor(const KeyViewType& keys_, const OffsetViewType& offsets_,
                    bool sortDescending_)
      : keys(keys_), offsets(offsets_), sortDescending(sortDescending_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem& t) const {
    int i = t.league_rank() * t.team_size() + t.team_rank();
    // Number of arrays to sort doesn't have to be divisible by team size, so
    // some threads may be idle.
    if (i < offsets.extent_int(0) - 1) {
      SizeType begin = offsets(i);
      SizeType end   = offsets(i + 1);
      if (sortDescending)
        Kokkos::Experimental::sort_thread(
            t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
            GreaterThan<KeyType>());
      else
        Kokkos::Experimental::sort_thread(
            t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)));
    }
  }
  KeyViewType keys;
  OffsetViewType offsets;
  bool sortDescending;
};

// Functor to test sort_by_key_thread
template <typename ExecSpace, typename KeyViewType, typename ValueViewType,
          typename OffsetViewType>
struct ThreadSortByKeyFunctor {
  using TeamMem  = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using SizeType = typename KeyViewType::size_type;
  using KeyType  = typename KeyViewType::non_const_value_type;
  ThreadSortByKeyFunctor(const KeyViewType& keys_, const ValueViewType& values_,
                         const OffsetViewType& offsets_, bool sortDescending_)
      : keys(keys_),
        values(values_),
        offsets(offsets_),
        sortDescending(sortDescending_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMem& t) const {
    int i = t.league_rank() * t.team_size() + t.team_rank();
    // Number of arrays to sort doesn't have to be divisible by team size, so
    // some threads may be idle.
    if (i < offsets.extent_int(0) - 1) {
      SizeType begin = offsets(i);
      SizeType end   = offsets(i + 1);
      if (sortDescending) {
        Kokkos::Experimental::sort_by_key_thread(
            t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
            Kokkos::subview(values, Kokkos::make_pair(begin, end)),
            GreaterThan<KeyType>());
      } else {
        Kokkos::Experimental::sort_by_key_thread(
            t, Kokkos::subview(keys, Kokkos::make_pair(begin, end)),
            Kokkos::subview(values, Kokkos::make_pair(begin, end)));
      }
    }
  }
  KeyViewType keys;
  ValueViewType values;
  OffsetViewType offsets;
  bool sortDescending;
};

// Generate the offsets view for a set of n packed arrays, each with uniform
// random length in [0,k]. Array i will occupy the indices [offsets(i),
// offsets(i+1)), like a row in a CRS graph. Returns the total length of all the
// arrays.
template <typename OffsetViewType>
size_t randomPackedArrayOffsets(unsigned n, unsigned k,
                                OffsetViewType& offsets) {
  offsets          = OffsetViewType("Offsets", n + 1);
  auto offsetsHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), offsets);
  std::mt19937 gen;
  std::uniform_int_distribution<> distrib(0, k);
  // This will leave offsetsHost(n) == 0.
  std::generate(offsetsHost.data(), offsetsHost.data() + n,
                [&]() { return distrib(gen); });
  // Exclusive prefix-sum to get offsets
  size_t accum = 0;
  for (unsigned i = 0; i <= n; i++) {
    size_t num     = offsetsHost(i);
    offsetsHost(i) = accum;
    accum += num;
  }
  Kokkos::deep_copy(offsets, offsetsHost);
  return offsetsHost(n);
}

template <typename ValueViewType>
ValueViewType uniformRandomViewFill(size_t totalLength,
                                    typename ValueViewType::value_type minVal,
                                    typename ValueViewType::value_type maxVal) {
  ValueViewType vals("vals", totalLength);
  Kokkos::Random_XorShift64_Pool<typename ValueViewType::execution_space> g(
      1931);
  Kokkos::fill_random(vals, g, minVal, maxVal);
  return vals;
}

template <class ExecutionSpace, typename KeyType>
void test_nested_sort_impl(unsigned narray, unsigned n, bool useTeams,
                           bool customCompare, KeyType minKey, KeyType maxKey) {
  using KeyViewType    = Kokkos::View<KeyType*, ExecutionSpace>;
  using OffsetViewType = Kokkos::View<unsigned*, ExecutionSpace>;
  using TeamPol        = Kokkos::TeamPolicy<ExecutionSpace>;
  OffsetViewType offsets;
  size_t totalLength = randomPackedArrayOffsets(narray, n, offsets);
  KeyViewType keys =
      uniformRandomViewFill<KeyViewType>(totalLength, minKey, maxKey);
  // note: doing create_mirror because we always want this to be a separate
  // copy, even if keys is already host-accessible. keysHost becomes the correct
  // result to compare against.
  auto keysHost = Kokkos::create_mirror(Kokkos::HostSpace(), keys);
  Kokkos::deep_copy(keysHost, keys);
  auto offsetsHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  // Sort the same arrays on host to compare against
  for (unsigned i = 0; i < narray; i++) {
    KeyType* begin = keysHost.data() + offsetsHost(i);
    KeyType* end   = keysHost.data() + offsetsHost(i + 1);
    if (customCompare)
      std::sort(begin, end,
                [](const KeyType& a, const KeyType& b) { return a > b; });
    else
      std::sort(begin, end);
  }
  if (useTeams) {
    int vectorLen = std::min<int>(4, TeamPol::vector_length_max());
    TeamPol policy(narray, Kokkos::AUTO(), vectorLen);
    Kokkos::parallel_for(
        policy, TeamSortFunctor<ExecutionSpace, KeyViewType, OffsetViewType>(
                    keys, offsets, customCompare));
  } else {
    ThreadSortFunctor<ExecutionSpace, KeyViewType, OffsetViewType> functor(
        keys, offsets, customCompare);
    int vectorLen = std::min<int>(4, TeamPol::vector_length_max());
    TeamPol dummy(1, Kokkos::AUTO(), vectorLen);
    int teamSize =
        dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
    int numTeams = (narray + teamSize - 1) / teamSize;
    Kokkos::parallel_for(TeamPol(numTeams, teamSize, vectorLen), functor);
  }
  auto keysOut = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keys);
  std::string testLabel = useTeams ? "sort_team" : "sort_thread";
  for (unsigned i = 0; i < keys.extent(0); i++) {
    EXPECT_EQ(keysOut(i), keysHost(i))
        << testLabel << ": after sorting, key at index " << i
        << " is incorrect.";
  }
}

template <class ExecutionSpace, typename KeyType, typename ValueType>
void test_nested_sort_by_key_impl(unsigned narray, unsigned n, bool useTeams,
                                  bool customCompare, KeyType minKey,
                                  KeyType maxKey, ValueType minVal,
                                  ValueType maxVal) {
  using KeyViewType    = Kokkos::View<KeyType*, ExecutionSpace>;
  using ValueViewType  = Kokkos::View<ValueType*, ExecutionSpace>;
  using OffsetViewType = Kokkos::View<unsigned*, ExecutionSpace>;
  using TeamPol        = Kokkos::TeamPolicy<ExecutionSpace>;
  OffsetViewType offsets;
  size_t totalLength = randomPackedArrayOffsets(narray, n, offsets);
  KeyViewType keys =
      uniformRandomViewFill<KeyViewType>(totalLength, minKey, maxKey);
  ValueViewType values =
      uniformRandomViewFill<ValueViewType>(totalLength, minVal, maxVal);
  // note: doing create_mirror because we always want this to be a separate
  // copy, even if keys/vals are already host-accessible. keysHost and valsHost
  // becomes the correct result to compare against.
  auto keysHost   = Kokkos::create_mirror(Kokkos::HostSpace(), keys);
  auto valuesHost = Kokkos::create_mirror(Kokkos::HostSpace(), values);
  Kokkos::deep_copy(keysHost, keys);
  Kokkos::deep_copy(valuesHost, values);
  auto offsetsHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsets);
  // Sort the same arrays on host to compare against
  for (unsigned i = 0; i < narray; i++) {
    // std:: doesn't have a sort_by_key, so sort a vector of key-value pairs
    // instead
    using KV = std::pair<KeyType, ValueType>;
    std::vector<KV> keysAndValues(offsetsHost(i + 1) - offsetsHost(i));
    for (unsigned j = 0; j < keysAndValues.size(); j++) {
      keysAndValues[j].first  = keysHost(offsetsHost(i) + j);
      keysAndValues[j].second = valuesHost(offsetsHost(i) + j);
    }
    if (customCompare) {
      std::sort(keysAndValues.begin(), keysAndValues.end(),
                [](const KV& a, const KV& b) { return a.first > b.first; });
    } else {
      std::sort(keysAndValues.begin(), keysAndValues.end(),
                [](const KV& a, const KV& b) { return a.first < b.first; });
    }
    // Copy back from pairs to views
    for (unsigned j = 0; j < keysAndValues.size(); j++) {
      keysHost(offsetsHost(i) + j)   = keysAndValues[j].first;
      valuesHost(offsetsHost(i) + j) = keysAndValues[j].second;
    }
  }
  if (useTeams) {
    int vectorLen = std::min<int>(4, TeamPol::vector_length_max());
    TeamPol policy(narray, Kokkos::AUTO(), vectorLen);
    Kokkos::parallel_for(
        policy, TeamSortByKeyFunctor<ExecutionSpace, KeyViewType, ValueViewType,
                                     OffsetViewType>(keys, values, offsets,
                                                     customCompare));
  } else {
    ThreadSortByKeyFunctor<ExecutionSpace, KeyViewType, ValueViewType,
                           OffsetViewType>
        functor(keys, values, offsets, customCompare);
    int vectorLen = std::min<int>(4, TeamPol::vector_length_max());
    TeamPol dummy(1, Kokkos::AUTO(), vectorLen);
    int teamSize =
        dummy.team_size_recommended(functor, Kokkos::ParallelForTag());
    int numTeams = (narray + teamSize - 1) / teamSize;
    Kokkos::parallel_for(TeamPol(numTeams, teamSize, vectorLen), functor);
  }
  auto keysOut = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), keys);
  auto valuesOut =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
  std::string testLabel = useTeams ? "sort_by_key_team" : "sort_by_key_thread";
  // First, compare keys since they will always match exactly
  for (unsigned i = 0; i < keys.extent(0); i++) {
    EXPECT_EQ(keysOut(i), keysHost(i))
        << testLabel << ": after sorting, key at index " << i
        << " is incorrect.";
  }
  // Kokkos::sort_by_key_X is not stable, so if a key happens to
  // appear more than once, the order of the values may not match exactly.
  // But the set of values for a given key should be identical.
  unsigned keyStart = 0;
  while (keyStart < keys.extent(0)) {
    KeyType key     = keysHost(keyStart);
    unsigned keyEnd = keyStart + 1;
    while (keyEnd < keys.extent(0) && keysHost(keyEnd) == key) keyEnd++;
    std::unordered_multiset<ValueType> correctVals;
    std::unordered_multiset<ValueType> outputVals;
    for (unsigned i = keyStart; i < keyEnd; i++) {
      correctVals.insert(valuesHost(i));
      outputVals.insert(valuesOut(i));
    }
    // Check one value at a time that they match
    for (auto it = correctVals.begin(); it != correctVals.end(); it++) {
      ValueType val = *it;
      EXPECT_TRUE(outputVals.find(val) != outputVals.end())
          << testLabel << ": after sorting, value " << val
          << " corresponding to key " << key << " is missing.";
      EXPECT_EQ(correctVals.count(val), outputVals.count(val))
          << testLabel << ": after sorting, the number of occurences of value "
          << val << " corresponding to key " << key << " changed.";
    }
    keyStart = keyEnd;
  }
}

template <class ExecutionSpace, typename KeyType>
void test_nested_sort(unsigned int N, KeyType minKey, KeyType maxKey) {
  // 2nd arg: true = team-level, false = thread-level.
  // 3rd arg: true = custom comparator, false = default comparator.
  test_nested_sort_impl<ExecutionSpace, KeyType>(N, N, true, false, minKey,
                                                 maxKey);
  test_nested_sort_impl<ExecutionSpace, KeyType>(N, N, true, true, minKey,
                                                 maxKey);
  test_nested_sort_impl<ExecutionSpace, KeyType>(N, N, false, false, minKey,
                                                 maxKey);
  test_nested_sort_impl<ExecutionSpace, KeyType>(N, N, false, true, minKey,
                                                 maxKey);
}

template <class ExecutionSpace, typename KeyType, typename ValueType>
void test_nested_sort_by_key(unsigned int N, KeyType minKey, KeyType maxKey,
                             ValueType minVal, ValueType maxVal) {
  // 2nd arg: true = team-level, false = thread-level.
  // 3rd arg: true = custom comparator, false = default comparator.
  test_nested_sort_by_key_impl<ExecutionSpace, KeyType, ValueType>(
      N, N, true, false, minKey, maxKey, minVal, maxVal);
  test_nested_sort_by_key_impl<ExecutionSpace, KeyType, ValueType>(
      N, N, true, true, minKey, maxKey, minVal, maxVal);
  test_nested_sort_by_key_impl<ExecutionSpace, KeyType, ValueType>(
      N, N, false, false, minKey, maxKey, minVal, maxVal);
  test_nested_sort_by_key_impl<ExecutionSpace, KeyType, ValueType>(
      N, N, false, true, minKey, maxKey, minVal, maxVal);
}
}  // namespace NestedSortImpl

TEST(TEST_CATEGORY, NestedSort) {
  using ExecutionSpace = TEST_EXECSPACE;
  NestedSortImpl::test_nested_sort<ExecutionSpace, unsigned>(171, 0U, UINT_MAX);
  NestedSortImpl::test_nested_sort<ExecutionSpace, float>(42, -1e6f, 1e6f);
  NestedSortImpl::test_nested_sort<ExecutionSpace, char>(67, CHAR_MIN,
                                                         CHAR_MAX);
}

TEST(TEST_CATEGORY, NestedSortByKey) {
  using ExecutionSpace = TEST_EXECSPACE;

  // Second/third template arguments are key and value respectively.
  // In sort_by_key_X functions, a key view and a value view are both permuted
  // to make the keys sorted. This means that the value type doesn't need to be
  // ordered, unlike key
  NestedSortImpl::test_nested_sort_by_key<ExecutionSpace, unsigned, unsigned>(
      161, 0U, UINT_MAX, 0U, UINT_MAX);
  NestedSortImpl::test_nested_sort_by_key<ExecutionSpace, float, char>(
      267, -1e6f, 1e6f, CHAR_MIN, CHAR_MAX);
  NestedSortImpl::test_nested_sort_by_key<ExecutionSpace, char, double>(
      11, CHAR_MIN, CHAR_MAX, 2.718, 3.14);
}

}  // namespace Test
#endif
