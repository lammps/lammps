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

#ifndef KOKKOS_KOKKOS_TUNERS_HPP
#define KOKKOS_KOKKOS_TUNERS_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

#include <array>
#include <utility>
#include <tuple>
#include <string>
#include <vector>
#include <map>
#include <cassert>

namespace Kokkos {
namespace Tools {

namespace Experimental {

// forward declarations
SetOrRange make_candidate_set(size_t size, int64_t* data);
bool have_tuning_tool();
size_t declare_output_type(const std::string&,
                           Kokkos::Tools::Experimental::VariableInfo);
void request_output_values(size_t, size_t,
                           Kokkos::Tools::Experimental::VariableValue*);
VariableValue make_variable_value(size_t, int64_t);
VariableValue make_variable_value(size_t, double);
SetOrRange make_candidate_range(double lower, double upper, double step,
                                bool openLower, bool openUpper);
size_t get_new_context_id();
void begin_context(size_t context_id);
void end_context(size_t context_id);
namespace Impl {

/** We're going to take in search space descriptions
 * as nested maps, which aren't efficient to
 * iterate across by index. These are very similar
 * to nested maps, but better for index-based lookup
 */
template <typename ValueType, typename ContainedType>
struct ValueHierarchyNode;

template <typename ValueType, typename ContainedType>
struct ValueHierarchyNode {
  std::vector<ValueType> root_values;
  std::vector<ContainedType> sub_values;
  void add_root_value(const ValueType& in) noexcept {
    root_values.push_back(in);
  }
  void add_sub_container(const ContainedType& in) { sub_values.push_back(in); }
  const ValueType& get_root_value(const size_t index) const {
    return root_values[index];
  }
  const ContainedType& get_sub_value(const size_t index) const {
    return sub_values[index];
  }
};

template <typename ValueType>
struct ValueHierarchyNode<ValueType, void> {
  std::vector<ValueType> root_values;
  explicit ValueHierarchyNode(std::vector<ValueType> rv)
      : root_values(std::move(rv)) {}
  void add_root_value(const ValueType& in) noexcept {
    root_values.push_back(in);
  }
  const ValueType& get_root_value(const size_t index) const {
    return root_values[index];
  }
};

/** For a given nested map type, we need a way to
 * declare the equivalent ValueHierarchyNode
 * structure
 */

template <class NestedMap>
struct MapTypeConverter;

// Vectors are our lowest-level, no nested values
template <class T>
struct MapTypeConverter<std::vector<T>> {
  using type = ValueHierarchyNode<T, void>;
};

// Maps contain both the "root" types and sub-vectors
template <class K, class V>
struct MapTypeConverter<std::map<K, V>> {
  using type = ValueHierarchyNode<K, typename MapTypeConverter<V>::type>;
};

/**
 * We also need to be able to construct a ValueHierarchyNode set from a
 * map
 */

template <class NestedMap>
struct ValueHierarchyConstructor;

// Vectors are our lowest-level, no nested values. Just fill in the fundamental
// values
template <class T>
struct ValueHierarchyConstructor<std::vector<T>> {
  using return_type = typename MapTypeConverter<std::vector<T>>::type;
  static return_type build(const std::vector<T>& in) { return return_type{in}; }
};

// For maps, we need to fill in the fundamental values, and construct child
// nodes
template <class K, class V>
struct ValueHierarchyConstructor<std::map<K, V>> {
  using return_type = typename MapTypeConverter<std::map<K, V>>::type;
  static return_type build(const std::map<K, V>& in) {
    return_type node_to_build;
    for (auto& entry : in) {
      node_to_build.add_root_value(entry.first);
      node_to_build.add_sub_container(
          ValueHierarchyConstructor<V>::build(entry.second));
    }
    return node_to_build;
  }
};

/**
 * We're going to be declaring a sparse multidimensional
 * tuning space as a set of nested maps. The innermost level
 * will be a vector. The dimensionality of such a space is the number of
 * maps + 1.
 *
 * The following templates implement such logic recursively
 */
template <class InspectForDepth>
struct get_space_dimensionality;

// The dimensionality of a vector is 1
template <class T>
struct get_space_dimensionality<std::vector<T>> {
  static constexpr int value = 1;
};

// The dimensionality of a map is 1 (the map) plus the dimensionality
// of the map's value type
template <class K, class V>
struct get_space_dimensionality<std::map<K, V>> {
  static constexpr int value = 1 + get_space_dimensionality<V>::value;
};

template <class T, int N>
struct n_dimensional_sparse_structure;

template <class T>
struct n_dimensional_sparse_structure<T, 1> {
  using type = std::vector<T>;
};

template <class T, int N>
struct n_dimensional_sparse_structure {
  using type =
      std::map<T, typename n_dimensional_sparse_structure<T, N - 1>::type>;
};

/**
 * This is the ugly part of this implementation: mapping a set of doubles in
 * [0.0,1.0) into a point in this multidimensional space. We're going to
 * implement this concept recursively, building up a tuple at each level.
 */

// First, a helper to get the value in one dimension
template <class Container>
struct DimensionValueExtractor;

// At any given level, just return your value at that level
template <class RootType, class Subtype>
struct DimensionValueExtractor<ValueHierarchyNode<RootType, Subtype>> {
  static RootType get(const ValueHierarchyNode<RootType, Subtype>& dimension,
                      double fraction_to_traverse) {
    size_t index = dimension.root_values.size() * fraction_to_traverse;
    return dimension.get_root_value(index);
  }
};

/** Now we're going to do the full "get a point in the space".
 * At a root level, we'll take in a ValueHierarchyNode and a set of doubles
 * representing the value in [0.0,1.0) we want to pick
 */

// At the bottom level, we have one double and a base-level ValueHierarchyNode

template <class HierarchyNode, class... InterpolationIndices>
struct GetMultidimensionalPoint;

template <class ValueType>
struct GetMultidimensionalPoint<ValueHierarchyNode<ValueType, void>, double> {
  using node_type   = ValueHierarchyNode<ValueType, void>;
  using return_type = std::tuple<ValueType>;
  static return_type build(const node_type& in, double index) {
    return std::make_tuple(DimensionValueExtractor<node_type>::get(in, index));
  }
};

// At levels above the bottom, we tuple_cat the result of our child on the end
// of our own tuple
template <class ValueType, class Subtype, class... Indices>
struct GetMultidimensionalPoint<ValueHierarchyNode<ValueType, Subtype>, double,
                                Indices...> {
  using node_type = ValueHierarchyNode<ValueType, Subtype>;
  using sub_tuple =
      typename GetMultidimensionalPoint<Subtype, Indices...>::return_type;
  using return_type = decltype(std::tuple_cat(
      std::declval<std::tuple<ValueType>>(), std::declval<sub_tuple>()));
  static return_type build(const node_type& in, double fraction_to_traverse,
                           Indices... indices) {
    size_t index         = in.sub_values.size() * fraction_to_traverse;
    auto dimension_value = std::make_tuple(
        DimensionValueExtractor<node_type>::get(in, fraction_to_traverse));
    return std::tuple_cat(dimension_value,
                          GetMultidimensionalPoint<Subtype, Indices...>::build(
                              in.get_sub_value(index), indices...));
  }
};

template <typename PointType, class ArrayType, size_t... Is>
auto get_point_helper(const PointType& in, const ArrayType& indices,
                      std::index_sequence<Is...>) {
  using helper = GetMultidimensionalPoint<
      PointType,
      decltype(std::get<Is>(std::declval<ArrayType>()).value.double_value)...>;
  return helper::build(in, std::get<Is>(indices).value.double_value...);
}

template <typename PointType, typename ArrayType>
struct GetPoint;

template <typename PointType, size_t X>
struct GetPoint<PointType,
                std::array<Kokkos::Tools::Experimental::VariableValue, X>> {
  using index_set_type =
      std::array<Kokkos::Tools::Experimental::VariableValue, X>;
  static auto build(const PointType& in, const index_set_type& indices) {
    return get_point_helper(in, indices, std::make_index_sequence<X>{});
  }
};

template <typename PointType, typename ArrayType>
auto get_point(const PointType& point, const ArrayType& indices) {
  return GetPoint<PointType, ArrayType>::build(point, indices);
}

}  // namespace Impl

template <template <class...> class Container, size_t MaxDimensionSize = 100,
          class... TemplateArguments>
class MultidimensionalSparseTuningProblem {
 public:
  using ProblemSpaceInput = Container<TemplateArguments...>;
  static constexpr int space_dimensionality =
      Impl::get_space_dimensionality<ProblemSpaceInput>::value;
  static constexpr size_t max_space_dimension_size = MaxDimensionSize;
  static constexpr double tuning_min               = 0.0;
  static constexpr double tuning_max               = 0.999;
  static constexpr double tuning_step = tuning_max / max_space_dimension_size;

  using StoredProblemSpace =
      typename Impl::MapTypeConverter<ProblemSpaceInput>::type;
  using HierarchyConstructor =
      typename Impl::ValueHierarchyConstructor<Container<TemplateArguments...>>;

  using ValueArray = std::array<Kokkos::Tools::Experimental::VariableValue,
                                space_dimensionality>;

 private:
  StoredProblemSpace m_space;
  std::array<size_t, space_dimensionality> variable_ids;
  size_t context;

 public:
  MultidimensionalSparseTuningProblem() = default;
  MultidimensionalSparseTuningProblem(ProblemSpaceInput space,
                                      const std::vector<std::string>& names)
      : m_space(HierarchyConstructor::build(space)) {
    assert(names.size() == space_dimensionality);
    for (unsigned long x = 0; x < names.size(); ++x) {
      VariableInfo info;
      info.type = Kokkos::Tools::Experimental::ValueType::kokkos_value_double;
      info.category = Kokkos::Tools::Experimental::StatisticalCategory::
          kokkos_value_interval;
      info.valueQuantity =
          Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_range;
      info.candidates = Kokkos::Tools::Experimental::make_candidate_range(
          tuning_min, tuning_max, tuning_step, true, true);
      variable_ids[x] = declare_output_type(names[x], info);
    }
  }

  auto begin() {
    context = Kokkos::Tools::Experimental::get_new_context_id();
    ValueArray values;
    for (int x = 0; x < space_dimensionality; ++x) {
      values[x] = Kokkos::Tools::Experimental::make_variable_value(
          variable_ids[x], 0.0);
    }
    begin_context(context);
    request_output_values(context, space_dimensionality, values.data());
    return get_point(m_space, values);
  }

  auto end() { end_context(context); }
};

template <size_t MaxDimensionSize = 100, template <class...> class Container,
          class... TemplateArguments>
auto make_multidimensional_sparse_tuning_problem(
    const Container<TemplateArguments...>& in, std::vector<std::string> names) {
  return MultidimensionalSparseTuningProblem<Container, MaxDimensionSize,
                                             TemplateArguments...>(in, names);
}
class TeamSizeTuner {
 private:
  using SpaceDescription = std::map<int64_t, std::vector<int64_t>>;
  using TunerType = decltype(make_multidimensional_sparse_tuning_problem<20>(
      std::declval<SpaceDescription>(),
      std::declval<std::vector<std::string>>()));
  TunerType tuner;

 public:
  TeamSizeTuner()        = default;
  TeamSizeTuner& operator=(const TeamSizeTuner& other) = default;
  TeamSizeTuner(const TeamSizeTuner& other)            = default;
  TeamSizeTuner& operator=(TeamSizeTuner&& other) = default;
  TeamSizeTuner(TeamSizeTuner&& other)            = default;
  template <typename ViableConfigurationCalculator, typename Functor,
            typename TagType, typename... Properties>
  TeamSizeTuner(const std::string& name,
                Kokkos::TeamPolicy<Properties...>& policy,
                const Functor& functor, const TagType& tag,
                ViableConfigurationCalculator calc) {
    using PolicyType           = Kokkos::TeamPolicy<Properties...>;
    auto initial_vector_length = policy.impl_vector_length();
    if (initial_vector_length < 1) {
      policy.impl_set_vector_length(1);
    }
    /**
     * Here we attempt to enumerate all of the possible configurations
     * to expose to an autotuner. There are three possibilities
     *
     * 1) We're tuning both vector length and team size
     * 2) We're tuning vector length but not team size
     * 3) We're tuning team size but not vector length
     *
     * (In the fourth case where nothing is tuned
     * this function won't be called)
     *
     * The set of valid team sizes is dependent on
     * a vector length, so this leads to three
     * algorithms
     *
     * 1) Loop over vector lengths to get the set
     *    of team sizes for each vector length,
     *    add it all to the set
     * 2) Loop over vector lengths to see if the
     *    provided team size is valid for that
     *    vector length. If so, add it
     * 3) A special case of (1) in which we only
     *    have one vector length
     *
     */
    SpaceDescription space_description;

    auto max_vector_length = PolicyType::vector_length_max();
    std::vector<int64_t> allowed_vector_lengths;

    if (policy.impl_auto_vector_length()) {  // case 1 or 2
      for (int vector_length = max_vector_length; vector_length >= 1;
           vector_length /= 2) {
        policy.impl_set_vector_length(vector_length);
        /**
         * Figuring out whether a vector length is valid depends
         * on whether we're in case 1 (tune everything) or 2 (just tune vector
         * length)
         *
         * If we're tuning everything, all legal vector lengths are valid.
         * If we're just tuning vector length, we need to check that if we
         * set this vector length, the team size provided will be valid.
         *
         * These are the left and right hand sides of the "or" in this
         * conditional, respectively.
         */
        auto max_team_size = calc.get_max_team_size(policy, functor, tag);
        if ((policy.impl_auto_team_size()) ||
            (policy.team_size() <= max_team_size)) {
          allowed_vector_lengths.push_back(vector_length);
        }
      }
    } else {  // case 3, there's only one vector length to care about
      allowed_vector_lengths.push_back(policy.impl_vector_length());
    }

    for (const auto vector_length : allowed_vector_lengths) {
      std::vector<int64_t> allowed_team_sizes;
      policy.impl_set_vector_length(vector_length);
      auto max_team_size = calc.get_max_team_size(policy, functor, tag);
      if (policy.impl_auto_team_size()) {  // case 1 or 3, try all legal team
                                           // sizes
        for (int team_size = max_team_size; team_size >= 1; team_size /= 2) {
          allowed_team_sizes.push_back(team_size);
        }
      } else {  // case 2, just try the provided team size
        allowed_team_sizes.push_back(policy.team_size());
      }
      space_description[vector_length] = allowed_team_sizes;
    }
    tuner = make_multidimensional_sparse_tuning_problem<20>(
        space_description, {std::string(name + "_vector_length"),
                            std::string(name + "_team_size")});
    policy.impl_set_vector_length(initial_vector_length);
  }

  template <typename... Properties>
  void tune(Kokkos::TeamPolicy<Properties...>& policy) {
    if (Kokkos::Tools::Experimental::have_tuning_tool()) {
      auto configuration = tuner.begin();
      auto team_size     = std::get<1>(configuration);
      auto vector_length = std::get<0>(configuration);
      if (vector_length > 0) {
        policy.impl_set_team_size(team_size);
        policy.impl_set_vector_length(vector_length);
      }
    }
  }
  void end() {
    if (Kokkos::Tools::Experimental::have_tuning_tool()) {
      tuner.end();
    }
  }

 private:
};

namespace Impl {

template <typename T>
void fill_tile(std::vector<T>& cont, int tile_size) {
  for (int x = 1; x < tile_size; x *= 2) {
    cont.push_back(x);
  }
}
template <typename T, typename Mapped>
void fill_tile(std::map<T, Mapped>& cont, int tile_size) {
  for (int x = 1; x < tile_size; x *= 2) {
    fill_tile(cont[x], tile_size / x);
  }
}
}  // namespace Impl

template <int MDRangeRank>
struct MDRangeTuner {
 private:
  static constexpr int rank       = MDRangeRank;
  static constexpr int max_slices = 15;
  using SpaceDescription =
      typename Impl::n_dimensional_sparse_structure<int, rank>::type;
  using TunerType =
      decltype(make_multidimensional_sparse_tuning_problem<max_slices>(
          std::declval<SpaceDescription>(),
          std::declval<std::vector<std::string>>()));
  TunerType tuner;

 public:
  MDRangeTuner() = default;
  template <typename Functor, typename TagType, typename Calculator,
            typename... Properties>
  MDRangeTuner(const std::string& name,
               const Kokkos::MDRangePolicy<Properties...>& policy,
               const Functor& functor, const TagType& tag, Calculator calc) {
    SpaceDescription desc;
    int max_tile_size =
        calc.get_mdrange_max_tile_size_product(policy, functor, tag);
    Impl::fill_tile(desc, max_tile_size);
    std::vector<std::string> feature_names;
    for (int x = 0; x < rank; ++x) {
      feature_names.push_back(name + "_tile_size_" + std::to_string(x));
    }
    tuner = make_multidimensional_sparse_tuning_problem<max_slices>(
        desc, feature_names);
  }
  template <typename Policy, typename Tuple, size_t... Indices>
  void set_policy_tile(Policy& policy, const Tuple& tuple,
                       const std::index_sequence<Indices...>&) {
    policy.impl_change_tile_size({std::get<Indices>(tuple)...});
  }
  template <typename... Properties>
  void tune(Kokkos::MDRangePolicy<Properties...>& policy) {
    if (Kokkos::Tools::Experimental::have_tuning_tool()) {
      auto configuration = tuner.begin();
      set_policy_tile(policy, configuration, std::make_index_sequence<rank>{});
    }
  }
  void end() {
    if (Kokkos::Tools::Experimental::have_tuning_tool()) {
      tuner.end();
    }
  }
};

}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

#endif
