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

#ifndef KOKKOS_SORT_BY_KEY_FREE_FUNCS_IMPL_HPP_
#define KOKKOS_SORT_BY_KEY_FREE_FUNCS_IMPL_HPP_

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA)

// Workaround for `Instruction 'shfl' without '.sync' is not supported on
// .target sm_70 and higher from PTX ISA version 6.4`.
// Also see https://github.com/NVIDIA/cub/pull/170.
#if !defined(CUB_USE_COOPERATIVE_GROUPS)
#define CUB_USE_COOPERATIVE_GROUPS
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"

#if defined(KOKKOS_COMPILER_CLANG)
// Some versions of Clang fail to compile Thrust, failing with errors like
// this:
//    <snip>/thrust/system/cuda/detail/core/agent_launcher.h:557:11:
//    error: use of undeclared identifier 'va_printf'
// The exact combination of versions for Clang and Thrust (or CUDA) for this
// failure was not investigated, however even very recent version combination
// (Clang 10.0.0 and Cuda 10.0) demonstrated failure.
//
// Defining _CubLog here locally allows us to avoid that code path, however
// disabling some debugging diagnostics
#pragma push_macro("_CubLog")
#ifdef _CubLog
#undef _CubLog
#endif
#define _CubLog
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#pragma pop_macro("_CubLog")
#else
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif

#pragma GCC diagnostic pop

#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif

#if defined(KOKKOS_ENABLE_ONEDPL) && \
    (ONEDPL_VERSION_MAJOR > 2022 ||  \
     (ONEDPL_VERSION_MAJOR == 2022 && ONEDPL_VERSION_MINOR >= 2))
#define KOKKOS_ONEDPL_HAS_SORT_BY_KEY
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#endif

namespace Kokkos::Impl {

template <typename T>
constexpr inline bool is_admissible_to_kokkos_sort_by_key =
    ::Kokkos::is_view<T>::value&& T::rank() == 1 &&
    (std::is_same<typename T::traits::array_layout,
                  Kokkos::LayoutLeft>::value ||
     std::is_same<typename T::traits::array_layout,
                  Kokkos::LayoutRight>::value ||
     std::is_same<typename T::traits::array_layout,
                  Kokkos::LayoutStride>::value);

template <class ViewType>
KOKKOS_INLINE_FUNCTION constexpr void
static_assert_is_admissible_to_kokkos_sort_by_key(const ViewType& /* view */) {
  static_assert(is_admissible_to_kokkos_sort_by_key<ViewType>,
                "Kokkos::sort_by_key only accepts 1D values View with "
                "LayoutRight, LayoutLeft or LayoutStride.");
}

// For the fallback implementation for sort_by_key using Kokkos::sort, we need
// to consider if Kokkos::sort defers to the fallback implementation that copies
// the array to the host and uses std::sort, see
// copy_to_host_run_stdsort_copy_back() in impl/Kokkos_SortImpl.hpp. If
// sort_on_device_v is true, we assume that std::sort doesn't copy data.
// Otherwise, we manually copy all data to the host and provide Kokkos::sort
// with a host execution space.
template <class ExecutionSpace, class Layout>
inline constexpr bool sort_on_device_v = false;

#if defined(KOKKOS_ENABLE_CUDA)
template <class Layout>
inline constexpr bool sort_on_device_v<Kokkos::Cuda, Layout> = true;

template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties, class... MaybeComparator>
void sort_by_key_cudathrust(
    const Kokkos::Cuda& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    MaybeComparator&&... maybeComparator) {
  const auto policy = thrust::cuda::par.on(exec.cuda_stream());
  auto keys_first   = ::Kokkos::Experimental::begin(keys);
  auto keys_last    = ::Kokkos::Experimental::end(keys);
  auto values_first = ::Kokkos::Experimental::begin(values);
  thrust::sort_by_key(policy, keys_first, keys_last, values_first,
                      std::forward<MaybeComparator>(maybeComparator)...);
}
#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
template <class Layout>
inline constexpr bool sort_on_device_v<Kokkos::HIP, Layout> = true;

template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties, class... MaybeComparator>
void sort_by_key_rocthrust(
    const Kokkos::HIP& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    MaybeComparator&&... maybeComparator) {
  const auto policy = thrust::hip::par.on(exec.hip_stream());
  auto keys_first   = ::Kokkos::Experimental::begin(keys);
  auto keys_last    = ::Kokkos::Experimental::end(keys);
  auto values_first = ::Kokkos::Experimental::begin(values);
  thrust::sort_by_key(policy, keys_first, keys_last, values_first,
                      std::forward<MaybeComparator>(maybeComparator)...);
}
#endif

#if defined(KOKKOS_ENABLE_ONEDPL)
template <class Layout>
inline constexpr bool sort_on_device_v<Kokkos::Experimental::SYCL, Layout> =
    std::is_same_v<Layout, Kokkos::LayoutLeft> ||
    std::is_same_v<Layout, Kokkos::LayoutRight>;

#ifdef KOKKOS_ONEDPL_HAS_SORT_BY_KEY
template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties, class... MaybeComparator>
void sort_by_key_onedpl(
    const Kokkos::Experimental::SYCL& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    MaybeComparator&&... maybeComparator) {
  if (keys.stride(0) != 1 && values.stride(0) != 1) {
    Kokkos::abort(
        "SYCL sort_by_key only supports rank-1 Views with stride(0) = 1.");
  }

  // Can't use Experimental::begin/end here since the oneDPL then assumes that
  // the data is on the host.
  auto queue  = exec.sycl_queue();
  auto policy = oneapi::dpl::execution::make_device_policy(queue);
  const int n = keys.extent(0);
  oneapi::dpl::sort_by_key(policy, keys.data(), keys.data() + n, values.data(),
                           std::forward<MaybeComparator>(maybeComparator)...);
}
#endif
#endif

template <typename ExecutionSpace, typename PermutationView, typename ViewType>
void applyPermutation(const ExecutionSpace& space,
                      const PermutationView& permutation,
                      const ViewType& view) {
  static_assert(std::is_integral<typename PermutationView::value_type>::value);

  auto view_copy = Kokkos::create_mirror(
      Kokkos::view_alloc(space, typename ExecutionSpace::memory_space{},
                         Kokkos::WithoutInitializing),
      view);
  Kokkos::deep_copy(space, view_copy, view);
  Kokkos::parallel_for(
      "Kokkos::sort_by_key_via_sort::permute_" + view.label(),
      Kokkos::RangePolicy<ExecutionSpace>(space, 0, view.extent(0)),
      KOKKOS_LAMBDA(int i) { view(i) = view_copy(permutation(i)); });
}

// FIXME_NVCC: nvcc has trouble compiling lambdas inside a function with
// variadic templates (sort_by_key_via_sort). Switch to using functors instead.
template <typename Permute>
struct IotaFunctor {
  Permute _permute;
  KOKKOS_FUNCTION void operator()(int i) const { _permute(i) = i; }
};
template <typename Keys>
struct LessFunctor {
  Keys _keys;
  KOKKOS_FUNCTION bool operator()(int i, int j) const {
    return _keys(i) < _keys(j);
  }
};

// FIXME_NVCC+MSVC: We can't use a lambda instead of a functor which gave us
// "For this host platform/dialect, an extended lambda cannot be defined inside
// the 'if' or 'else' block of a constexpr if statement"
template <typename Keys, typename Comparator>
struct KeyComparisonFunctor {
  Keys m_keys;
  Comparator m_comparator;
  KOKKOS_FUNCTION bool operator()(int i, int j) const {
    return m_comparator(m_keys(i), m_keys(j));
  }
};

template <class ExecutionSpace, class KeysDataType, class... KeysProperties,
          class ValuesDataType, class... ValuesProperties,
          class... MaybeComparator>
void sort_by_key_via_sort(
    const ExecutionSpace& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    MaybeComparator&&... maybeComparator) {
  static_assert(sizeof...(MaybeComparator) <= 1);

  auto const n = keys.size();

  Kokkos::View<unsigned int*, ExecutionSpace> permute(
      Kokkos::view_alloc(exec, Kokkos::WithoutInitializing,
                         "Kokkos::sort_by_key_via_sort::permute"),
      n);

  // iota
  Kokkos::parallel_for("Kokkos::sort_by_key_via_sort::iota",
                       Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                       IotaFunctor<decltype(permute)>{permute});

  using Layout =
      typename Kokkos::View<unsigned int*, ExecutionSpace>::array_layout;
  if constexpr (!sort_on_device_v<ExecutionSpace, Layout>) {
    auto host_keys = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::HostSpace{}, Kokkos::WithoutInitializing),
        keys);
    auto host_permute = Kokkos::create_mirror_view(
        Kokkos::view_alloc(Kokkos::HostSpace{}, Kokkos::WithoutInitializing),
        permute);
    Kokkos::deep_copy(exec, host_keys, keys);
    Kokkos::deep_copy(exec, host_permute, permute);

    exec.fence("Kokkos::Impl::sort_by_key_via_sort: before host sort");
    Kokkos::DefaultHostExecutionSpace host_exec;

    if constexpr (sizeof...(MaybeComparator) == 0) {
      Kokkos::sort(host_exec, host_permute,
                   LessFunctor<decltype(host_keys)>{host_keys});
    } else {
      auto keys_comparator =
          std::get<0>(std::tuple<MaybeComparator...>(maybeComparator...));
      Kokkos::sort(
          host_exec, host_permute,
          KeyComparisonFunctor<decltype(host_keys), decltype(keys_comparator)>{
              host_keys, keys_comparator});
    }
    host_exec.fence("Kokkos::Impl::sort_by_key_via_sort: after host sort");
    Kokkos::deep_copy(exec, permute, host_permute);
  } else {
#ifdef KOKKOS_ENABLE_SYCL
    auto* raw_keys_in_comparator = keys.data();
    auto stride                  = keys.stride(0);
    if constexpr (sizeof...(MaybeComparator) == 0) {
      Kokkos::sort(
          exec, permute, KOKKOS_LAMBDA(int i, int j) {
            return raw_keys_in_comparator[i * stride] <
                   raw_keys_in_comparator[j * stride];
          });
    } else {
      auto keys_comparator =
          std::get<0>(std::tuple<MaybeComparator...>(maybeComparator...));
      Kokkos::sort(
          exec, permute, KOKKOS_LAMBDA(int i, int j) {
            return keys_comparator(raw_keys_in_comparator[i * stride],
                                   raw_keys_in_comparator[j * stride]);
          });
    }
#else
    if constexpr (sizeof...(MaybeComparator) == 0) {
      Kokkos::sort(exec, permute, LessFunctor<decltype(keys)>{keys});
    } else {
      auto keys_comparator =
          std::get<0>(std::tuple<MaybeComparator...>(maybeComparator...));
      Kokkos::sort(
          exec, permute,
          KeyComparisonFunctor<decltype(keys), decltype(keys_comparator)>{
              keys, keys_comparator});
    }
#endif
  }

  applyPermutation(exec, permute, keys);
  applyPermutation(exec, permute, values);
}

// ------------------------------------------------------
//
// specialize cases for sorting by key without comparator
//
// ------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties>
void sort_by_key_device_view_without_comparator(
    const Kokkos::Cuda& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values) {
  sort_by_key_cudathrust(exec, keys, values);
}
#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties>
void sort_by_key_device_view_without_comparator(
    const Kokkos::HIP& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values) {
  sort_by_key_rocthrust(exec, keys, values);
}
#endif

#if defined(KOKKOS_ENABLE_ONEDPL)
template <class KeysDataType, class... KeysProperties, class ValuesDataType,
          class... ValuesProperties>
void sort_by_key_device_view_without_comparator(
    const Kokkos::Experimental::SYCL& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values) {
#ifdef KOKKOS_ONEDPL_HAS_SORT_BY_KEY
  if (keys.stride(0) == 1 && values.stride(0) == 1)
    sort_by_key_onedpl(exec, keys, values);
  else
#endif
    sort_by_key_via_sort(exec, keys, values);
}
#endif

// fallback case
template <class ExecutionSpace, class KeysDataType, class... KeysProperties,
          class ValuesDataType, class... ValuesProperties>
std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value>
sort_by_key_device_view_without_comparator(
    const ExecutionSpace& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values) {
  sort_by_key_via_sort(exec, keys, values);
}

// ---------------------------------------------------
//
// specialize cases for sorting by key with comparator
//
// ---------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
template <class ComparatorType, class KeysDataType, class... KeysProperties,
          class ValuesDataType, class... ValuesProperties>
void sort_by_key_device_view_with_comparator(
    const Kokkos::Cuda& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    const ComparatorType& comparator) {
  sort_by_key_cudathrust(exec, keys, values, comparator);
}
#endif

#if defined(KOKKOS_ENABLE_ROCTHRUST)
template <class ComparatorType, class KeysDataType, class... KeysProperties,
          class ValuesDataType, class... ValuesProperties>
void sort_by_key_device_view_with_comparator(
    const Kokkos::HIP& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    const ComparatorType& comparator) {
  sort_by_key_rocthrust(exec, keys, values, comparator);
}
#endif

#if defined(KOKKOS_ENABLE_ONEDPL)
template <class ComparatorType, class KeysDataType, class... KeysProperties,
          class ValuesDataType, class... ValuesProperties>
void sort_by_key_device_view_with_comparator(
    const Kokkos::Experimental::SYCL& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    const ComparatorType& comparator) {
#ifdef KOKKOS_ONEDPL_HAS_SORT_BY_KEY
  if (keys.stride(0) == 1 && values.stride(0) == 1)
    sort_by_key_onedpl(exec, keys, values, comparator);
  else
#endif
    sort_by_key_via_sort(exec, keys, values, comparator);
}
#endif

// fallback case
template <class ComparatorType, class ExecutionSpace, class KeysDataType,
          class... KeysProperties, class ValuesDataType,
          class... ValuesProperties>
std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value>
sort_by_key_device_view_with_comparator(
    const ExecutionSpace& exec,
    const Kokkos::View<KeysDataType, KeysProperties...>& keys,
    const Kokkos::View<ValuesDataType, ValuesProperties...>& values,
    const ComparatorType& comparator) {
  sort_by_key_via_sort(exec, keys, values, comparator);
}

#undef KOKKOS_ONEDPL_HAS_SORT_BY_KEY

}  // namespace Kokkos::Impl
#endif
