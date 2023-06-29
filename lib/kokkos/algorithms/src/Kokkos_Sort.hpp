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

#ifndef KOKKOS_SORT_HPP_
#define KOKKOS_SORT_HPP_
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_SORT
#endif

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <algorithm>

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

#if defined(KOKKOS_ENABLE_ONEDPL)
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#endif

namespace Kokkos {

namespace Impl {

template <class DstViewType, class SrcViewType, int Rank = DstViewType::rank>
struct CopyOp;

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 1> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    dst(i_dst) = src(i_src);
  }
};

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 2> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    for (int j = 0; j < (int)dst.extent(1); j++) dst(i_dst, j) = src(i_src, j);
  }
};

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 3> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    for (int j = 0; j < dst.extent(1); j++)
      for (int k = 0; k < dst.extent(2); k++)
        dst(i_dst, j, k) = src(i_src, j, k);
  }
};
}  // namespace Impl

//----------------------------------------------------------------------------

template <class KeyViewType, class BinSortOp,
          class Space    = typename KeyViewType::device_type,
          class SizeType = typename KeyViewType::memory_space::size_type>
class BinSort {
 public:
  template <class DstViewType, class SrcViewType>
  struct copy_functor {
    using src_view_type = typename SrcViewType::const_type;

    using copy_op = Impl::CopyOp<DstViewType, src_view_type>;

    DstViewType dst_values;
    src_view_type src_values;
    int dst_offset;

    copy_functor(DstViewType const& dst_values_, int const& dst_offset_,
                 SrcViewType const& src_values_)
        : dst_values(dst_values_),
          src_values(src_values_),
          dst_offset(dst_offset_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i) const {
      copy_op::copy(dst_values, i + dst_offset, src_values, i);
    }
  };

  template <class DstViewType, class PermuteViewType, class SrcViewType>
  struct copy_permute_functor {
    // If a Kokkos::View then can generate constant random access
    // otherwise can only use the constant type.

    using src_view_type = std::conditional_t<
        Kokkos::is_view<SrcViewType>::value,
        Kokkos::View<typename SrcViewType::const_data_type,
                     typename SrcViewType::array_layout,
                     typename SrcViewType::device_type
#if !defined(KOKKOS_COMPILER_NVHPC)  // FIXME_NVHPC
                     ,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess>
#endif
                     >,
        typename SrcViewType::const_type>;

    using perm_view_type = typename PermuteViewType::const_type;

    using copy_op = Impl::CopyOp<DstViewType, src_view_type>;

    DstViewType dst_values;
    perm_view_type sort_order;
    src_view_type src_values;
    int src_offset;

    copy_permute_functor(DstViewType const& dst_values_,
                         PermuteViewType const& sort_order_,
                         SrcViewType const& src_values_, int const& src_offset_)
        : dst_values(dst_values_),
          sort_order(sort_order_),
          src_values(src_values_),
          src_offset(src_offset_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i) const {
      copy_op::copy(dst_values, i, src_values, src_offset + sort_order(i));
    }
  };

  // Naming this alias "execution_space" would be problematic since it would be
  // considered as execution space for the various functors which might use
  // another execution space through sort() or create_permute_vector().
  using exec_space  = typename Space::execution_space;
  using bin_op_type = BinSortOp;

  struct bin_count_tag {};
  struct bin_offset_tag {};
  struct bin_binning_tag {};
  struct bin_sort_bins_tag {};

 public:
  using size_type  = SizeType;
  using value_type = size_type;

  using offset_type    = Kokkos::View<size_type*, Space>;
  using bin_count_type = Kokkos::View<const int*, Space>;

  using const_key_view_type = typename KeyViewType::const_type;

  // If a Kokkos::View then can generate constant random access
  // otherwise can only use the constant type.

  using const_rnd_key_view_type = std::conditional_t<
      Kokkos::is_view<KeyViewType>::value,
      Kokkos::View<typename KeyViewType::const_data_type,
                   typename KeyViewType::array_layout,
                   typename KeyViewType::device_type,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess> >,
      const_key_view_type>;

  using non_const_key_scalar = typename KeyViewType::non_const_value_type;
  using const_key_scalar     = typename KeyViewType::const_value_type;

  using bin_count_atomic_type =
      Kokkos::View<int*, Space, Kokkos::MemoryTraits<Kokkos::Atomic> >;

 private:
  const_key_view_type keys;
  const_rnd_key_view_type keys_rnd;

 public:
  BinSortOp bin_op;
  offset_type bin_offsets;
  bin_count_atomic_type bin_count_atomic;
  bin_count_type bin_count_const;
  offset_type sort_order;

  int range_begin;
  int range_end;
  bool sort_within_bins;

 public:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED BinSort() = default;
#else
  BinSort() = delete;
#endif

  //----------------------------------------
  // Constructor: takes the keys, the binning_operator and optionally whether to
  // sort within bins (default false)
  template <typename ExecutionSpace>
  BinSort(const ExecutionSpace& exec, const_key_view_type keys_,
          int range_begin_, int range_end_, BinSortOp bin_op_,
          bool sort_within_bins_ = false)
      : keys(keys_),
        keys_rnd(keys_),
        bin_op(bin_op_),
        bin_offsets(),
        bin_count_atomic(),
        bin_count_const(),
        sort_order(),
        range_begin(range_begin_),
        range_end(range_end_),
        sort_within_bins(sort_within_bins_) {
    static_assert(
        Kokkos::SpaceAccessibility<ExecutionSpace,
                                   typename Space::memory_space>::accessible,
        "The provided execution space must be able to access the memory space "
        "BinSort was initialized with!");
    if (bin_op.max_bins() <= 0)
      Kokkos::abort(
          "The number of bins in the BinSortOp object must be greater than 0!");
    bin_count_atomic = Kokkos::View<int*, Space>(
        "Kokkos::SortImpl::BinSortFunctor::bin_count", bin_op.max_bins());
    bin_count_const = bin_count_atomic;
    bin_offsets =
        offset_type(view_alloc(exec, WithoutInitializing,
                               "Kokkos::SortImpl::BinSortFunctor::bin_offsets"),
                    bin_op.max_bins());
    sort_order =
        offset_type(view_alloc(exec, WithoutInitializing,
                               "Kokkos::SortImpl::BinSortFunctor::sort_order"),
                    range_end - range_begin);
  }

  BinSort(const_key_view_type keys_, int range_begin_, int range_end_,
          BinSortOp bin_op_, bool sort_within_bins_ = false)
      : BinSort(exec_space{}, keys_, range_begin_, range_end_, bin_op_,
                sort_within_bins_) {}

  template <typename ExecutionSpace>
  BinSort(const ExecutionSpace& exec, const_key_view_type keys_,
          BinSortOp bin_op_, bool sort_within_bins_ = false)
      : BinSort(exec, keys_, 0, keys_.extent(0), bin_op_, sort_within_bins_) {}

  BinSort(const_key_view_type keys_, BinSortOp bin_op_,
          bool sort_within_bins_ = false)
      : BinSort(exec_space{}, keys_, bin_op_, sort_within_bins_) {}

  //----------------------------------------
  // Create the permutation vector, the bin_offset array and the bin_count
  // array. Can be called again if keys changed
  template <class ExecutionSpace>
  void create_permute_vector(const ExecutionSpace& exec) {
    static_assert(
        Kokkos::SpaceAccessibility<ExecutionSpace,
                                   typename Space::memory_space>::accessible,
        "The provided execution space must be able to access the memory space "
        "BinSort was initialized with!");

    const size_t len = range_end - range_begin;
    Kokkos::parallel_for(
        "Kokkos::Sort::BinCount",
        Kokkos::RangePolicy<ExecutionSpace, bin_count_tag>(exec, 0, len),
        *this);
    Kokkos::parallel_scan("Kokkos::Sort::BinOffset",
                          Kokkos::RangePolicy<ExecutionSpace, bin_offset_tag>(
                              exec, 0, bin_op.max_bins()),
                          *this);

    Kokkos::deep_copy(exec, bin_count_atomic, 0);
    Kokkos::parallel_for(
        "Kokkos::Sort::BinBinning",
        Kokkos::RangePolicy<ExecutionSpace, bin_binning_tag>(exec, 0, len),
        *this);

    if (sort_within_bins)
      Kokkos::parallel_for(
          "Kokkos::Sort::BinSort",
          Kokkos::RangePolicy<ExecutionSpace, bin_sort_bins_tag>(
              exec, 0, bin_op.max_bins()),
          *this);
  }

  // Create the permutation vector, the bin_offset array and the bin_count
  // array. Can be called again if keys changed
  void create_permute_vector() {
    Kokkos::fence("Kokkos::Binsort::create_permute_vector: before");
    exec_space e{};
    create_permute_vector(e);
    e.fence("Kokkos::Binsort::create_permute_vector: after");
  }

  // Sort a subset of a view with respect to the first dimension using the
  // permutation array
  template <class ExecutionSpace, class ValuesViewType>
  void sort(const ExecutionSpace& exec, ValuesViewType const& values,
            int values_range_begin, int values_range_end) const {
    if (values.extent(0) == 0) {
      return;
    }

    static_assert(
        Kokkos::SpaceAccessibility<ExecutionSpace,
                                   typename Space::memory_space>::accessible,
        "The provided execution space must be able to access the memory space "
        "BinSort was initialized with!");
    static_assert(
        Kokkos::SpaceAccessibility<
            ExecutionSpace, typename ValuesViewType::memory_space>::accessible,
        "The provided execution space must be able to access the memory space "
        "of the View argument!");

    const size_t len        = range_end - range_begin;
    const size_t values_len = values_range_end - values_range_begin;
    if (len != values_len) {
      Kokkos::abort(
          "BinSort::sort: values range length != permutation vector length");
    }

    using scratch_view_type =
        Kokkos::View<typename ValuesViewType::data_type,
                     typename ValuesViewType::device_type>;
    scratch_view_type sorted_values(
        view_alloc(exec, WithoutInitializing,
                   "Kokkos::SortImpl::BinSortFunctor::sorted_values"),
        values.rank_dynamic > 0 ? len : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 1 ? values.extent(1)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 2 ? values.extent(2)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 3 ? values.extent(3)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 4 ? values.extent(4)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 5 ? values.extent(5)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 6 ? values.extent(6)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        values.rank_dynamic > 7 ? values.extent(7)
                                : KOKKOS_IMPL_CTOR_DEFAULT_ARG);

    {
      copy_permute_functor<scratch_view_type /* DstViewType */
                           ,
                           offset_type /* PermuteViewType */
                           ,
                           ValuesViewType /* SrcViewType */
                           >
          functor(sorted_values, sort_order, values,
                  values_range_begin - range_begin);

      parallel_for("Kokkos::Sort::CopyPermute",
                   Kokkos::RangePolicy<ExecutionSpace>(exec, 0, len), functor);
    }

    {
      copy_functor<ValuesViewType, scratch_view_type> functor(
          values, range_begin, sorted_values);

      parallel_for("Kokkos::Sort::Copy",
                   Kokkos::RangePolicy<ExecutionSpace>(exec, 0, len), functor);
    }
  }

  // Sort a subset of a view with respect to the first dimension using the
  // permutation array
  template <class ValuesViewType>
  void sort(ValuesViewType const& values, int values_range_begin,
            int values_range_end) const {
    Kokkos::fence("Kokkos::Binsort::sort: before");
    exec_space exec;
    sort(exec, values, values_range_begin, values_range_end);
    exec.fence("Kokkos::BinSort:sort: after");
  }

  template <class ExecutionSpace, class ValuesViewType>
  void sort(ExecutionSpace const& exec, ValuesViewType const& values) const {
    this->sort(exec, values, 0, /*values.extent(0)*/ range_end - range_begin);
  }

  template <class ValuesViewType>
  void sort(ValuesViewType const& values) const {
    this->sort(values, 0, /*values.extent(0)*/ range_end - range_begin);
  }

  // Get the permutation vector
  KOKKOS_INLINE_FUNCTION
  offset_type get_permute_vector() const { return sort_order; }

  // Get the start offsets for each bin
  KOKKOS_INLINE_FUNCTION
  offset_type get_bin_offsets() const { return bin_offsets; }

  // Get the count for each bin
  KOKKOS_INLINE_FUNCTION
  bin_count_type get_bin_count() const { return bin_count_const; }

 public:
  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_count_tag& /*tag*/, const int i) const {
    const int j = range_begin + i;
    bin_count_atomic(bin_op.bin(keys, j))++;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_offset_tag& /*tag*/, const int i,
                  value_type& offset, const bool& final) const {
    if (final) {
      bin_offsets(i) = offset;
    }
    offset += bin_count_const(i);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_binning_tag& /*tag*/, const int i) const {
    const int j     = range_begin + i;
    const int bin   = bin_op.bin(keys, j);
    const int count = bin_count_atomic(bin)++;

    sort_order(bin_offsets(bin) + count) = j;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_sort_bins_tag& /*tag*/, const int i) const {
    auto bin_size = bin_count_const(i);
    if (bin_size <= 1) return;
    constexpr bool use_std_sort =
        std::is_same_v<typename exec_space::memory_space, HostSpace>;
    int lower_bound = bin_offsets(i);
    int upper_bound = lower_bound + bin_size;
    // Switching to std::sort for more than 10 elements has been found
    // reasonable experimentally.
    if (use_std_sort && bin_size > 10) {
      if constexpr (use_std_sort) {
        std::sort(&sort_order(lower_bound), &sort_order(upper_bound),
                  [this](int p, int q) { return bin_op(keys_rnd, p, q); });
      }
    } else {
      for (int k = lower_bound + 1; k < upper_bound; ++k) {
        int old_idx = sort_order(k);
        int j       = k - 1;
        while (j >= lower_bound) {
          int new_idx = sort_order(j);
          if (!bin_op(keys_rnd, old_idx, new_idx)) break;
          sort_order(j + 1) = new_idx;
          --j;
        }
        sort_order(j + 1) = old_idx;
      }
    }
  }
};

//----------------------------------------------------------------------------

template <class KeyViewType>
struct BinOp1D {
  int max_bins_ = {};
  double mul_   = {};
  double min_   = {};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED BinOp1D() = default;
#else
  BinOp1D() = delete;
#endif

  // Construct BinOp with number of bins, minimum value and maximum value
  BinOp1D(int max_bins__, typename KeyViewType::const_value_type min,
          typename KeyViewType::const_value_type max)
      : max_bins_(max_bins__ + 1),
        // Cast to double to avoid possible overflow when using integer
        mul_(static_cast<double>(max_bins__) /
             (static_cast<double>(max) - static_cast<double>(min))),
        min_(static_cast<double>(min)) {
    // For integral types the number of bins may be larger than the range
    // in which case we can exactly have one unique value per bin
    // and then don't need to sort bins.
    if (std::is_integral<typename KeyViewType::const_value_type>::value &&
        (static_cast<double>(max) - static_cast<double>(min)) <=
            static_cast<double>(max_bins__)) {
      mul_ = 1.;
    }
  }

  // Determine bin index from key value
  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    return static_cast<int>(mul_ * (static_cast<double>(keys(i)) - min_));
  }

  // Return maximum bin index + 1
  KOKKOS_INLINE_FUNCTION
  int max_bins() const { return max_bins_; }

  // Compare to keys within a bin if true new_val will be put before old_val
  template <class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION bool operator()(ViewType& keys, iType1& i1,
                                         iType2& i2) const {
    return keys(i1) < keys(i2);
  }
};

template <class KeyViewType>
struct BinOp3D {
  int max_bins_[3] = {};
  double mul_[3]   = {};
  double min_[3]   = {};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED BinOp3D() = default;
#else
  BinOp3D() = delete;
#endif

  BinOp3D(int max_bins__[], typename KeyViewType::const_value_type min[],
          typename KeyViewType::const_value_type max[]) {
    max_bins_[0] = max_bins__[0];
    max_bins_[1] = max_bins__[1];
    max_bins_[2] = max_bins__[2];
    mul_[0]      = static_cast<double>(max_bins__[0]) /
              (static_cast<double>(max[0]) - static_cast<double>(min[0]));
    mul_[1] = static_cast<double>(max_bins__[1]) /
              (static_cast<double>(max[1]) - static_cast<double>(min[1]));
    mul_[2] = static_cast<double>(max_bins__[2]) /
              (static_cast<double>(max[2]) - static_cast<double>(min[2]));
    min_[0] = static_cast<double>(min[0]);
    min_[1] = static_cast<double>(min[1]);
    min_[2] = static_cast<double>(min[2]);
  }

  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    return int((((int(mul_[0] * (keys(i, 0) - min_[0])) * max_bins_[1]) +
                 int(mul_[1] * (keys(i, 1) - min_[1]))) *
                max_bins_[2]) +
               int(mul_[2] * (keys(i, 2) - min_[2])));
  }

  KOKKOS_INLINE_FUNCTION
  int max_bins() const { return max_bins_[0] * max_bins_[1] * max_bins_[2]; }

  template <class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION bool operator()(ViewType& keys, iType1& i1,
                                         iType2& i2) const {
    if (keys(i1, 0) > keys(i2, 0))
      return true;
    else if (keys(i1, 0) == keys(i2, 0)) {
      if (keys(i1, 1) > keys(i2, 1))
        return true;
      else if (keys(i1, 1) == keys(i2, 1)) {
        if (keys(i1, 2) > keys(i2, 2)) return true;
      }
    }
    return false;
  }
};

namespace Impl {

template <class ViewType>
struct min_max_functor {
  using minmax_scalar =
      Kokkos::MinMaxScalar<typename ViewType::non_const_value_type>;

  ViewType view;
  min_max_functor(const ViewType& view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t& i, minmax_scalar& minmax) const {
    if (view(i) < minmax.min_val) minmax.min_val = view(i);
    if (view(i) > minmax.max_val) minmax.max_val = view(i);
  }
};

}  // namespace Impl

template <class ExecutionSpace, class DataType, class... Properties>
std::enable_if_t<(Kokkos::is_execution_space<ExecutionSpace>::value) &&
                 (!SpaceAccessibility<
                     HostSpace, typename Kokkos::View<DataType, Properties...>::
                                    memory_space>::accessible)>
sort(const ExecutionSpace& exec,
     const Kokkos::View<DataType, Properties...>& view) {
  if (view.extent(0) == 0) {
    return;
  }

  using ViewType = Kokkos::View<DataType, Properties...>;
  using CompType = BinOp1D<ViewType>;

  Kokkos::MinMaxScalar<typename ViewType::non_const_value_type> result;
  Kokkos::MinMax<typename ViewType::non_const_value_type> reducer(result);
  parallel_reduce("Kokkos::Sort::FindExtent",
                  Kokkos::RangePolicy<typename ViewType::execution_space>(
                      exec, 0, view.extent(0)),
                  Impl::min_max_functor<ViewType>(view), reducer);
  if (result.min_val == result.max_val) return;
  // For integral types the number of bins may be larger than the range
  // in which case we can exactly have one unique value per bin
  // and then don't need to sort bins.
  bool sort_in_bins = true;
  // TODO: figure out better max_bins then this ...
  int64_t max_bins = view.extent(0) / 2;
  if (std::is_integral<typename ViewType::non_const_value_type>::value) {
    // Cast to double to avoid possible overflow when using integer
    auto const max_val = static_cast<double>(result.max_val);
    auto const min_val = static_cast<double>(result.min_val);
    // using 10M as the cutoff for special behavior (roughly 40MB for the count
    // array)
    if ((max_val - min_val) < 10000000) {
      max_bins     = max_val - min_val + 1;
      sort_in_bins = false;
    }
  }
  if (std::is_floating_point<typename ViewType::non_const_value_type>::value) {
    KOKKOS_ASSERT(std::isfinite(static_cast<double>(result.max_val) -
                                static_cast<double>(result.min_val)));
  }

  BinSort<ViewType, CompType> bin_sort(
      view, CompType(max_bins, result.min_val, result.max_val), sort_in_bins);
  bin_sort.create_permute_vector(exec);
  bin_sort.sort(exec, view);
}

#if defined(KOKKOS_ENABLE_ONEDPL)
template <class DataType, class... Properties>
void sort(const Experimental::SYCL& space,
          const Kokkos::View<DataType, Properties...>& view) {
  if (view.extent(0) == 0) {
    return;
  }

  using ViewType = Kokkos::View<DataType, Properties...>;
  static_assert(SpaceAccessibility<Experimental::SYCL,
                                   typename ViewType::memory_space>::accessible,
                "SYCL execution space is not able to access the memory space "
                "of the View argument!");

  auto queue  = space.sycl_queue();
  auto policy = oneapi::dpl::execution::make_device_policy(queue);

  // Can't use Experimental::begin/end here since the oneDPL then assumes that
  // the data is on the host.
  static_assert(
      ViewType::rank == 1 &&
          (std::is_same<typename ViewType::array_layout, LayoutRight>::value ||
           std::is_same<typename ViewType::array_layout, LayoutLeft>::value),
      "SYCL sort only supports contiguous 1D Views.");
  const int n = view.extent(0);
  oneapi::dpl::sort(policy, view.data(), view.data() + n);
}
#endif

template <class ExecutionSpace, class DataType, class... Properties>
std::enable_if_t<(Kokkos::is_execution_space<ExecutionSpace>::value) &&
                 (SpaceAccessibility<
                     HostSpace, typename Kokkos::View<DataType, Properties...>::
                                    memory_space>::accessible)>
sort(const ExecutionSpace&, const Kokkos::View<DataType, Properties...>& view) {
  if (view.extent(0) == 0) {
    return;
  }
  auto first = Experimental::begin(view);
  auto last  = Experimental::end(view);
  std::sort(first, last);
}

#if defined(KOKKOS_ENABLE_CUDA)
template <class DataType, class... Properties>
void sort(const Cuda& space,
          const Kokkos::View<DataType, Properties...>& view) {
  if (view.extent(0) == 0) {
    return;
  }
  const auto exec = thrust::cuda::par.on(space.cuda_stream());
  auto first      = Experimental::begin(view);
  auto last       = Experimental::end(view);
  thrust::sort(exec, first, last);
}
#endif

template <class ViewType>
void sort(ViewType const& view) {
  Kokkos::fence("Kokkos::sort: before");

  if (view.extent(0) == 0) {
    return;
  }

  typename ViewType::execution_space exec;
  sort(exec, view);
  exec.fence("Kokkos::sort: fence after sorting");
}

template <class ExecutionSpace, class ViewType>
std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value> sort(
    const ExecutionSpace& exec, ViewType view, size_t const begin,
    size_t const end) {
  if (view.extent(0) == 0) {
    return;
  }

  using range_policy = Kokkos::RangePolicy<typename ViewType::execution_space>;
  using CompType     = BinOp1D<ViewType>;

  Kokkos::MinMaxScalar<typename ViewType::non_const_value_type> result;
  Kokkos::MinMax<typename ViewType::non_const_value_type> reducer(result);

  parallel_reduce("Kokkos::Sort::FindExtent", range_policy(exec, begin, end),
                  Impl::min_max_functor<ViewType>(view), reducer);

  if (result.min_val == result.max_val) return;

  BinSort<ViewType, CompType> bin_sort(
      exec, view, begin, end,
      CompType((end - begin) / 2, result.min_val, result.max_val), true);

  bin_sort.create_permute_vector(exec);
  bin_sort.sort(exec, view, begin, end);
}

template <class ViewType>
void sort(ViewType view, size_t const begin, size_t const end) {
  Kokkos::fence("Kokkos::sort: before");

  if (view.extent(0) == 0) {
    return;
  }

  typename ViewType::execution_space exec;
  sort(exec, view, begin, end);
  exec.fence("Kokkos::Sort: fence after sorting");
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_SORT
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_SORT
#endif
#endif
