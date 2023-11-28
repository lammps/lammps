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

#ifndef KOKKOS_BIN_SORT_PUBLIC_API_HPP_
#define KOKKOS_BIN_SORT_PUBLIC_API_HPP_

#include "Kokkos_BinOpsPublicAPI.hpp"
#include "impl/Kokkos_CopyOpsForBinSortImpl.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>

namespace Kokkos {

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
#if !defined(KOKKOS_COMPILER_NVHPC) || (KOKKOS_COMPILER_NVHPC >= 230700)
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
      KOKKOS_IF_ON_HOST(
          (std::sort(&sort_order(lower_bound), &sort_order(upper_bound),
                     [this](int p, int q) { return bin_op(keys_rnd, p, q); });))
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

}  // namespace Kokkos
#endif
