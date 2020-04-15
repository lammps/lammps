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

#ifndef KOKKOS_SORT_HPP_
#define KOKKOS_SORT_HPP_

#include <Kokkos_Core.hpp>

#include <algorithm>

namespace Kokkos {

namespace Impl {

template <class DstViewType, class SrcViewType, int Rank = DstViewType::Rank>
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
    typedef typename SrcViewType::const_type src_view_type;

    typedef Impl::CopyOp<DstViewType, src_view_type> copy_op;

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

    typedef typename std::conditional<
        Kokkos::is_view<SrcViewType>::value,
        Kokkos::View<typename SrcViewType::const_data_type,
                     typename SrcViewType::array_layout,
                     typename SrcViewType::device_type,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess> >,
        typename SrcViewType::const_type>::type src_view_type;

    typedef typename PermuteViewType::const_type perm_view_type;

    typedef Impl::CopyOp<DstViewType, src_view_type> copy_op;

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

  typedef typename Space::execution_space execution_space;
  typedef BinSortOp bin_op_type;

  struct bin_count_tag {};
  struct bin_offset_tag {};
  struct bin_binning_tag {};
  struct bin_sort_bins_tag {};

 public:
  typedef SizeType size_type;
  typedef size_type value_type;

  typedef Kokkos::View<size_type*, Space> offset_type;
  typedef Kokkos::View<const int*, Space> bin_count_type;

  typedef typename KeyViewType::const_type const_key_view_type;

  // If a Kokkos::View then can generate constant random access
  // otherwise can only use the constant type.

  typedef typename std::conditional<
      Kokkos::is_view<KeyViewType>::value,
      Kokkos::View<typename KeyViewType::const_data_type,
                   typename KeyViewType::array_layout,
                   typename KeyViewType::device_type,
                   Kokkos::MemoryTraits<Kokkos::RandomAccess> >,
      const_key_view_type>::type const_rnd_key_view_type;

  typedef typename KeyViewType::non_const_value_type non_const_key_scalar;
  typedef typename KeyViewType::const_value_type const_key_scalar;

  typedef Kokkos::View<int*, Space, Kokkos::MemoryTraits<Kokkos::Atomic> >
      bin_count_atomic_type;

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
  BinSort() {}

  //----------------------------------------
  // Constructor: takes the keys, the binning_operator and optionally whether to
  // sort within bins (default false)
  BinSort(const_key_view_type keys_, int range_begin_, int range_end_,
          BinSortOp bin_op_, bool sort_within_bins_ = false)
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
    bin_count_atomic = Kokkos::View<int*, Space>(
        "Kokkos::SortImpl::BinSortFunctor::bin_count", bin_op.max_bins());
    bin_count_const = bin_count_atomic;
    bin_offsets =
        offset_type(ViewAllocateWithoutInitializing(
                        "Kokkos::SortImpl::BinSortFunctor::bin_offsets"),
                    bin_op.max_bins());
    sort_order =
        offset_type(ViewAllocateWithoutInitializing(
                        "Kokkos::SortImpl::BinSortFunctor::sort_order"),
                    range_end - range_begin);
  }

  BinSort(const_key_view_type keys_, BinSortOp bin_op_,
          bool sort_within_bins_ = false)
      : BinSort(keys_, 0, keys_.extent(0), bin_op_, sort_within_bins_) {}

  //----------------------------------------
  // Create the permutation vector, the bin_offset array and the bin_count
  // array. Can be called again if keys changed
  void create_permute_vector() {
    const size_t len = range_end - range_begin;
    Kokkos::parallel_for(
        "Kokkos::Sort::BinCount",
        Kokkos::RangePolicy<execution_space, bin_count_tag>(0, len), *this);
    Kokkos::parallel_scan("Kokkos::Sort::BinOffset",
                          Kokkos::RangePolicy<execution_space, bin_offset_tag>(
                              0, bin_op.max_bins()),
                          *this);

    Kokkos::deep_copy(bin_count_atomic, 0);
    Kokkos::parallel_for(
        "Kokkos::Sort::BinBinning",
        Kokkos::RangePolicy<execution_space, bin_binning_tag>(0, len), *this);

    if (sort_within_bins)
      Kokkos::parallel_for(
          "Kokkos::Sort::BinSort",
          Kokkos::RangePolicy<execution_space, bin_sort_bins_tag>(
              0, bin_op.max_bins()),
          *this);
  }

  // Sort a subset of a view with respect to the first dimension using the
  // permutation array
  template <class ValuesViewType>
  void sort(ValuesViewType const& values, int values_range_begin,
            int values_range_end) const {
    typedef Kokkos::View<typename ValuesViewType::data_type,
                         typename ValuesViewType::array_layout,
                         typename ValuesViewType::device_type>
        scratch_view_type;

    const size_t len        = range_end - range_begin;
    const size_t values_len = values_range_end - values_range_begin;
    if (len != values_len) {
      Kokkos::abort(
          "BinSort::sort: values range length != permutation vector length");
    }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    scratch_view_type sorted_values(
        ViewAllocateWithoutInitializing(
            "Kokkos::SortImpl::BinSortFunctor::sorted_values"),
        len, values.extent(1), values.extent(2), values.extent(3),
        values.extent(4), values.extent(5), values.extent(6), values.extent(7));
#else
    scratch_view_type sorted_values(
        ViewAllocateWithoutInitializing(
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
#endif

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
                   Kokkos::RangePolicy<execution_space>(0, len), functor);
    }

    {
      copy_functor<ValuesViewType, scratch_view_type> functor(
          values, range_begin, sorted_values);

      parallel_for("Kokkos::Sort::Copy",
                   Kokkos::RangePolicy<execution_space>(0, len), functor);
    }

    Kokkos::fence();
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
  void operator()(const bin_count_tag& tag, const int& i) const {
    const int j = range_begin + i;
    bin_count_atomic(bin_op.bin(keys, j))++;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_offset_tag& tag, const int& i, value_type& offset,
                  const bool& final) const {
    if (final) {
      bin_offsets(i) = offset;
    }
    offset += bin_count_const(i);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_binning_tag& tag, const int& i) const {
    const int j     = range_begin + i;
    const int bin   = bin_op.bin(keys, j);
    const int count = bin_count_atomic(bin)++;

    sort_order(bin_offsets(bin) + count) = j;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const bin_sort_bins_tag& tag, const int& i) const {
    auto bin_size = bin_count_const(i);
    if (bin_size <= 1) return;
    int upper_bound = bin_offsets(i) + bin_size;
    bool sorted     = false;
    while (!sorted) {
      sorted      = true;
      int old_idx = sort_order(bin_offsets(i));
      int new_idx;
      for (int k = bin_offsets(i) + 1; k < upper_bound; k++) {
        new_idx = sort_order(k);

        if (!bin_op(keys_rnd, old_idx, new_idx)) {
          sort_order(k - 1) = new_idx;
          sort_order(k)     = old_idx;
          sorted            = false;
        } else {
          old_idx = new_idx;
        }
      }
      upper_bound--;
    }
  }
};

//----------------------------------------------------------------------------

template <class KeyViewType>
struct BinOp1D {
  int max_bins_;
  double mul_;
  typename KeyViewType::const_value_type range_;
  typename KeyViewType::const_value_type min_;

  BinOp1D()
      : max_bins_(0),
        mul_(0.0),
        range_(typename KeyViewType::const_value_type()),
        min_(typename KeyViewType::const_value_type()) {}

  // Construct BinOp with number of bins, minimum value and maxuimum value
  BinOp1D(int max_bins__, typename KeyViewType::const_value_type min,
          typename KeyViewType::const_value_type max)
      : max_bins_(max_bins__ + 1),
        mul_(1.0 * max_bins__ / (max - min)),
        range_(max - min),
        min_(min) {}

  // Determine bin index from key value
  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    return int(mul_ * (keys(i) - min_));
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
  int max_bins_[3];
  double mul_[3];
  typename KeyViewType::non_const_value_type range_[3];
  typename KeyViewType::non_const_value_type min_[3];

  BinOp3D() {}

  BinOp3D(int max_bins__[], typename KeyViewType::const_value_type min[],
          typename KeyViewType::const_value_type max[]) {
    max_bins_[0] = max_bins__[0];
    max_bins_[1] = max_bins__[1];
    max_bins_[2] = max_bins__[2];
    mul_[0]      = 1.0 * max_bins__[0] / (max[0] - min[0]);
    mul_[1]      = 1.0 * max_bins__[1] / (max[1] - min[1]);
    mul_[2]      = 1.0 * max_bins__[2] / (max[2] - min[2]);
    range_[0]    = max[0] - min[0];
    range_[1]    = max[1] - min[1];
    range_[2]    = max[2] - min[2];
    min_[0]      = min[0];
    min_[1]      = min[1];
    min_[2]      = min[2];
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
bool try_std_sort(ViewType view) {
  bool possible    = true;
  size_t stride[8] = {view.stride_0(), view.stride_1(), view.stride_2(),
                      view.stride_3(), view.stride_4(), view.stride_5(),
                      view.stride_6(), view.stride_7()};
  possible         = possible &&
             std::is_same<typename ViewType::memory_space, HostSpace>::value;
  possible = possible && (ViewType::Rank == 1);
  possible = possible && (stride[0] == 1);
  if (possible) {
    std::sort(view.data(), view.data() + view.extent(0));
  }
  return possible;
}

template <class ViewType>
struct min_max_functor {
  typedef Kokkos::MinMaxScalar<typename ViewType::non_const_value_type>
      minmax_scalar;

  ViewType view;
  min_max_functor(const ViewType& view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t& i, minmax_scalar& minmax) const {
    if (view(i) < minmax.min_val) minmax.min_val = view(i);
    if (view(i) > minmax.max_val) minmax.max_val = view(i);
  }
};

}  // namespace Impl

template <class ViewType>
void sort(ViewType const& view, bool const always_use_kokkos_sort = false) {
  if (!always_use_kokkos_sort) {
    if (Impl::try_std_sort(view)) return;
  }
  typedef BinOp1D<ViewType> CompType;

  Kokkos::MinMaxScalar<typename ViewType::non_const_value_type> result;
  Kokkos::MinMax<typename ViewType::non_const_value_type> reducer(result);
  parallel_reduce("Kokkos::Sort::FindExtent",
                  Kokkos::RangePolicy<typename ViewType::execution_space>(
                      0, view.extent(0)),
                  Impl::min_max_functor<ViewType>(view), reducer);
  if (result.min_val == result.max_val) return;
  BinSort<ViewType, CompType> bin_sort(
      view, CompType(view.extent(0) / 2, result.min_val, result.max_val), true);
  bin_sort.create_permute_vector();
  bin_sort.sort(view);
}

template <class ViewType>
void sort(ViewType view, size_t const begin, size_t const end) {
  typedef Kokkos::RangePolicy<typename ViewType::execution_space> range_policy;
  typedef BinOp1D<ViewType> CompType;

  Kokkos::MinMaxScalar<typename ViewType::non_const_value_type> result;
  Kokkos::MinMax<typename ViewType::non_const_value_type> reducer(result);

  parallel_reduce("Kokkos::Sort::FindExtent", range_policy(begin, end),
                  Impl::min_max_functor<ViewType>(view), reducer);

  if (result.min_val == result.max_val) return;

  BinSort<ViewType, CompType> bin_sort(
      view, begin, end,
      CompType((end - begin) / 2, result.min_val, result.max_val), true);

  bin_sort.create_permute_vector();
  bin_sort.sort(view, begin, end);
}

}  // namespace Kokkos

#endif
