/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


#ifndef KOKKOS_SORT_HPP_
#define KOKKOS_SORT_HPP_

#include <Kokkos_Core.hpp>

#include <algorithm>

namespace Kokkos {

  namespace SortImpl {

  template<class ValuesViewType, int Rank=ValuesViewType::Rank>
  struct CopyOp;

  template<class ValuesViewType>
  struct CopyOp<ValuesViewType,1> {
    template<class DstType, class SrcType>
    KOKKOS_INLINE_FUNCTION
    static void copy(DstType& dst, size_t i_dst,
                     SrcType& src, size_t i_src ) {
      dst(i_dst) = src(i_src);
    }
  };

  template<class ValuesViewType>
  struct CopyOp<ValuesViewType,2> {
    template<class DstType, class SrcType>
    KOKKOS_INLINE_FUNCTION
    static void copy(DstType& dst, size_t i_dst,
                     SrcType& src, size_t i_src ) {
      for(int j = 0;j< (int) dst.dimension_1(); j++)
        dst(i_dst,j) = src(i_src,j);
    }
  };

  template<class ValuesViewType>
  struct CopyOp<ValuesViewType,3> {
    template<class DstType, class SrcType>
    KOKKOS_INLINE_FUNCTION
    static void copy(DstType& dst, size_t i_dst,
                     SrcType& src, size_t i_src ) {
      for(int j = 0; j<dst.dimension_1(); j++)
        for(int k = 0; k<dst.dimension_2(); k++)
          dst(i_dst,j,k) = src(i_src,j,k);
    }
  };
  }

template<class KeyViewType, class BinSortOp, class ExecutionSpace = typename KeyViewType::execution_space,
         class SizeType = typename KeyViewType::memory_space::size_type>
class BinSort {


public:
  template<class ValuesViewType, class PermuteViewType, class CopyOp>
  struct bin_sort_sort_functor {
    typedef ExecutionSpace execution_space;
    typedef typename ValuesViewType::non_const_type values_view_type;
    typedef typename ValuesViewType::const_type const_values_view_type;
    Kokkos::View<typename values_view_type::const_data_type,typename values_view_type::array_layout,
                 typename values_view_type::memory_space,Kokkos::MemoryTraits<Kokkos::RandomAccess> > values;
    values_view_type sorted_values;
    typename PermuteViewType::const_type sort_order;
    bin_sort_sort_functor(const_values_view_type values_, values_view_type  sorted_values_, PermuteViewType sort_order_):
       values(values_),sorted_values(sorted_values_),sort_order(sort_order_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int& i)  const {
      //printf("Sort: %i %i\n",i,sort_order(i));
      CopyOp::copy(sorted_values,i,values,sort_order(i));
    }
  };

  typedef ExecutionSpace execution_space;
  typedef BinSortOp bin_op_type;

  struct bin_count_tag {};
  struct bin_offset_tag {};
  struct bin_binning_tag {};
  struct bin_sort_bins_tag {};

public:
  typedef SizeType size_type;
  typedef size_type value_type;

  typedef Kokkos::View<size_type*, execution_space> offset_type;
  typedef Kokkos::View<const int*, execution_space> bin_count_type;


  typedef Kokkos::View<typename KeyViewType::const_data_type,
                       typename KeyViewType::array_layout,
                       typename KeyViewType::memory_space> const_key_view_type;
  typedef Kokkos::View<typename KeyViewType::const_data_type,
                       typename KeyViewType::array_layout,
                       typename KeyViewType::memory_space,
                       Kokkos::MemoryTraits<Kokkos::RandomAccess> > const_rnd_key_view_type;

  typedef typename KeyViewType::non_const_value_type non_const_key_scalar;
  typedef typename KeyViewType::const_value_type     const_key_scalar;

private:
  const_key_view_type keys;
  const_rnd_key_view_type keys_rnd;

public:
  BinSortOp bin_op;

  offset_type bin_offsets;

  Kokkos::View<int*, ExecutionSpace, Kokkos::MemoryTraits<Kokkos::Atomic> > bin_count_atomic;
  bin_count_type bin_count_const;

  offset_type sort_order;

  bool sort_within_bins;

public:

  // Constructor: takes the keys, the binning_operator and optionally whether to sort within bins (default false)
  BinSort(const_key_view_type keys_, BinSortOp bin_op_,
          bool sort_within_bins_ = false)
     :keys(keys_),keys_rnd(keys_), bin_op(bin_op_) {

    bin_count_atomic = Kokkos::View<int*, ExecutionSpace >("Kokkos::SortImpl::BinSortFunctor::bin_count",bin_op.max_bins());
    bin_count_const =  bin_count_atomic;
    bin_offsets =      offset_type("Kokkos::SortImpl::BinSortFunctor::bin_offsets",bin_op.max_bins());
    sort_order =       offset_type("PermutationVector",keys.dimension_0());
    sort_within_bins = sort_within_bins_;
  }

  // Create the permutation vector, the bin_offset array and the bin_count array. Can be called again if keys changed
  void create_permute_vector() {
    Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_count_tag>    (0,keys.dimension_0()),*this);
    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecutionSpace,bin_offset_tag>   (0,bin_op.max_bins()) ,*this);

    Kokkos::deep_copy(bin_count_atomic,0);
    Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_binning_tag>  (0,keys.dimension_0()),*this);

    if(sort_within_bins)
      Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_sort_bins_tag>(0,bin_op.max_bins()) ,*this);
  }

  // Sort a view with respect ot the first dimension using the permutation array
  template<class ValuesViewType>
  void sort(ValuesViewType values) {
    ValuesViewType sorted_values = ValuesViewType("Copy",
           values.dimension_0(),
           values.dimension_1(),
           values.dimension_2(),
           values.dimension_3(),
           values.dimension_4(),
           values.dimension_5(),
           values.dimension_6(),
           values.dimension_7());

    parallel_for(values.dimension_0(),
        bin_sort_sort_functor<ValuesViewType, offset_type,
                              SortImpl::CopyOp<ValuesViewType> >(values,sorted_values,sort_order));

    deep_copy(values,sorted_values);
  }

  // Get the permutation vector
  KOKKOS_INLINE_FUNCTION
  offset_type get_permute_vector() const { return sort_order;}

  // Get the start offsets for each bin
  KOKKOS_INLINE_FUNCTION
  offset_type get_bin_offsets() const { return bin_offsets;}

  // Get the count for each bin
  KOKKOS_INLINE_FUNCTION
  bin_count_type get_bin_count() const {return bin_count_const;}

public:
  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_count_tag& tag, const int& i) const {
    bin_count_atomic(bin_op.bin(keys,i))++;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_offset_tag& tag, const int& i, value_type& offset, const bool& final)  const {
    if(final) {
      bin_offsets(i) = offset;
    }
    offset+=bin_count_const(i);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_binning_tag& tag, const int& i)  const {
    const int bin = bin_op.bin(keys,i);
    const int count = bin_count_atomic(bin)++;

    sort_order(bin_offsets(bin) + count) = i;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_sort_bins_tag& tag, const int&i )  const {
    bool sorted = false;
    int upper_bound = bin_offsets(i)+bin_count_const(i);
    while(!sorted) {
      sorted = true;
      int old_idx = sort_order(bin_offsets(i));
      int new_idx;
      for(int k=bin_offsets(i)+1; k<upper_bound; k++) {
        new_idx = sort_order(k);

        if(!bin_op(keys_rnd,old_idx,new_idx)) {
          sort_order(k-1) = new_idx;
          sort_order(k) = old_idx;
          sorted = false;
        } else {
          old_idx = new_idx;
        }
      }
      upper_bound--;
    }
  }
};

namespace SortImpl {

template<class KeyViewType>
struct DefaultBinOp1D {
  const int max_bins_;
  const double mul_;
  typename KeyViewType::const_value_type range_;
  typename KeyViewType::const_value_type min_;

  //Construct BinOp with number of bins, minimum value and maxuimum value
  DefaultBinOp1D(int max_bins, typename KeyViewType::const_value_type min,
                               typename KeyViewType::const_value_type max )
     :max_bins_(max_bins+1),mul_(1.0*max_bins/(max-min)),range_(max-min),min_(min) {}

  //Determine bin index from key value
  template<class ViewType>
  KOKKOS_INLINE_FUNCTION
  int bin(ViewType& keys, const int& i) const {
    return int(mul_*(keys(i)-min_));
  }

  //Return maximum bin index + 1
  KOKKOS_INLINE_FUNCTION
  int max_bins() const {
    return max_bins_;
  }

  //Compare to keys within a bin if true new_val will be put before old_val
  template<class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION
  bool operator()(ViewType& keys, iType1& i1, iType2& i2) const {
    return keys(i1)<keys(i2);
  }
};

template<class KeyViewType>
struct DefaultBinOp3D {
  int max_bins_[3];
  double mul_[3];
  typename KeyViewType::non_const_value_type range_[3];
  typename KeyViewType::non_const_value_type min_[3];

  DefaultBinOp3D(int max_bins[], typename KeyViewType::const_value_type min[],
                               typename KeyViewType::const_value_type max[] )
  {
    max_bins_[0] = max_bins[0]+1;
    max_bins_[1] = max_bins[1]+1;
    max_bins_[2] = max_bins[2]+1;
    mul_[0] = 1.0*max_bins[0]/(max[0]-min[0]);
    mul_[1] = 1.0*max_bins[1]/(max[1]-min[1]);
    mul_[2] = 1.0*max_bins[2]/(max[2]-min[2]);
    range_[0] = max[0]-min[0];
    range_[1] = max[1]-min[1];
    range_[2] = max[2]-min[2];
    min_[0] = min[0];
    min_[1] = min[1];
    min_[2] = min[2];
  }

  template<class ViewType>
  KOKKOS_INLINE_FUNCTION
  int bin(ViewType& keys, const int& i) const {
    return int( (((int(mul_[0]*(keys(i,0)-min_[0]))*max_bins_[1]) +
                   int(mul_[1]*(keys(i,1)-min_[1])))*max_bins_[2]) +
                   int(mul_[2]*(keys(i,2)-min_[2])));
  }

  KOKKOS_INLINE_FUNCTION
  int max_bins() const {
    return max_bins_[0]*max_bins_[1]*max_bins_[2];
  }

  template<class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION
  bool operator()(ViewType& keys, iType1& i1 , iType2& i2) const {
    if (keys(i1,0)>keys(i2,0)) return true;
    else if (keys(i1,0)==keys(i2,0)) {
      if (keys(i1,1)>keys(i2,1)) return true;
      else if (keys(i1,1)==keys(i2,2)) {
        if (keys(i1,2)>keys(i2,2)) return true;
      }
    }
    return false;
  }
};

template<typename Scalar>
struct min_max {
  Scalar min;
  Scalar max;
  bool init;

  KOKKOS_INLINE_FUNCTION
  min_max() {
    min = 0;
    max = 0;
    init = 0;
  }

  KOKKOS_INLINE_FUNCTION
  min_max (const min_max& val) {
    min = val.min;
    max = val.max;
    init = val.init;
  }

  KOKKOS_INLINE_FUNCTION
  min_max operator = (const min_max& val) {
    min = val.min;
    max = val.max;
    init = val.init;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+= (const Scalar& val) {
    if(init) {
      min = min<val?min:val;
      max = max>val?max:val;
    } else {
      min = val;
      max = val;
      init = 1;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator+= (const min_max& val) {
    if(init && val.init) {
      min = min<val.min?min:val.min;
      max = max>val.max?max:val.max;
    } else {
      if(val.init) {
        min = val.min;
        max = val.max;
        init = 1;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator+= (volatile const Scalar& val) volatile {
    if(init) {
      min = min<val?min:val;
      max = max>val?max:val;
    } else {
      min = val;
      max = val;
      init = 1;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator+= (volatile const min_max& val) volatile {
    if(init && val.init) {
      min = min<val.min?min:val.min;
      max = max>val.max?max:val.max;
    } else {
      if(val.init) {
        min = val.min;
        max = val.max;
        init = 1;
      }
    }
  }
};


template<class ViewType>
struct min_max_functor {
  typedef typename ViewType::execution_space execution_space;
  ViewType view;
  typedef min_max<typename ViewType::non_const_value_type> value_type;
  min_max_functor (const ViewType view_):view(view_) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t& i, value_type& val) const {
    val += view(i);
  }
};

template<class ViewType>
bool try_std_sort(ViewType view) {
  bool possible = true;
  size_t stride[8];
  view.stride(stride);
  possible  = possible && Impl::is_same<typename ViewType::memory_space, HostSpace>::value;
  possible  = possible && (ViewType::Rank == 1);
  possible  = possible && (stride[0] == 1);
  if(possible)  {
   std::sort(view.ptr_on_device(),view.ptr_on_device()+view.dimension_0());
  }
  return possible;
}

}

template<class ViewType>
void sort(ViewType view, bool always_use_kokkos_sort = false) {
  if(!always_use_kokkos_sort) {
    if(SortImpl::try_std_sort(view)) return;
  }

  typedef SortImpl::DefaultBinOp1D<ViewType> CompType;
  SortImpl::min_max<typename ViewType::non_const_value_type> val;
  parallel_reduce(view.dimension_0(),SortImpl::min_max_functor<ViewType>(view),val);
  BinSort<ViewType, CompType> bin_sort(view,CompType(view.dimension_0()/2,val.min,val.max),true);
  bin_sort.create_permute_vector();
  bin_sort.sort(view);
}

/*template<class ViewType, class Comparator>
void sort(ViewType view, Comparator comp, bool always_use_kokkos_sort = false) {

}*/

}

#endif
