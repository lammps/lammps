/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_SEGMENTED_VIEW_HPP_
#define KOKKOS_SEGMENTED_VIEW_HPP_

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <cstdio>

namespace Kokkos {

namespace Impl {

template<class DataType, class Arg1Type, class Arg2Type, class Arg3Type>
struct delete_segmented_view;

template<class MemorySpace>
inline
void DeviceSetAllocatableMemorySize(size_t) {}

#if defined( KOKKOS_HAVE_CUDA )

template<>
inline
void DeviceSetAllocatableMemorySize<Kokkos::CudaSpace>(size_t size) {
#ifdef __CUDACC__
  size_t size_limit;
  cudaDeviceGetLimit(&size_limit,cudaLimitMallocHeapSize);
  if(size_limit<size)
    cudaDeviceSetLimit(cudaLimitMallocHeapSize,2*size);
  cudaDeviceGetLimit(&size_limit,cudaLimitMallocHeapSize);
#endif
}

template<>
inline
void DeviceSetAllocatableMemorySize<Kokkos::CudaUVMSpace>(size_t size) {
#ifdef __CUDACC__
  size_t size_limit;
  cudaDeviceGetLimit(&size_limit,cudaLimitMallocHeapSize);
  if(size_limit<size)
    cudaDeviceSetLimit(cudaLimitMallocHeapSize,2*size);
  cudaDeviceGetLimit(&size_limit,cudaLimitMallocHeapSize);
#endif
}

#endif /* #if defined( KOKKOS_HAVE_CUDA ) */

}

template< class DataType ,
          class Arg1Type = void ,
          class Arg2Type = void ,
          class Arg3Type = void>
class SegmentedView : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:
  //! \name Typedefs for device types and various Kokkos::View specializations.
  //@{
  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

  //! The type of a Kokkos::View on the device.
  typedef View< typename traits::data_type ,
                typename traits::array_layout ,
                typename traits::memory_space ,
                Kokkos::MemoryUnmanaged > t_dev ;


private:
  Kokkos::View<t_dev*,typename traits::memory_space> segments_;

  Kokkos::View<int,typename traits::memory_space> realloc_lock;
  Kokkos::View<int,typename traits::memory_space> nsegments_;

  size_t segment_length_;
  size_t segment_length_m1_;
  int max_segments_;

  int segment_length_log2;

  // Dimensions, cardinality, capacity, and offset computation for
  // multidimensional array view of contiguous memory.
  // Inherits from Impl::Shape
  typedef Impl::ViewOffset< typename traits::shape_type
                          , typename traits::array_layout
                          > offset_map_type ;

  offset_map_type               m_offset_map ;

  typedef View< typename traits::array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::memory_space ,
                typename traits::memory_traits > array_type ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::memory_space ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::memory_space ,
                typename traits::memory_traits > non_const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                HostSpace ,
                void > HostMirror ;

  template< bool Accessible >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if< Accessible , typename traits::size_type >::type
  dimension_0_intern() const { return nsegments_() * segment_length_ ; }

  template< bool Accessible >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if< ! Accessible , typename traits::size_type >::type
  dimension_0_intern() const
  {
    // In Host space
    int n = 0 ;
#if ! defined( __CUDA_ARCH__ )
    Impl::DeepCopy< HostSpace , typename traits::memory_space >( & n , nsegments_.ptr_on_device() , sizeof(int) );
#endif

    return n * segment_length_ ;
  }

public:

  enum { Rank = traits::rank };

  KOKKOS_INLINE_FUNCTION offset_map_type shape() const { return m_offset_map ; }

  /* \brief return (current) size of dimension 0 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const {
    enum { Accessible = Impl::VerifyExecutionCanAccessMemorySpace<
             Impl::ActiveExecutionMemorySpace, typename traits::memory_space >::value };
    int n = SegmentedView::dimension_0_intern< Accessible >();
    return n ;
  }

  /* \brief return size of dimension 1 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_offset_map.N1 ; }
  /* \brief return size of dimension 2 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_offset_map.N2 ; }
  /* \brief return size of dimension 3 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_offset_map.N3 ; }
  /* \brief return size of dimension 4 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_offset_map.N4 ; }
  /* \brief return size of dimension 5 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_offset_map.N5 ; }
  /* \brief return size of dimension 6 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_offset_map.N6 ; }
  /* \brief return size of dimension 7 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_offset_map.N7 ; }

  /* \brief return size of dimension 2 */
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const {
    return dimension_0() *
        m_offset_map.N1 * m_offset_map.N2 * m_offset_map.N3 * m_offset_map.N4 *
        m_offset_map.N5 * m_offset_map.N6 * m_offset_map.N7 ;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const {
    if(i==0)
      return dimension_0();
    else
      return Impl::dimension( m_offset_map , i );
  }

  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() {
    return segments_.dimension_0() *
        m_offset_map.N1 * m_offset_map.N2 * m_offset_map.N3 * m_offset_map.N4 *
        m_offset_map.N5 * m_offset_map.N6 * m_offset_map.N7;
  }

  KOKKOS_INLINE_FUNCTION
  typename traits::size_type get_num_segments() {
    enum { Accessible = Impl::VerifyExecutionCanAccessMemorySpace<
             Impl::ActiveExecutionMemorySpace, typename traits::memory_space >::value };
    int n = SegmentedView::dimension_0_intern< Accessible >();
    return n/segment_length_ ;
  }

  KOKKOS_INLINE_FUNCTION
  typename traits::size_type get_max_segments() {
    return max_segments_;
  }

  /// \brief Constructor that allocates View objects with an initial length of 0.
  ///
  /// This constructor works mostly like the analogous constructor of View.
  /// The first argument is a string label, which is entirely for your
  /// benefit.  (Different SegmentedView objects may have the same label if
  /// you like.)  The second argument 'view_length' is the size of the segments.
  /// This number must be a power of two. The third argument n0 is the maximum
  /// value for the first dimension of the segmented view. The maximal allocatable
  /// number of Segments is thus: (n0+view_length-1)/view_length.
  /// The arguments that follow are the other dimensions of the (1-7) of the
  /// View objects.  For example, for a View with 3 runtime dimensions,
  /// the first 4 integer arguments will be nonzero:
  /// SegmentedView("Name",32768,10000000,8,4). This allocates a SegmentedView
  /// with a maximum of 306 segments of dimension (32768,8,4). The logical size of
  /// the segmented view is (n,8,4) with n between 0 and 10000000.
  /// You may omit the integer arguments that follow.
  template< class LabelType >
  SegmentedView(const LabelType & label ,
      const size_t view_length ,
      const size_t n0 ,
      const size_t n1 = 0 ,
      const size_t n2 = 0 ,
      const size_t n3 = 0 ,
      const size_t n4 = 0 ,
      const size_t n5 = 0 ,
      const size_t n6 = 0 ,
      const size_t n7 = 0
      ): segment_length_(view_length),segment_length_m1_(view_length-1)
  {
    segment_length_log2 = -1;
    size_t l = segment_length_;
    while(l>0) {
      l>>=1;
      segment_length_log2++;
    }
    l = 1<<segment_length_log2;
    if(l!=segment_length_)
      Impl::throw_runtime_exception("Kokkos::SegmentedView requires a 'power of 2' segment length");

    max_segments_ = (n0+segment_length_m1_)/segment_length_;

    Impl::DeviceSetAllocatableMemorySize<typename traits::memory_space>(segment_length_*max_segments_*sizeof(typename traits::value_type));

    segments_ = Kokkos::View<t_dev*,typename traits::execution_space>(label , max_segments_);
    realloc_lock = Kokkos::View<int,typename traits::execution_space>("Lock");
    nsegments_ = Kokkos::View<int,typename traits::execution_space>("nviews");
    m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7, n0*n1*n2*n3*n4*n5*n6*n7 );

  }

  KOKKOS_INLINE_FUNCTION
  SegmentedView(const SegmentedView& src):
    segments_(src.segments_),
    realloc_lock (src.realloc_lock),
    nsegments_ (src.nsegments_),
    segment_length_(src.segment_length_),
    segment_length_m1_(src.segment_length_m1_),
    max_segments_ (src.max_segments_),
    segment_length_log2(src.segment_length_log2),
    m_offset_map (src.m_offset_map)
  {}

  KOKKOS_INLINE_FUNCTION
  SegmentedView& operator= (const SegmentedView& src) {
    segments_ = src.segments_;
    realloc_lock = src.realloc_lock;
    nsegments_ = src.nsegments_;
    segment_length_= src.segment_length_;
    segment_length_m1_= src.segment_length_m1_;
    max_segments_ = src.max_segments_;
    segment_length_log2= src.segment_length_log2;
    m_offset_map = src.m_offset_map;
    return *this;
  }

  ~SegmentedView() {
    if (traits::execution_space::in_parallel()) return;
    int count = traits::memory_space::count(segments_.ptr_on_device());
    if(count == 1) {
      Kokkos::fence();
      typename Kokkos::View<int,typename traits::execution_space>::HostMirror h_nviews("h_nviews");
      Kokkos::deep_copy(h_nviews,nsegments_);
      Kokkos::parallel_for(h_nviews(),Impl::delete_segmented_view<DataType , Arg1Type , Arg2Type, Arg3Type>(*this));
    }
  }

  KOKKOS_INLINE_FUNCTION
  t_dev get_segment(const int& i) const {
    return segments_[i];
  }

  template< class MemberType>
  KOKKOS_INLINE_FUNCTION
  void grow (MemberType& team_member, const size_t& growSize) const {
    if (growSize>max_segments_*segment_length_) {
      printf ("Exceeding maxSize: %lu %lu\n", growSize, max_segments_*segment_length_);
      return;
    }
    if(team_member.team_rank()==0) {
      bool too_small = growSize > segment_length_ * nsegments_();
      while(too_small && Kokkos::atomic_compare_exchange(&realloc_lock(),0,1) ) {
        too_small = growSize > segment_length_ * nsegments_();
      }
      if(too_small) {
        while(too_small) {
          const size_t alloc_size = segment_length_*m_offset_map.N1*m_offset_map.N2*m_offset_map.N3*
                              m_offset_map.N4*m_offset_map.N5*m_offset_map.N6*m_offset_map.N7;
          typename traits::non_const_value_type* const ptr = new typename traits::non_const_value_type[alloc_size];

          segments_(nsegments_()) =
            t_dev(ptr,segment_length_,m_offset_map.N1,m_offset_map.N2,m_offset_map.N3,m_offset_map.N4,m_offset_map.N5,m_offset_map.N6,m_offset_map.N7);
          nsegments_()++;
          too_small = growSize > segment_length_ * nsegments_();
        }
        realloc_lock() = 0;
      }
    }
    team_member.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void grow_non_thread_safe (const size_t& growSize) const {
    if (growSize>max_segments_*segment_length_) {
      printf ("Exceeding maxSize: %lu %lu\n", growSize, max_segments_*segment_length_);
      return;
    }
    bool too_small = growSize > segment_length_ * nsegments_();
    if(too_small) {
      while(too_small) {
        const size_t alloc_size = segment_length_*m_offset_map.N1*m_offset_map.N2*m_offset_map.N3*
                            m_offset_map.N4*m_offset_map.N5*m_offset_map.N6*m_offset_map.N7;
        typename traits::non_const_value_type* const ptr =
          new typename traits::non_const_value_type[alloc_size];

        segments_(nsegments_()) =
          t_dev (ptr, segment_length_, m_offset_map.N1, m_offset_map.N2,
                 m_offset_map.N3, m_offset_map.N4, m_offset_map.N5,
                 m_offset_map.N6, m_offset_map.N7);
        nsegments_()++;
        too_small = growSize > segment_length_ * nsegments_();
      }
    }
  }

  template< typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 1, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_));
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 2,
               iType0 , iType1>::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1);
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 3,
               iType0 , iType1 , iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2);
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 4,
               iType0 , iType1 , iType2 , iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2,i3);
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 5,
               iType0 , iType1 , iType2 , iType3 , iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2,i3,i4);
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 6,
               iType0 , iType1 , iType2 , iType3 , iType4 , iType5>::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2,i3,i4,i5);
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 7,
               iType0 , iType1 , iType2 , iType3 , iType4 , iType5 , iType6>::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2,i3,i4,i5,i6);
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type & , traits, typename traits::array_layout, 8,
               iType0 , iType1 , iType2 , iType3 , iType4 , iType5 , iType6 , iType7>::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      return segments_[i0>>segment_length_log2](i0&(segment_length_m1_),i1,i2,i3,i4,i5,i6,i7);
    }
};

namespace Impl {
template<class DataType, class Arg1Type, class Arg2Type, class Arg3Type>
struct delete_segmented_view {
  typedef SegmentedView<DataType , Arg1Type , Arg2Type, Arg3Type> view_type;
  typedef typename view_type::execution_space device_type;

  view_type view_;
  delete_segmented_view(view_type view):view_(view) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    delete [] view_.get_segment(i).ptr_on_device();
  }
};

}
}

#endif
