/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef EXAMPLE_SORT_ARRAY
#define EXAMPLE_SORT_ARRAY

#include <cstdlib>
#include <algorithm>

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Timer.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Example {

template< class Device >
struct SortView {

  template< typename ValueType >
  SortView( const Kokkos::View<ValueType*,Device> v , int begin , int end )
    {
      std::sort( v.ptr_on_device() + begin , v.ptr_on_device() + end );
    }
};

}

#if defined(KOKKOS_ENABLE_CUDA)

#include <thrust/device_ptr.h>
#include <thrust/sort.h>

namespace Example {

template<>
struct SortView< Kokkos::Cuda > {
  template< typename ValueType >
  SortView( const Kokkos::View<ValueType*,Kokkos::Cuda> v , int begin , int end )
    {
      thrust::sort( thrust::device_ptr<ValueType>( v.ptr_on_device() + begin )
                  , thrust::device_ptr<ValueType>( v.ptr_on_device() + end ) );
    }
};

}

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Example {

template< class Device >
void sort_array( const size_t array_length /* length of spans of array to sort */
               , const size_t total_length /* total length of array */
               , const int print = 1 )
{
  typedef Device execution_space ;
  typedef Kokkos::View<int*,Device>  device_array_type ;

#if defined( KOKKOS_ENABLE_CUDA )

  typedef typename
    Kokkos::Impl::if_c< std::is_same< Device , Kokkos::Cuda >::value
                      , Kokkos::View<int*,Kokkos::Cuda::array_layout,Kokkos::CudaHostPinnedSpace>
                      , typename device_array_type::HostMirror
                      >::type  host_array_type ;

#else

  typedef typename device_array_type::HostMirror  host_array_type ;

#endif

  Kokkos::Timer timer;

  const device_array_type  work_array("work_array" , array_length );
  const host_array_type    host_array("host_array" , total_length );

  std::cout << "sort_array length( " << total_length << " )"
            << " in chunks( " << array_length << " )"
            << std::endl ;

  double sec = timer.seconds();
  std::cout << "declaring Views took "
            << sec << " seconds" << std::endl;
  timer.reset();

  for ( size_t i = 0 ; i < total_length ; ++i ) {
    host_array(i) = ( lrand48() * total_length ) >> 31 ;
  }

  sec = timer.seconds();
  std::cout << "initializing " << total_length << " elements on host took "
            << sec << " seconds" << std::endl;
  timer.reset();

  double sec_copy_in  = 0 ;
  double sec_sort     = 0 ;
  double sec_copy_out = 0 ;
  double sec_error    = 0 ;
  size_t error_count  = 0 ;

  for ( size_t begin = 0 ; begin < total_length ; begin += array_length ) {

    const size_t end = begin + array_length < total_length
                     ? begin + array_length : total_length ;

    const std::pair<size_t,size_t> host_range(begin,end);

    const host_array_type host_subarray = Kokkos::subview( host_array , host_range );

    timer.reset();

    Kokkos::deep_copy( work_array , host_subarray );

    sec_copy_in += timer.seconds(); timer.reset();

    SortView< execution_space >( work_array , 0 , end - begin );

    sec_sort += timer.seconds(); timer.reset();

    Kokkos::deep_copy( host_subarray , work_array );

    sec_copy_out += timer.seconds(); timer.reset();

    for ( size_t i = begin + 1 ; i < end ; ++i ) {
      if ( host_array(i) < host_array(i-1) ) ++error_count ;
    }

    sec_error += timer.seconds(); timer.reset();
  }

  std::cout << "copy to   device " << sec_copy_in  << " seconds" << std::endl
            << "sort on   device " << sec_sort     << " seconds" << std::endl
            << "copy from device " << sec_copy_out << " seconds" << std::endl
            << "errors " << error_count << " took " << sec_error << " seconds" << std::endl
            ;
}

} // namespace Example

//----------------------------------------------------------------------------

#endif /* #ifndef EXAMPLE_SORT_ARRAY */

