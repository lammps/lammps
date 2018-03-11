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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CORE_CONCEPTS_HPP
#define KOKKOS_CORE_CONCEPTS_HPP

#include <type_traits>

// Needed for 'is_space<S>::host_mirror_space
#include <Kokkos_Core_fwd.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//Schedules for Execution Policies
struct Static {};
struct Dynamic {};

//Schedule Wrapper Type
template<class T>
struct Schedule
{
  static_assert(  std::is_same<T,Static>::value
               || std::is_same<T,Dynamic>::value
               , "Kokkos: Invalid Schedule<> type."
               );
  using schedule_type = Schedule ;
  using type = T;
};

//Specify Iteration Index Type
template<typename T>
struct IndexType
{
  static_assert(std::is_integral<T>::value,"Kokkos: Invalid IndexType<>.");
  using index_type = IndexType ;
  using type = T;
};

/**\brief Specify Launch Bounds for CUDA execution.
 *
 *  If no launch bounds specified then do not set launch bounds.
 */
template< unsigned int maxT = 0 /* Max threads per block */
        , unsigned int minB = 0 /* Min blocks per SM */
        >
struct LaunchBounds
{
  using launch_bounds = LaunchBounds;
  using type = LaunchBounds<maxT,minB>;
  static unsigned int constexpr maxTperB {maxT};
  static unsigned int constexpr minBperSM {minB};
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

#define KOKKOS_IMPL_IS_CONCEPT( CONCEPT ) \
  template< typename T > struct is_ ## CONCEPT { \
  private: \
    template< typename , typename = std::true_type > struct have : std::false_type {}; \
    template< typename U > struct have<U,typename std::is_same<U,typename U:: CONCEPT >::type> : std::true_type {}; \
  public: \
    enum { value = is_ ## CONCEPT::template have<T>::value }; \
  };

// Public concept:

KOKKOS_IMPL_IS_CONCEPT( memory_space )
KOKKOS_IMPL_IS_CONCEPT( memory_traits )
KOKKOS_IMPL_IS_CONCEPT( execution_space )
KOKKOS_IMPL_IS_CONCEPT( execution_policy )
KOKKOS_IMPL_IS_CONCEPT( array_layout )
KOKKOS_IMPL_IS_CONCEPT( reducer )

namespace Impl {

// For backward compatibility:

using Kokkos::is_memory_space ;
using Kokkos::is_memory_traits ;
using Kokkos::is_execution_space ;
using Kokkos::is_execution_policy ;
using Kokkos::is_array_layout ;

// Implementation concept:

KOKKOS_IMPL_IS_CONCEPT( iteration_pattern )
KOKKOS_IMPL_IS_CONCEPT( schedule_type )
KOKKOS_IMPL_IS_CONCEPT( index_type )
KOKKOS_IMPL_IS_CONCEPT( launch_bounds )

}

#undef KOKKOS_IMPL_IS_CONCEPT

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

template< class ExecutionSpace , class MemorySpace >
struct Device {
  static_assert( Kokkos::is_execution_space<ExecutionSpace>::value
               , "Execution space is not valid" );
  static_assert( Kokkos::is_memory_space<MemorySpace>::value
               , "Memory space is not valid" );
  typedef ExecutionSpace                        execution_space;
  typedef MemorySpace                           memory_space;
  typedef Device<execution_space,memory_space>  device_type;
};


template< typename T >
struct is_space {
private:

  template< typename , typename = void >
  struct exe : std::false_type { typedef void space ; };

  template< typename , typename = void >
  struct mem : std::false_type { typedef void space ; };

  template< typename , typename = void >
  struct dev : std::false_type { typedef void space ; };

  template< typename U >
  struct exe<U,typename std::conditional<true,void,typename U::execution_space>::type>
    : std::is_same<U,typename U::execution_space>::type
    { typedef typename U::execution_space space ; };

  template< typename U >
  struct mem<U,typename std::conditional<true,void,typename U::memory_space>::type>
    : std::is_same<U,typename U::memory_space>::type
    { typedef typename U::memory_space space ; };

  template< typename U >
  struct dev<U,typename std::conditional<true,void,typename U::device_type>::type>
    : std::is_same<U,typename U::device_type>::type
    { typedef typename U::device_type space ; };

  typedef typename is_space::template exe<T> is_exe ;
  typedef typename is_space::template mem<T> is_mem ;
  typedef typename is_space::template dev<T> is_dev ;

public:

  enum { value = is_exe::value || is_mem::value || is_dev::value };

  typedef typename is_exe::space execution_space ;
  typedef typename is_mem::space memory_space ;

  // For backward compatibility, deprecated in favor of
  // Kokkos::Impl::HostMirror<S>::host_mirror_space

  typedef typename std::conditional
    < std::is_same< memory_space , Kokkos::HostSpace >::value
#if defined( KOKKOS_ENABLE_CUDA )
      || std::is_same< memory_space , Kokkos::CudaUVMSpace >::value
      || std::is_same< memory_space , Kokkos::CudaHostPinnedSpace >::value
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
    , memory_space
    , Kokkos::HostSpace
    >::type  host_memory_space ;

#if defined( KOKKOS_ENABLE_CUDA )
  typedef typename std::conditional
    < std::is_same< execution_space , Kokkos::Cuda >::value
    , Kokkos::DefaultHostExecutionSpace , execution_space
    >::type  host_execution_space ;
#else
  #if defined( KOKKOS_ENABLE_OPENMPTARGET )
    typedef typename std::conditional
      < std::is_same< execution_space , Kokkos::Experimental::OpenMPTarget >::value
      , Kokkos::DefaultHostExecutionSpace , execution_space
      >::type  host_execution_space ;
  #else
    typedef execution_space  host_execution_space ;
  #endif
#endif

  typedef typename std::conditional
    < std::is_same< execution_space , host_execution_space >::value &&
      std::is_same< memory_space ,    host_memory_space    >::value
    , T , Kokkos::Device< host_execution_space , host_memory_space >
    >::type  host_mirror_space ;
};

// For backward compatiblity

namespace Impl {

using Kokkos::is_space ;

}

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Access relationship between DstMemorySpace and SrcMemorySpace
 *
 *  The default case can assume accessibility for the same space.
 *  Specializations must be defined for different memory spaces.
 */
template< typename DstMemorySpace , typename SrcMemorySpace >
struct MemorySpaceAccess {

  static_assert( Kokkos::is_memory_space< DstMemorySpace >::value &&
                 Kokkos::is_memory_space< SrcMemorySpace >::value
               , "template arguments must be memory spaces" );

  /**\brief  Can a View (or pointer) to memory in SrcMemorySpace
   *         be assigned to a View (or pointer) to memory marked DstMemorySpace.
   *
   *  1. DstMemorySpace::execution_space == SrcMemorySpace::execution_space
   *  2. All execution spaces that can access DstMemorySpace can also access
   *     SrcMemorySpace.
   */
  enum { assignable = std::is_same<DstMemorySpace,SrcMemorySpace>::value };

  /**\brief  For all DstExecSpace::memory_space == DstMemorySpace
   *         DstExecSpace can access SrcMemorySpace.
   */
  enum { accessible = assignable };

  /**\brief  Does a DeepCopy capability exist
   *         to DstMemorySpace from SrcMemorySpace
   */
  enum { deepcopy = assignable };
};

}} // namespace Kokkos::Impl

namespace Kokkos {

/**\brief  Can AccessSpace access MemorySpace ?
 *
 *   Requires:
 *     Kokkos::is_space< AccessSpace >::value
 *     Kokkos::is_memory_space< MemorySpace >::value
 *
 *   Can AccessSpace::execution_space access MemorySpace ?
 *     enum : bool { accessible };
 *
 *   Is View<AccessSpace::memory_space> assignable from View<MemorySpace> ?
 *     enum : bool { assignable };
 *
 *   If ! accessible then through which intercessory memory space
 *   should a be used to deep copy memory for
 *     AccessSpace::execution_space
 *   to get access.
 *   When AccessSpace::memory_space == Kokkos::HostSpace
 *   then space is the View host mirror space.
 */
template< typename AccessSpace , typename MemorySpace >
struct SpaceAccessibility {
private:

  static_assert( Kokkos::is_space< AccessSpace >::value
               , "template argument #1 must be a Kokkos space" );

  static_assert( Kokkos::is_memory_space< MemorySpace >::value
               , "template argument #2 must be a Kokkos memory space" );

  // The input AccessSpace may be a Device<ExecSpace,MemSpace>
  // verify that it is a valid combination of spaces.
  static_assert( Kokkos::Impl::MemorySpaceAccess
                   < typename AccessSpace::execution_space::memory_space
                   , typename AccessSpace::memory_space
                   >::accessible
               , "template argument #1 is an invalid space" );

  typedef Kokkos::Impl::MemorySpaceAccess
    < typename AccessSpace::execution_space::memory_space , MemorySpace >
      exe_access ;

  typedef Kokkos::Impl::MemorySpaceAccess
    < typename AccessSpace::memory_space , MemorySpace >
      mem_access ;

public:

  /**\brief  Can AccessSpace::execution_space access MemorySpace ?
   *
   *  Default based upon memory space accessibility.
   *  Specialization required for other relationships.
   */
  enum { accessible = exe_access::accessible };

  /**\brief  Can assign to AccessSpace from MemorySpace ?
   *
   *  Default based upon memory space accessibility.
   *  Specialization required for other relationships.
   */
  enum { assignable =
    is_memory_space< AccessSpace >::value && mem_access::assignable };

  /**\brief  Can deep copy to AccessSpace::memory_Space from MemorySpace ?  */
  enum { deepcopy = mem_access::deepcopy };

  // What intercessory space for AccessSpace::execution_space
  // to be able to access MemorySpace?
  // If same memory space or not accessible use the AccessSpace
  // else construct a device with execution space and memory space.
  typedef typename std::conditional
    < std::is_same<typename AccessSpace::memory_space,MemorySpace>::value ||
      ! exe_access::accessible
    , AccessSpace
    , Kokkos::Device< typename AccessSpace::execution_space , MemorySpace >
    >::type  space ;
};

} // namespace Kokkos

namespace Kokkos {
namespace Impl {

using Kokkos::SpaceAccessibility ; // For backward compatibility

}} // namespace Kokkos::Impl

//----------------------------------------------------------------------------

#endif // KOKKOS_CORE_CONCEPTS_HPP

