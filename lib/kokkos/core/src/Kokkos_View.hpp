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

#ifndef KOKKOS_VIEW_HPP
#define KOKKOS_VIEW_HPP

#include <type_traits>
#include <string>
#include <algorithm>
#include <initializer_list>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_ExecPolicy.hpp>

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DataType >
struct ViewArrayAnalysis ;

template< class DataType , class ArrayLayout
        , typename ValueType =
          typename ViewArrayAnalysis< DataType >::non_const_value_type
        >
struct ViewDataAnalysis ;

template< class , class ... >
class ViewMapping { public: enum { is_assignable = false }; };

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \class ViewTraits
 *  \brief Traits class for accessing attributes of a View.
 *
 * This is an implementation detail of View.  It is only of interest
 * to developers implementing a new specialization of View.
 *
 * Template argument options:
 *   - View< DataType >
 *   - View< DataType , Space >
 *   - View< DataType , Space , MemoryTraits >
 *   - View< DataType , ArrayLayout >
 *   - View< DataType , ArrayLayout , Space >
 *   - View< DataType , ArrayLayout , MemoryTraits >
 *   - View< DataType , ArrayLayout , Space , MemoryTraits >
 *   - View< DataType , MemoryTraits >
 */

template< class DataType , class ... Properties >
struct ViewTraits ;

template<>
struct ViewTraits< void >
{
  typedef void  execution_space ;
  typedef void  memory_space ;
  typedef void  HostMirrorSpace ;
  typedef void  array_layout ;
  typedef void  memory_traits ;
};

template< class ... Prop >
struct ViewTraits< void , void , Prop ... >
{
  // Ignore an extraneous 'void'
  typedef typename ViewTraits<void,Prop...>::execution_space  execution_space ;
  typedef typename ViewTraits<void,Prop...>::memory_space     memory_space ;
  typedef typename ViewTraits<void,Prop...>::HostMirrorSpace  HostMirrorSpace ;
  typedef typename ViewTraits<void,Prop...>::array_layout     array_layout ;
  typedef typename ViewTraits<void,Prop...>::memory_traits    memory_traits ;
};

template< class ArrayLayout , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_array_layout<ArrayLayout>::value >::type , ArrayLayout , Prop ... >
{
  // Specify layout, keep subsequent space and memory traits arguments

  typedef typename ViewTraits<void,Prop...>::execution_space  execution_space ;
  typedef typename ViewTraits<void,Prop...>::memory_space     memory_space ;
  typedef typename ViewTraits<void,Prop...>::HostMirrorSpace  HostMirrorSpace ;
  typedef          ArrayLayout                                array_layout ;
  typedef typename ViewTraits<void,Prop...>::memory_traits    memory_traits ;
};

template< class Space , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_space<Space>::value >::type , Space , Prop ... >
{
  // Specify Space, memory traits should be the only subsequent argument.

  static_assert( std::is_same< typename ViewTraits<void,Prop...>::execution_space , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::memory_space    , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::HostMirrorSpace , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::array_layout    , void >::value
               , "Only one View Execution or Memory Space template argument" );

  typedef typename Space::execution_space                   execution_space ;
  typedef typename Space::memory_space                      memory_space ;
  typedef typename Kokkos::Impl::HostMirror< Space >::Space HostMirrorSpace ;
  typedef typename execution_space::array_layout            array_layout ;
  typedef typename ViewTraits<void,Prop...>::memory_traits  memory_traits ;
};

template< class MemoryTraits , class ... Prop >
struct ViewTraits< typename std::enable_if< Kokkos::Impl::is_memory_traits<MemoryTraits>::value >::type , MemoryTraits , Prop ... >
{
  // Specify memory trait, should not be any subsequent arguments

  static_assert( std::is_same< typename ViewTraits<void,Prop...>::execution_space , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::memory_space    , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::array_layout    , void >::value &&
                 std::is_same< typename ViewTraits<void,Prop...>::memory_traits   , void >::value
               , "MemoryTrait is the final optional template argument for a View" );

  typedef void          execution_space ;
  typedef void          memory_space ;
  typedef void          HostMirrorSpace ;
  typedef void          array_layout ;
  typedef MemoryTraits  memory_traits ;
};


template< class DataType , class ... Properties >
struct ViewTraits {
private:

  // Unpack the properties arguments
  typedef ViewTraits< void , Properties ... >  prop ;

  typedef typename
    std::conditional< ! std::is_same< typename prop::execution_space , void >::value
                    , typename prop::execution_space
                    , Kokkos::DefaultExecutionSpace
                    >::type
      ExecutionSpace ;

  typedef typename
    std::conditional< ! std::is_same< typename prop::memory_space , void >::value
                    , typename prop::memory_space
                    , typename ExecutionSpace::memory_space
                    >::type
      MemorySpace ;

  typedef typename
    std::conditional< ! std::is_same< typename prop::array_layout , void >::value
                    , typename prop::array_layout
                    , typename ExecutionSpace::array_layout
                    >::type
      ArrayLayout ;

  typedef typename
    std::conditional
      < ! std::is_same< typename prop::HostMirrorSpace , void >::value
      , typename prop::HostMirrorSpace
      , typename Kokkos::Impl::HostMirror< ExecutionSpace >::Space
      >::type
      HostMirrorSpace ;

  typedef typename
    std::conditional< ! std::is_same< typename prop::memory_traits , void >::value
                    , typename prop::memory_traits
                    , typename Kokkos::MemoryManaged
                    >::type
      MemoryTraits ;

  // Analyze data type's properties,
  // May be specialized based upon the layout and value type
  typedef Kokkos::Impl::ViewDataAnalysis< DataType , ArrayLayout > data_analysis ;

public:

  //------------------------------------
  // Data type traits:

  typedef typename data_analysis::type            data_type ;
  typedef typename data_analysis::const_type      const_data_type ;
  typedef typename data_analysis::non_const_type  non_const_data_type ;

  //------------------------------------
  // Compatible array of trivial type traits:

  typedef typename data_analysis::scalar_array_type            scalar_array_type ;
  typedef typename data_analysis::const_scalar_array_type      const_scalar_array_type ;
  typedef typename data_analysis::non_const_scalar_array_type  non_const_scalar_array_type ;

  //------------------------------------
  // Value type traits:

  typedef typename data_analysis::value_type            value_type ;
  typedef typename data_analysis::const_value_type      const_value_type ;
  typedef typename data_analysis::non_const_value_type  non_const_value_type ;

  //------------------------------------
  // Mapping traits:

  typedef ArrayLayout                         array_layout ;
  typedef typename data_analysis::dimension   dimension ;
  typedef typename data_analysis::specialize  specialize /* mapping specialization tag */ ;

  enum { rank         = dimension::rank };
  enum { rank_dynamic = dimension::rank_dynamic };

  //------------------------------------
  // Execution space, memory space, memory access traits, and host mirror space.

  typedef ExecutionSpace                              execution_space ;
  typedef MemorySpace                                 memory_space ;
  typedef Kokkos::Device<ExecutionSpace,MemorySpace>  device_type ;
  typedef MemoryTraits                                memory_traits ;
  typedef HostMirrorSpace                             host_mirror_space ;

  typedef typename MemorySpace::size_type  size_type ;

  enum { is_hostspace      = std::is_same< MemorySpace , HostSpace >::value };
  enum { is_managed        = MemoryTraits::Unmanaged    == 0 };
  enum { is_random_access  = MemoryTraits::RandomAccess == 1 };

  //------------------------------------
};

/** \class View
 *  \brief View to an array of data.
 *
 * A View represents an array of one or more dimensions.
 * For details, please refer to Kokkos' tutorial materials.
 *
 * \section Kokkos_View_TemplateParameters Template parameters
 *
 * This class has both required and optional template parameters.  The
 * \c DataType parameter must always be provided, and must always be
 * first. The parameters \c Arg1Type, \c Arg2Type, and \c Arg3Type are
 * placeholders for different template parameters.  The default value
 * of the fifth template parameter \c Specialize suffices for most use
 * cases.  When explaining the template parameters, we won't refer to
 * \c Arg1Type, \c Arg2Type, and \c Arg3Type; instead, we will refer
 * to the valid categories of template parameters, in whatever order
 * they may occur.
 *
 * Valid ways in which template arguments may be specified:
 *   - View< DataType >
 *   - View< DataType , Layout >
 *   - View< DataType , Layout , Space >
 *   - View< DataType , Layout , Space , MemoryTraits >
 *   - View< DataType , Space >
 *   - View< DataType , Space , MemoryTraits >
 *   - View< DataType , MemoryTraits >
 *
 * \tparam DataType (required) This indicates both the type of each
 *   entry of the array, and the combination of compile-time and
 *   run-time array dimension(s).  For example, <tt>double*</tt>
 *   indicates a one-dimensional array of \c double with run-time
 *   dimension, and <tt>int*[3]</tt> a two-dimensional array of \c int
 *   with run-time first dimension and compile-time second dimension
 *   (of 3).  In general, the run-time dimensions (if any) must go
 *   first, followed by zero or more compile-time dimensions.  For
 *   more examples, please refer to the tutorial materials.
 *
 * \tparam Space (required) The memory space.
 *
 * \tparam Layout (optional) The array's layout in memory.  For
 *   example, LayoutLeft indicates a column-major (Fortran style)
 *   layout, and LayoutRight a row-major (C style) layout.  If not
 *   specified, this defaults to the preferred layout for the
 *   <tt>Space</tt>.
 *
 * \tparam MemoryTraits (optional) Assertion of the user's intended
 *   access behavior.  For example, RandomAccess indicates read-only
 *   access with limited spatial locality, and Unmanaged lets users
 *   wrap externally allocated memory in a View without automatic
 *   deallocation.
 *
 * \section Kokkos_View_MT MemoryTraits discussion
 *
 * \subsection Kokkos_View_MT_Interp MemoryTraits interpretation depends on Space
 *
 * Some \c MemoryTraits options may have different interpretations for
 * different \c Space types.  For example, with the Cuda device,
 * \c RandomAccess tells Kokkos to fetch the data through the texture
 * cache, whereas the non-GPU devices have no such hardware construct.
 *
 * \subsection Kokkos_View_MT_PrefUse Preferred use of MemoryTraits
 *
 * Users should defer applying the optional \c MemoryTraits parameter
 * until the point at which they actually plan to rely on it in a
 * computational kernel.  This minimizes the number of template
 * parameters exposed in their code, which reduces the cost of
 * compilation.  Users may always assign a View without specified
 * \c MemoryTraits to a compatible View with that specification.
 * For example:
 * \code
 * // Pass in the simplest types of View possible.
 * void
 * doSomething (View<double*, Cuda> out,
 *              View<const double*, Cuda> in)
 * {
 *   // Assign the "generic" View in to a RandomAccess View in_rr.
 *   // Note that RandomAccess View objects must have const data.
 *   View<const double*, Cuda, RandomAccess> in_rr = in;
 *   // ... do something with in_rr and out ...
 * }
 * \endcode
 */
template< class DataType , class ... Properties >
class View ;

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/Kokkos_ViewMapping.hpp>
#include <impl/Kokkos_ViewArray.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

namespace {

constexpr Kokkos::Impl::ALL_t
  ALL = Kokkos::Impl::ALL_t();

constexpr Kokkos::Impl::WithoutInitializing_t
  WithoutInitializing = Kokkos::Impl::WithoutInitializing_t();

constexpr Kokkos::Impl::AllowPadding_t
  AllowPadding        = Kokkos::Impl::AllowPadding_t();

}

/** \brief  Create View allocation parameter bundle from argument list.
 *
 *  Valid argument list members are:
 *    1) label as a "string" or std::string
 *    2) memory space instance of the View::memory_space type
 *    3) execution space instance compatible with the View::memory_space
 *    4) Kokkos::WithoutInitializing to bypass initialization
 *    4) Kokkos::AllowPadding to allow allocation to pad dimensions for memory alignment
 */
template< class ... Args >
inline
Impl::ViewCtorProp< typename Impl::ViewCtorProp< void , Args >::type ... >
view_alloc( Args const & ... args )
{
  typedef
    Impl::ViewCtorProp< typename Impl::ViewCtorProp< void , Args >::type ... >
      return_type ;

  static_assert( ! return_type::has_pointer
               , "Cannot give pointer-to-memory for view allocation" );

  return return_type( args... );
}

template< class ... Args >
KOKKOS_INLINE_FUNCTION
Impl::ViewCtorProp< typename Impl::ViewCtorProp< void , Args >::type ... >
view_wrap( Args const & ... args )
{
  typedef
    Impl::ViewCtorProp< typename Impl::ViewCtorProp< void , Args >::type ... >
      return_type ;

  static_assert( ! return_type::has_memory_space &&
                 ! return_type::has_execution_space &&
                 ! return_type::has_label &&
                 return_type::has_pointer
               , "Must only give pointer-to-memory for view wrapping" );

  return return_type( args... );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class DataType , class ... Properties >
class View ;

template< class > struct is_view : public std::false_type {};

template< class D, class ... P >
struct is_view< View<D,P...> > : public std::true_type {};

template< class D, class ... P >
struct is_view< const View<D,P...> > : public std::true_type {};

template< class DataType , class ... Properties >
class View : public ViewTraits< DataType , Properties ... > {
private:

  template< class , class ... > friend class View ;
  template< class , class ... > friend class Kokkos::Impl::ViewMapping ;

public:

  typedef ViewTraits< DataType , Properties ... > traits ;

private:

  typedef Kokkos::Impl::ViewMapping< traits , void > map_type ;
  typedef Kokkos::Impl::SharedAllocationTracker      track_type ;

  track_type  m_track ;
  map_type    m_map ;

public:

  //----------------------------------------
  /** \brief  Compatible view of array of scalar types */
  typedef View< typename traits::scalar_array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits >
    array_type ;

  /** \brief  Compatible view of const data type */
  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits >
    const_type ;

  /** \brief  Compatible view of non-const data type */
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits >
    non_const_type ;

  /** \brief  Compatible HostMirror view */
  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space >
    HostMirror ;

  //----------------------------------------
  // Domain rank and extents

  enum { Rank = map_type::Rank };

 /** \brief rank() to be implemented
  */
  //KOKKOS_INLINE_FUNCTION
  //static
  //constexpr unsigned rank() { return map_type::Rank; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr
  typename std::enable_if< std::is_integral<iType>::value , size_t >::type
  extent( const iType & r ) const
    { return m_map.extent(r); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr
  typename std::enable_if< std::is_integral<iType>::value , int >::type
  extent_int( const iType & r ) const
    { return static_cast<int>(m_map.extent(r)); }

  KOKKOS_INLINE_FUNCTION constexpr
  typename traits::array_layout layout() const
    { return m_map.layout(); }

  //----------------------------------------
  /*  Deprecate all 'dimension' functions in favor of
   *  ISO/C++ vocabulary 'extent'.
   */

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr
  typename std::enable_if< std::is_integral<iType>::value , size_t >::type
  dimension( const iType & r ) const { return extent( r ); }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const { return m_map.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return m_map.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return m_map.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return m_map.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return m_map.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return m_map.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return m_map.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return m_map.dimension_7(); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_t size() const { return m_map.dimension_0() *
                                                                m_map.dimension_1() *
                                                                m_map.dimension_2() *
                                                                m_map.dimension_3() *
                                                                m_map.dimension_4() *
                                                                m_map.dimension_5() *
                                                                m_map.dimension_6() *
                                                                m_map.dimension_7(); }

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return m_map.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return m_map.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return m_map.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return m_map.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return m_map.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return m_map.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return m_map.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return m_map.stride_7(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const { m_map.stride(s); }

  //----------------------------------------
  // Range span is the span which contains all members.

  typedef typename map_type::reference_type  reference_type ;
  typedef typename map_type::pointer_type    pointer_type ;

  enum { reference_type_is_lvalue_reference = std::is_lvalue_reference< reference_type >::value };

  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return m_map.span(); }
  // Deprecated, use 'span()' instead
  KOKKOS_INLINE_FUNCTION constexpr size_t capacity() const { return m_map.span(); }
  KOKKOS_INLINE_FUNCTION constexpr bool   span_is_contiguous() const { return m_map.span_is_contiguous(); }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const { return m_map.data(); }

  // Deprecated, use 'span_is_contigous()' instead
  KOKKOS_INLINE_FUNCTION constexpr bool   is_contiguous() const { return m_map.span_is_contiguous(); }
  // Deprecated, use 'data()' instead
  KOKKOS_INLINE_FUNCTION constexpr pointer_type ptr_on_device() const { return m_map.data(); }

  //----------------------------------------
  // Allow specializations to query their specialized map

  KOKKOS_INLINE_FUNCTION
  const Kokkos::Impl::ViewMapping< traits , void > &
  implementation_map() const { return m_map ; }

  //----------------------------------------

private:

  enum {
    is_layout_left = std::is_same< typename traits::array_layout
                                  , Kokkos::LayoutLeft >::value ,

    is_layout_right = std::is_same< typename traits::array_layout
                                  , Kokkos::LayoutRight >::value ,

    is_layout_stride = std::is_same< typename traits::array_layout
                                   , Kokkos::LayoutStride >::value ,

    is_default_map =
      std::is_same< typename traits::specialize , void >::value &&
      ( is_layout_left || is_layout_right || is_layout_stride )
  };

  template< class Space , bool = Kokkos::Impl::MemorySpaceAccess< Space , typename traits::memory_space >::accessible > struct verify_space
    { KOKKOS_FORCEINLINE_FUNCTION static void check() {} };

  template< class Space > struct verify_space<Space,false>
    { KOKKOS_FORCEINLINE_FUNCTION static void check()
        { Kokkos::abort("Kokkos::View ERROR: attempt to access inaccessible memory space"); };
    };

#if defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )

#define KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( ARG ) \
  View::template verify_space< Kokkos::Impl::ActiveExecutionMemorySpace >::check(); \
  Kokkos::Impl::view_verify_operator_bounds< typename traits::memory_space > ARG ;

#else

#define KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( ARG ) \
  View::template verify_space< Kokkos::Impl::ActiveExecutionMemorySpace >::check();

#endif

public:

  //------------------------------
  // Rank 0 operator()

  template< class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<( Kokkos::Impl::are_integral<Args...>::value
                            && ( 0 == Rank )
                          ), reference_type >::type
  operator()( Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,args...) )
      return m_map.reference();
    }

  //------------------------------
  // Rank 1 operator()

  template< typename I0
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,Args...>::value
      && ( 1 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,args...) )
      return m_map.reference(i0);
    }

  template< typename I0
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,Args...>::value
      && ( 1 == Rank )
      && is_default_map
      && ! is_layout_stride
    ), reference_type >::type
  operator()( const I0 & i0
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,args...) )
      return m_map.m_handle[ i0 ];
    }

  template< typename I0
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,Args...>::value
      && ( 1 == Rank )
      && is_default_map
      && is_layout_stride
    ), reference_type >::type
  operator()( const I0 & i0
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,args...) )
      return m_map.m_handle[ m_map.m_offset.m_stride.S0 * i0 ];
    }

  //------------------------------
  // Rank 1 operator[]

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0>::value
      && ( 1 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator[]( const I0 & i0 ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0) )
      return m_map.reference(i0);
    }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0>::value
      && ( 1 == Rank )
      && is_default_map
      && ! is_layout_stride
    ), reference_type >::type
  operator[]( const I0 & i0 ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0) )
      return m_map.m_handle[ i0 ];
    }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0>::value
      && ( 1 == Rank )
      && is_default_map
      && is_layout_stride
    ), reference_type >::type
  operator[]( const I0 & i0 ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0) )
      return m_map.m_handle[ m_map.m_offset.m_stride.S0 * i0 ];
    }

  //------------------------------
  // Rank 2

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.reference(i0,i1);
    }

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && is_default_map
      && is_layout_left && ( traits::rank_dynamic == 0 )
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.m_handle[ i0 + m_map.m_offset.m_dim.N0 * i1 ];
    }

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && is_default_map
      && is_layout_left && ( traits::rank_dynamic != 0 )
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.m_handle[ i0 + m_map.m_offset.m_stride * i1 ];
    }

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && is_default_map
      && is_layout_right && ( traits::rank_dynamic == 0 )
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.m_handle[ i1 + m_map.m_offset.m_dim.N1 * i0 ];
    }

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && is_default_map
      && is_layout_right && ( traits::rank_dynamic != 0 )
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.m_handle[ i1 + m_map.m_offset.m_stride * i0 ];
    }

  template< typename I0 , typename I1
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,Args...>::value
      && ( 2 == Rank )
      && is_default_map
      && is_layout_stride
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,args...) )
      return m_map.m_handle[ i0 * m_map.m_offset.m_stride.S0 +
                             i1 * m_map.m_offset.m_stride.S1 ];
    }

  //------------------------------
  // Rank 3

  template< typename I0 , typename I1 , typename I2
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,Args...>::value
      && ( 3 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2) ];
    }

  template< typename I0 , typename I1 , typename I2
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,Args...>::value
      && ( 3 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,args...) )
      return m_map.reference(i0,i1,i2);
    }

  //------------------------------
  // Rank 4

  template< typename I0 , typename I1 , typename I2 , typename I3
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,Args...>::value
      && ( 4 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2,i3) ];
    }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,Args...>::value
      && ( 4 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,args...) )
      return m_map.reference(i0,i1,i2,i3);
    }

  //------------------------------
  // Rank 5

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,Args...>::value
      && ( 5 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2,i3,i4) ];
    }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,Args...>::value
      && ( 5 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,args...) )
      return m_map.reference(i0,i1,i2,i3,i4);
    }

  //------------------------------
  // Rank 6

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,Args...>::value
      && ( 6 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2,i3,i4,i5) ];
    }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,Args...>::value
      && ( 6 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,args...) )
      return m_map.reference(i0,i1,i2,i3,i4,i5);
    }

  //------------------------------
  // Rank 7

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,Args...>::value
      && ( 7 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5 , const I6 & i6
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,i6,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2,i3,i4,i5,i6) ];
    }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,Args...>::value
      && ( 7 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5 , const I6 & i6
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,i6,args...) )
      return m_map.reference(i0,i1,i2,i3,i4,i5,i6);
    }

  //------------------------------
  // Rank 8

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,I7,Args...>::value
      && ( 8 == Rank )
      && is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,i6,i7,args...) )
      return m_map.m_handle[ m_map.m_offset(i0,i1,i2,i3,i4,i5,i6,i7) ];
    }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7
          , class ... Args >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<
    ( Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,I7,Args...>::value
      && ( 8 == Rank )
      && ! is_default_map
    ), reference_type >::type
  operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
            , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7
            , Args ... args ) const
    {
      KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (m_track,m_map,i0,i1,i2,i3,i4,i5,i6,i7,args...) )
      return m_map.reference(i0,i1,i2,i3,i4,i5,i6,i7);
    }

#undef KOKKOS_IMPL_VIEW_OPERATOR_VERIFY

  //----------------------------------------
  // Standard destructor, constructors, and assignment operators

  KOKKOS_INLINE_FUNCTION
  ~View() {}

  KOKKOS_INLINE_FUNCTION
  View() : m_track(), m_map() {}

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs ) : m_track( rhs.m_track ), m_map( rhs.m_map ) {}

  KOKKOS_INLINE_FUNCTION
  View( View && rhs ) : m_track( rhs.m_track ), m_map( rhs.m_map ) {}

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { m_track = rhs.m_track ; m_map = rhs.m_map ; return *this ; }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( View && rhs ) { m_track = rhs.m_track ; m_map = rhs.m_map ; return *this ; }

  //----------------------------------------
  // Compatible view copy constructor and assignment
  // may assign unmanaged from managed.

  template< class RT , class ... RP >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RP...> & rhs )
    : m_track( rhs.m_track , traits::is_managed )
    , m_map()
    {
      typedef typename View<RT,RP...>::traits  SrcTraits ;
      typedef Kokkos::Impl::ViewMapping< traits , SrcTraits , void >  Mapping ;
      static_assert( Mapping::is_assignable , "Incompatible View copy construction" );
      Mapping::assign( m_map , rhs.m_map , rhs.m_track );
    }

  template< class RT , class ... RP >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RP...> & rhs )
    {
      typedef typename View<RT,RP...>::traits  SrcTraits ;
      typedef Kokkos::Impl::ViewMapping< traits , SrcTraits , void >  Mapping ;
      static_assert( Mapping::is_assignable , "Incompatible View copy assignment" );
      Mapping::assign( m_map , rhs.m_map , rhs.m_track );
      m_track.assign( rhs.m_track , traits::is_managed );
      return *this ;
    }

  //----------------------------------------
  // Compatible subview constructor
  // may assign unmanaged from managed.

  template< class RT , class ... RP , class Arg0 , class ... Args >
  KOKKOS_INLINE_FUNCTION
  View( const View< RT , RP... > & src_view
      , const Arg0 & arg0 , Args ... args )
    : m_track( src_view.m_track , traits::is_managed )
    , m_map()
    {
      typedef View< RT , RP... > SrcType ;

      typedef Kokkos::Impl::ViewMapping
        < void /* deduce destination view type from source view traits */
        , typename SrcType::traits
        , Arg0 , Args... > Mapping ;

      typedef typename Mapping::type DstType ;

      static_assert( Kokkos::Impl::ViewMapping< traits , typename DstType::traits , void >::is_assignable
        , "Subview construction requires compatible view and subview arguments" );

      Mapping::assign( m_map, src_view.m_map, arg0 , args... );
    }

  //----------------------------------------
  // Allocation tracking properties

  KOKKOS_INLINE_FUNCTION
  int use_count() const
    { return m_track.use_count(); }

  inline
  const std::string label() const
    { return m_track.template get_label< typename traits::memory_space >(); }

  //----------------------------------------
  // Allocation according to allocation properties and array layout

  template< class ... P >
  explicit inline
  View( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< ! Impl::ViewCtorProp< P... >::has_pointer
                               , typename traits::array_layout
                               >::type const & arg_layout
      )
    : m_track()
    , m_map()
    {
      // Append layout and spaces if not input
      typedef Impl::ViewCtorProp< P ... > alloc_prop_input ;

      // use 'std::integral_constant<unsigned,I>' for non-types
      // to avoid duplicate class error.
      typedef Impl::ViewCtorProp
        < P ...
        , typename std::conditional
            < alloc_prop_input::has_label
            , std::integral_constant<unsigned,0>
            , typename std::string
            >::type
        , typename std::conditional
            < alloc_prop_input::has_memory_space
            , std::integral_constant<unsigned,1>
            , typename traits::device_type::memory_space
            >::type
        , typename std::conditional
            < alloc_prop_input::has_execution_space
            , std::integral_constant<unsigned,2>
            , typename traits::device_type::execution_space
            >::type
        > alloc_prop ;

      static_assert( traits::is_managed
                   , "View allocation constructor requires managed memory" );

      if ( alloc_prop::initialize &&
           ! alloc_prop::execution_space::is_initialized() ) {
        // If initializing view data then
        // the execution space must be initialized.
        Kokkos::Impl::throw_runtime_exception("Constructing View and initializing data with uninitialized execution space");
      }

      // Copy the input allocation properties with possibly defaulted properties
      alloc_prop prop( arg_prop );

//------------------------------------------------------------
#if defined( KOKKOS_ENABLE_CUDA )
      // If allocating in CudaUVMSpace must fence before and after
      // the allocation to protect against possible concurrent access
      // on the CPU and the GPU.
      // Fence using the trait's executon space (which will be Kokkos::Cuda)
      // to avoid incomplete type errors from usng Kokkos::Cuda directly.
      if ( std::is_same< Kokkos::CudaUVMSpace , typename traits::device_type::memory_space >::value ) {
        traits::device_type::memory_space::execution_space::fence();
      }
#endif
//------------------------------------------------------------

      Kokkos::Impl::SharedAllocationRecord<> *
        record = m_map.allocate_shared( prop , arg_layout );

//------------------------------------------------------------
#if defined( KOKKOS_ENABLE_CUDA )
      if ( std::is_same< Kokkos::CudaUVMSpace , typename traits::device_type::memory_space >::value ) {
        traits::device_type::memory_space::execution_space::fence();
      }
#endif
//------------------------------------------------------------

      // Setup and initialization complete, start tracking
      m_track.assign_allocated_record_to_uninitialized( record );
    }

  KOKKOS_INLINE_FUNCTION
  void assign_data( pointer_type arg_data )
    {
      m_track.clear();
      m_map.assign_data( arg_data );
    }

  // Wrap memory according to properties and array layout
  template< class ... P >
  explicit KOKKOS_INLINE_FUNCTION
  View( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< Impl::ViewCtorProp< P... >::has_pointer
                               , typename traits::array_layout
                               >::type const & arg_layout
      )
    : m_track() // No memory tracking
    , m_map( arg_prop , arg_layout )
    {
      static_assert(
        std::is_same< pointer_type
                    , typename Impl::ViewCtorProp< P... >::pointer_type
                    >::value ,
        "Constructing View to wrap user memory must supply matching pointer type" );
    }

  // Simple dimension-only layout
  template< class ... P >
  explicit inline
  View( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< ! Impl::ViewCtorProp< P... >::has_pointer
                               , size_t
                               >::type const arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : View( arg_prop
          , typename traits::array_layout
              ( arg_N0 , arg_N1 , arg_N2 , arg_N3
              , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  template< class ... P >
  explicit KOKKOS_INLINE_FUNCTION
  View( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< Impl::ViewCtorProp< P... >::has_pointer
                               , size_t
                               >::type const arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : View( arg_prop
          , typename traits::array_layout
              ( arg_N0 , arg_N1 , arg_N2 , arg_N3
              , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  // Allocate with label and layout
  template< typename Label >
  explicit inline
  View( const Label & arg_label
      , typename std::enable_if<
          Kokkos::Impl::is_view_label<Label>::value ,
          typename traits::array_layout >::type const & arg_layout
      )
    : View( Impl::ViewCtorProp< std::string >( arg_label ) , arg_layout )
    {}

  // Allocate label and layout, must disambiguate from subview constructor.
  template< typename Label >
  explicit inline
  View( const Label & arg_label
      , typename std::enable_if<
          Kokkos::Impl::is_view_label<Label>::value ,
        const size_t >::type arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : View( Impl::ViewCtorProp< std::string >( arg_label )
          , typename traits::array_layout
              ( arg_N0 , arg_N1 , arg_N2 , arg_N3
              , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  // For backward compatibility
  explicit inline
  View( const ViewAllocateWithoutInitializing & arg_prop
      , const typename traits::array_layout & arg_layout
      )
    : View( Impl::ViewCtorProp< std::string , Kokkos::Impl::WithoutInitializing_t >( arg_prop.label , Kokkos::WithoutInitializing )
          , arg_layout
          )
    {}

  explicit inline
  View( const ViewAllocateWithoutInitializing & arg_prop
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : View( Impl::ViewCtorProp< std::string , Kokkos::Impl::WithoutInitializing_t >( arg_prop.label , Kokkos::WithoutInitializing )
          , typename traits::array_layout
              ( arg_N0 , arg_N1 , arg_N2 , arg_N3
              , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  //----------------------------------------
  // Memory span required to wrap these dimensions.
  static constexpr size_t required_allocation_size(
                                       const size_t arg_N0 = 0
                                     , const size_t arg_N1 = 0
                                     , const size_t arg_N2 = 0
                                     , const size_t arg_N3 = 0
                                     , const size_t arg_N4 = 0
                                     , const size_t arg_N5 = 0
                                     , const size_t arg_N6 = 0
                                     , const size_t arg_N7 = 0
                                     )
    {
      return map_type::memory_span(
        typename traits::array_layout
          ( arg_N0 , arg_N1 , arg_N2 , arg_N3
          , arg_N4 , arg_N5 , arg_N6 , arg_N7 ) );
    }

  explicit KOKKOS_INLINE_FUNCTION
  View( pointer_type arg_ptr
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : View( Impl::ViewCtorProp<pointer_type>(arg_ptr)
          , typename traits::array_layout
             ( arg_N0 , arg_N1 , arg_N2 , arg_N3
             , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  explicit KOKKOS_INLINE_FUNCTION
  View( pointer_type arg_ptr
      , const typename traits::array_layout & arg_layout
      )
    : View( Impl::ViewCtorProp<pointer_type>(arg_ptr) , arg_layout )
    {}

  //----------------------------------------
  // Shared scratch memory constructor

  static inline
  size_t shmem_size( const size_t arg_N0 = ~size_t(0) ,
                     const size_t arg_N1 = ~size_t(0) ,
                     const size_t arg_N2 = ~size_t(0) ,
                     const size_t arg_N3 = ~size_t(0) ,
                     const size_t arg_N4 = ~size_t(0) ,
                     const size_t arg_N5 = ~size_t(0) ,
                     const size_t arg_N6 = ~size_t(0) ,
                     const size_t arg_N7 = ~size_t(0) )
  {
    const size_t num_passed_args =
      ( arg_N0 != ~size_t(0) ) + ( arg_N1 != ~size_t(0) ) + ( arg_N2 != ~size_t(0) ) +
      ( arg_N3 != ~size_t(0) ) + ( arg_N4 != ~size_t(0) ) + ( arg_N5 != ~size_t(0) ) +
      ( arg_N6 != ~size_t(0) ) + ( arg_N7 != ~size_t(0) );

    if ( std::is_same<typename traits::specialize,void>::value && num_passed_args != traits::rank_dynamic ) {
      Kokkos::abort( "Kokkos::View::shmem_size() rank_dynamic != number of arguments.\n" );
    }

    return map_type::memory_span(
           typename traits::array_layout
            ( arg_N0 , arg_N1 , arg_N2 , arg_N3
            , arg_N4 , arg_N5 , arg_N6 , arg_N7 ) );
  }

  explicit KOKKOS_INLINE_FUNCTION
  View( const typename traits::execution_space::scratch_memory_space & arg_space
      , const typename traits::array_layout & arg_layout )
    : View( Impl::ViewCtorProp<pointer_type>(
              reinterpret_cast<pointer_type>(
                arg_space.get_shmem( map_type::memory_span( arg_layout ) ) ) )
         , arg_layout )
    {}

  explicit KOKKOS_INLINE_FUNCTION
  View( const typename traits::execution_space::scratch_memory_space & arg_space
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0 )
    : View( Impl::ViewCtorProp<pointer_type>(
              reinterpret_cast<pointer_type>(
                arg_space.get_shmem(
                  map_type::memory_span(
                    typename traits::array_layout
                     ( arg_N0 , arg_N1 , arg_N2 , arg_N3
                     , arg_N4 , arg_N5 , arg_N6 , arg_N7 ) ) ) ) )
          , typename traits::array_layout
             ( arg_N0 , arg_N1 , arg_N2 , arg_N3
             , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
       )
    {}
};


 /** \brief Temporary free function rank()
  *         until rank() is implemented
  *         in the View
  */
  template < typename D , class ... P >
  KOKKOS_INLINE_FUNCTION
  constexpr unsigned rank( const View<D , P...> & V ) { return V.Rank; } //Temporary until added to view

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class V , class ... Args >
using Subview =
  typename Kokkos::Impl::ViewMapping
    < void /* deduce subview type from source view traits */
    , typename V::traits
    , Args ...
    >::type ;

template< class D, class ... P , class ... Args >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::ViewMapping
  < void /* deduce subview type from source view traits */
  , ViewTraits< D , P... >
  , Args ...
  >::type
subview( const View< D, P... > & src , Args ... args )
{
  static_assert( View< D , P... >::Rank == sizeof...(Args) ,
    "subview requires one argument for each source View rank" );

  return typename
    Kokkos::Impl::ViewMapping
      < void /* deduce subview type from source view traits */
      , ViewTraits< D , P ... >
      , Args ... >::type( src , args ... );
}

template< class MemoryTraits , class D, class ... P , class ... Args >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::ViewMapping
  < void /* deduce subview type from source view traits */
  , ViewTraits< D , P... >
  , Args ...
  >::template apply< MemoryTraits >::type
subview( const View< D, P... > & src , Args ... args )
{
  static_assert( View< D , P... >::Rank == sizeof...(Args) ,
    "subview requires one argument for each source View rank" );

  return typename
    Kokkos::Impl::ViewMapping
      < void /* deduce subview type from source view traits */
      , ViewTraits< D , P ... >
      , Args ... >
      ::template apply< MemoryTraits >
      ::type( src , args ... );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class LT , class ... LP , class RT , class ... RP >
KOKKOS_INLINE_FUNCTION
bool operator == ( const View<LT,LP...> & lhs ,
                   const View<RT,RP...> & rhs )
{
  // Same data, layout, dimensions
  typedef ViewTraits<LT,LP...>  lhs_traits ;
  typedef ViewTraits<RT,RP...>  rhs_traits ;

  return
    std::is_same< typename lhs_traits::const_value_type ,
                  typename rhs_traits::const_value_type >::value &&
    std::is_same< typename lhs_traits::array_layout ,
                  typename rhs_traits::array_layout >::value &&
    std::is_same< typename lhs_traits::memory_space ,
                  typename rhs_traits::memory_space >::value &&
    unsigned(lhs_traits::rank) == unsigned(rhs_traits::rank) &&
    lhs.data()        == rhs.data() &&
    lhs.span()        == rhs.span() &&
    lhs.dimension_0() == rhs.dimension_0() &&
    lhs.dimension_1() == rhs.dimension_1() &&
    lhs.dimension_2() == rhs.dimension_2() &&
    lhs.dimension_3() == rhs.dimension_3() &&
    lhs.dimension_4() == rhs.dimension_4() &&
    lhs.dimension_5() == rhs.dimension_5() &&
    lhs.dimension_6() == rhs.dimension_6() &&
    lhs.dimension_7() == rhs.dimension_7();
}

template< class LT , class ... LP , class RT , class ... RP >
KOKKOS_INLINE_FUNCTION
bool operator != ( const View<LT,LP...> & lhs ,
                   const View<RT,RP...> & rhs )
{
  return ! ( operator==(lhs,rhs) );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

inline
void shared_allocation_tracking_disable()
{ Kokkos::Impl::SharedAllocationRecord<void,void>::tracking_disable(); }

inline
void shared_allocation_tracking_enable()
{ Kokkos::Impl::SharedAllocationRecord<void,void>::tracking_enable(); }

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class OutputView , typename Enable = void >
struct ViewFill {

  typedef typename OutputView::const_value_type  const_value_type ;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    const size_t n1 = output.dimension_1();
    const size_t n2 = output.dimension_2();
    const size_t n3 = output.dimension_3();
    const size_t n4 = output.dimension_4();
    const size_t n5 = output.dimension_5();
    const size_t n6 = output.dimension_6();
    const size_t n7 = output.dimension_7();

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_t i7 = 0 ; i7 < n7 ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      typedef typename OutputView::execution_space  execution_space ;
      typedef Kokkos::RangePolicy< execution_space > Policy ;

      const Kokkos::Impl::ParallelFor< ViewFill , Policy > closure( *this , Policy( 0 , output.dimension_0() ) );

      closure.execute();

      execution_space::fence();
    }
};

template< class OutputView >
struct ViewFill< OutputView , typename std::enable_if< OutputView::Rank == 0 >::type > {
  ViewFill( const OutputView & dst , const typename OutputView::const_value_type & src )
    {
      Kokkos::Impl::DeepCopy< typename OutputView::memory_space , Kokkos::HostSpace >
        ( dst.data() , & src , sizeof(typename OutputView::const_value_type) );
    }
};

template< class OutputView , class InputView , class ExecSpace = typename OutputView::execution_space >
struct ViewRemap {

  const OutputView output ;
  const InputView  input ;
  const size_t n0 ;
  const size_t n1 ;
  const size_t n2 ;
  const size_t n3 ;
  const size_t n4 ;
  const size_t n5 ;
  const size_t n6 ;
  const size_t n7 ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    , n0( std::min( (size_t)arg_out.dimension_0() , (size_t)arg_in.dimension_0() ) )
    , n1( std::min( (size_t)arg_out.dimension_1() , (size_t)arg_in.dimension_1() ) )
    , n2( std::min( (size_t)arg_out.dimension_2() , (size_t)arg_in.dimension_2() ) )
    , n3( std::min( (size_t)arg_out.dimension_3() , (size_t)arg_in.dimension_3() ) )
    , n4( std::min( (size_t)arg_out.dimension_4() , (size_t)arg_in.dimension_4() ) )
    , n5( std::min( (size_t)arg_out.dimension_5() , (size_t)arg_in.dimension_5() ) )
    , n6( std::min( (size_t)arg_out.dimension_6() , (size_t)arg_in.dimension_6() ) )
    , n7( std::min( (size_t)arg_out.dimension_7() , (size_t)arg_in.dimension_7() ) )
    {
      typedef Kokkos::RangePolicy< ExecSpace > Policy ;
      const Kokkos::Impl::ParallelFor< ViewRemap , Policy > closure( *this , Policy( 0 , n0 ) );
      closure.execute();
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_t i7 = 0 ; i7 < n7 ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Deep copy a value from Host memory into a view.  */
template< class DT , class ... DP >
inline
void deep_copy
  ( const View<DT,DP...> & dst
  , typename ViewTraits<DT,DP...>::const_value_type & value
  , typename std::enable_if<
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value
    >::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::non_const_value_type ,
                  typename ViewTraits<DT,DP...>::value_type >::value
    , "deep_copy requires non-const type" );

  Kokkos::Impl::ViewFill< View<DT,DP...> >( dst , value );
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template< class ST , class ... SP >
inline
void deep_copy
  ( typename ViewTraits<ST,SP...>::non_const_value_type & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value
    >::type * = 0 )
{
  static_assert( ViewTraits<ST,SP...>::rank == 0
               , "ERROR: Non-rank-zero view in deep_copy( value , View )" );

  typedef ViewTraits<ST,SP...>               src_traits ;
  typedef typename src_traits::memory_space  src_memory_space ;
  Kokkos::Impl::DeepCopy< HostSpace , src_memory_space >( & dst , src.data() , sizeof(ST) );
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of compatible type, and rank zero.  */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy
  ( const View<DT,DP...> & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<(
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value &&
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value &&
    ( unsigned(ViewTraits<DT,DP...>::rank) == unsigned(0) &&
      unsigned(ViewTraits<ST,SP...>::rank) == unsigned(0) )
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<ST,SP...>::non_const_value_type >::value
    , "deep_copy requires matching non-const destination type" );

  typedef View<DT,DP...>  dst_type ;
  typedef View<ST,SP...>  src_type ;

  typedef typename dst_type::value_type    value_type ;
  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( dst.data() != src.data() ) {
    Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , sizeof(value_type) );
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the default specialization, compatible type,
 *          same non-zero rank, same contiguous layout.
 */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy
  ( const View<DT,DP...> & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<(
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value &&
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value &&
    ( unsigned(ViewTraits<DT,DP...>::rank) != 0 ||
      unsigned(ViewTraits<ST,SP...>::rank) != 0 )
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "deep_copy requires non-const destination type" );

  static_assert(
    ( unsigned(ViewTraits<DT,DP...>::rank) ==
      unsigned(ViewTraits<ST,SP...>::rank) )
    , "deep_copy requires Views of equal rank" );

  typedef View<DT,DP...>  dst_type ;
  typedef View<ST,SP...>  src_type ;

  typedef typename dst_type::execution_space  dst_execution_space ;
  typedef typename src_type::execution_space  src_execution_space ;
  typedef typename dst_type::memory_space     dst_memory_space ;
  typedef typename src_type::memory_space     src_memory_space ;

  enum { DstExecCanAccessSrc =
   Kokkos::Impl::SpaceAccessibility< dst_execution_space , src_memory_space >::accessible };

  enum { SrcExecCanAccessDst =
   Kokkos::Impl::SpaceAccessibility< src_execution_space , dst_memory_space >::accessible };


  if ( (void *) dst.data() != (void*) src.data() ) {

#if defined(KOKKOS_ENABLE_PROFILING)
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();
      Kokkos::Profiling::beginDeepCopy(
          Kokkos::Profiling::SpaceHandle(dst_memory_space::name()),
          dst.label(),
          dst.data(),
          Kokkos::Profiling::SpaceHandle(src_memory_space::name()),
          src.label(),
          src.data(),
          nbytes);
    }
#endif

    // Concern: If overlapping views then a parallel copy will be erroneous.
    // ...

    // If same type, equal layout, equal dimensions, equal span, and contiguous memory then can byte-wise copy

    if ( std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                       typename ViewTraits<ST,SP...>::non_const_value_type >::value &&
         (
           ( std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                           typename ViewTraits<ST,SP...>::array_layout >::value
             &&
             ( std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                             typename Kokkos::LayoutLeft>::value
             ||
               std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                             typename Kokkos::LayoutRight>::value
             )
           )
           ||
           ( ViewTraits<DT,DP...>::rank == 1 &&
             ViewTraits<ST,SP...>::rank == 1 )
         ) &&
         dst.span_is_contiguous() &&
         src.span_is_contiguous() &&
         dst.span() == src.span() &&
         dst.dimension_0() == src.dimension_0() &&
         dst.dimension_1() == src.dimension_1() &&
         dst.dimension_2() == src.dimension_2() &&
         dst.dimension_3() == src.dimension_3() &&
         dst.dimension_4() == src.dimension_4() &&
         dst.dimension_5() == src.dimension_5() &&
         dst.dimension_6() == src.dimension_6() &&
         dst.dimension_7() == src.dimension_7() ) {

      const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();

      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
    }
    else if ( std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                            typename ViewTraits<ST,SP...>::non_const_value_type >::value &&
         (
           ( std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                           typename ViewTraits<ST,SP...>::array_layout >::value
             &&
             std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                          typename Kokkos::LayoutStride>::value
           )
           ||
           ( ViewTraits<DT,DP...>::rank == 1 &&
             ViewTraits<ST,SP...>::rank == 1 )
         ) &&
         dst.span_is_contiguous() &&
         src.span_is_contiguous() &&
         dst.span() == src.span() &&
         dst.dimension_0() == src.dimension_0() &&
         dst.dimension_1() == src.dimension_1() &&
         dst.dimension_2() == src.dimension_2() &&
         dst.dimension_3() == src.dimension_3() &&
         dst.dimension_4() == src.dimension_4() &&
         dst.dimension_5() == src.dimension_5() &&
         dst.dimension_6() == src.dimension_6() &&
         dst.dimension_7() == src.dimension_7() &&
         dst.stride_0() == src.stride_0() &&
         dst.stride_1() == src.stride_1() &&
         dst.stride_2() == src.stride_2() &&
         dst.stride_3() == src.stride_3() &&
         dst.stride_4() == src.stride_4() &&
         dst.stride_5() == src.stride_5() &&
         dst.stride_6() == src.stride_6() &&
         dst.stride_7() == src.stride_7()
         ) {

      const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();

      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space >( dst.data() , src.data() , nbytes );
    }
    else if ( DstExecCanAccessSrc ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::Impl::ViewRemap< dst_type , src_type >( dst , src );
    }
    else if ( SrcExecCanAccessDst ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::Impl::ViewRemap< dst_type , src_type , src_execution_space >( dst , src );
    }
    else {
      Kokkos::Impl::throw_runtime_exception("deep_copy given views that would require a temporary allocation");
    }

#if defined(KOKKOS_ENABLE_PROFILING)
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
#endif

  } // ( (void *) dst.data() != (void*) src.data() )
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Deep copy a value from Host memory into a view.  */
template< class ExecSpace ,class DT , class ... DP >
inline
void deep_copy
  ( const ExecSpace &
  , const View<DT,DP...> & dst
  , typename ViewTraits<DT,DP...>::const_value_type & value
  , typename std::enable_if<
    Kokkos::Impl::is_execution_space< ExecSpace >::value &&
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value
    >::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::non_const_value_type ,
                  typename ViewTraits<DT,DP...>::value_type >::value
    , "deep_copy requires non-const type" );

  Kokkos::Impl::ViewFill< View<DT,DP...> >( dst , value );
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template< class ExecSpace , class ST , class ... SP >
inline
void deep_copy
  ( const ExecSpace & exec_space
  , typename ViewTraits<ST,SP...>::non_const_value_type & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<
    Kokkos::Impl::is_execution_space< ExecSpace >::value &&
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value
    >::type * = 0 )
{
  static_assert( ViewTraits<ST,SP...>::rank == 0
               , "ERROR: Non-rank-zero view in deep_copy( value , View )" );

  typedef ViewTraits<ST,SP...>               src_traits ;
  typedef typename src_traits::memory_space  src_memory_space ;
  Kokkos::Impl::DeepCopy< HostSpace , src_memory_space , ExecSpace >
    ( exec_space , & dst , src.data() , sizeof(ST) );
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of compatible type, and rank zero.  */
template< class ExecSpace , class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy
  ( const ExecSpace & exec_space
  , const View<DT,DP...> & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<(
    Kokkos::Impl::is_execution_space< ExecSpace >::value &&
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value &&
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value &&
    ( unsigned(ViewTraits<DT,DP...>::rank) == unsigned(0) &&
      unsigned(ViewTraits<ST,SP...>::rank) == unsigned(0) )
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<ST,SP...>::non_const_value_type >::value
    , "deep_copy requires matching non-const destination type" );

  typedef View<DT,DP...>  dst_type ;
  typedef View<ST,SP...>  src_type ;

  typedef typename dst_type::value_type    value_type ;
  typedef typename dst_type::memory_space  dst_memory_space ;
  typedef typename src_type::memory_space  src_memory_space ;

  if ( dst.data() != src.data() ) {
    Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space , ExecSpace >
      ( exec_space , dst.data() , src.data() , sizeof(value_type) );
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the default specialization, compatible type,
 *          same non-zero rank, same contiguous layout.
 */
template< class ExecSpace , class DT, class ... DP, class ST, class ... SP >
inline
void deep_copy
  ( const ExecSpace & exec_space
  , const View<DT,DP...> & dst
  , const View<ST,SP...> & src
  , typename std::enable_if<(
    Kokkos::Impl::is_execution_space< ExecSpace >::value &&
    std::is_same< typename ViewTraits<DT,DP...>::specialize , void >::value &&
    std::is_same< typename ViewTraits<ST,SP...>::specialize , void >::value &&
    ( unsigned(ViewTraits<DT,DP...>::rank) != 0 ||
      unsigned(ViewTraits<ST,SP...>::rank) != 0 )
  )>::type * = 0 )
{
  static_assert(
    std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                  typename ViewTraits<DT,DP...>::non_const_value_type >::value
    , "deep_copy requires non-const destination type" );

  static_assert(
    ( unsigned(ViewTraits<DT,DP...>::rank) ==
      unsigned(ViewTraits<ST,SP...>::rank) )
    , "deep_copy requires Views of equal rank" );

  typedef View<DT,DP...>  dst_type ;
  typedef View<ST,SP...>  src_type ;

  typedef typename dst_type::execution_space  dst_execution_space ;
  typedef typename src_type::execution_space  src_execution_space ;
  typedef typename dst_type::memory_space     dst_memory_space ;
  typedef typename src_type::memory_space     src_memory_space ;

  enum { DstExecCanAccessSrc =
   Kokkos::Impl::SpaceAccessibility< dst_execution_space , src_memory_space >::accessible };

  enum { SrcExecCanAccessDst =
   Kokkos::Impl::SpaceAccessibility< src_execution_space , dst_memory_space >::accessible };

  if ( (void *) dst.data() != (void*) src.data() ) {

    // Concern: If overlapping views then a parallel copy will be erroneous.
    // ...

    // If same type, equal layout, equal dimensions, equal span, and contiguous memory then can byte-wise copy

    if ( std::is_same< typename ViewTraits<DT,DP...>::value_type ,
                       typename ViewTraits<ST,SP...>::non_const_value_type >::value &&
         (
           std::is_same< typename ViewTraits<DT,DP...>::array_layout ,
                         typename ViewTraits<ST,SP...>::array_layout >::value
           ||
           ( ViewTraits<DT,DP...>::rank == 1 &&
             ViewTraits<ST,SP...>::rank == 1 )
         ) &&
         dst.span_is_contiguous() &&
         src.span_is_contiguous() &&
         dst.span() == src.span() &&
         dst.dimension_0() == src.dimension_0() &&
         dst.dimension_1() == src.dimension_1() &&
         dst.dimension_2() == src.dimension_2() &&
         dst.dimension_3() == src.dimension_3() &&
         dst.dimension_4() == src.dimension_4() &&
         dst.dimension_5() == src.dimension_5() &&
         dst.dimension_6() == src.dimension_6() &&
         dst.dimension_7() == src.dimension_7() ) {

      const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();

      Kokkos::Impl::DeepCopy< dst_memory_space , src_memory_space , ExecSpace >
        ( exec_space , dst.data() , src.data() , nbytes );
    }
    else if ( DstExecCanAccessSrc ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::Impl::ViewRemap< dst_type , src_type >( dst , src );
    }
    else if ( SrcExecCanAccessDst ) {
      // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
      Kokkos::Impl::ViewRemap< dst_type , src_type , src_execution_space >( dst , src );
    }
    else {
      Kokkos::Impl::throw_runtime_exception("deep_copy given views that would require a temporary allocation");
    }
  }
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Deduce Mirror Types
template<class Space, class T, class ... P>
struct MirrorViewType {
  // The incoming view_type
  typedef typename Kokkos::View<T,P...> src_view_type;
  // The memory space for the mirror view
  typedef typename Space::memory_space memory_space;
  // Check whether it is the same memory space
  enum { is_same_memspace = std::is_same<memory_space,typename src_view_type::memory_space>::value };
  // The array_layout
  typedef typename src_view_type::array_layout array_layout;
  // The data type (we probably want it non-const since otherwise we can't even deep_copy to it.
  typedef typename src_view_type::non_const_data_type data_type;
  // The destination view type if it is not the same memory space
  typedef Kokkos::View<data_type,array_layout,Space> dest_view_type;
  // If it is the same memory_space return the existsing view_type
  // This will also keep the unmanaged trait if necessary
  typedef typename std::conditional<is_same_memspace,src_view_type,dest_view_type>::type view_type;
};

template<class Space, class T, class ... P>
struct MirrorType {
  // The incoming view_type
  typedef typename Kokkos::View<T,P...> src_view_type;
  // The memory space for the mirror view
  typedef typename Space::memory_space memory_space;
  // Check whether it is the same memory space
  enum { is_same_memspace = std::is_same<memory_space,typename src_view_type::memory_space>::value };
  // The array_layout
  typedef typename src_view_type::array_layout array_layout;
  // The data type (we probably want it non-const since otherwise we can't even deep_copy to it.
  typedef typename src_view_type::non_const_data_type data_type;
  // The destination view type if it is not the same memory space
  typedef Kokkos::View<data_type,array_layout,Space> view_type;
};

}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror( const Kokkos::View<T,P...> & src
             , typename std::enable_if<
                 ! std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout
                               , Kokkos::LayoutStride >::value
               >::type * = 0
             )
{
  typedef View<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  return dst_type( std::string( src.label() ).append("_mirror")
                 , src.dimension_0()
                 , src.dimension_1()
                 , src.dimension_2()
                 , src.dimension_3()
                 , src.dimension_4()
                 , src.dimension_5()
                 , src.dimension_6()
                 , src.dimension_7() );
}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror( const Kokkos::View<T,P...> & src
             , typename std::enable_if<
                 std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout
                             , Kokkos::LayoutStride >::value
               >::type * = 0
             )
{
  typedef View<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  Kokkos::LayoutStride layout ;

  layout.dimension[0] = src.dimension_0();
  layout.dimension[1] = src.dimension_1();
  layout.dimension[2] = src.dimension_2();
  layout.dimension[3] = src.dimension_3();
  layout.dimension[4] = src.dimension_4();
  layout.dimension[5] = src.dimension_5();
  layout.dimension[6] = src.dimension_6();
  layout.dimension[7] = src.dimension_7();

  layout.stride[0] = src.stride_0();
  layout.stride[1] = src.stride_1();
  layout.stride[2] = src.stride_2();
  layout.stride[3] = src.stride_3();
  layout.stride[4] = src.stride_4();
  layout.stride[5] = src.stride_5();
  layout.stride[6] = src.stride_6();
  layout.stride[7] = src.stride_7();

  return dst_type( std::string( src.label() ).append("_mirror") , layout );
}


// Create a mirror in a new space (specialization for different space)
template<class Space, class T, class ... P>
typename Impl::MirrorType<Space,T,P ...>::view_type create_mirror(const Space& , const Kokkos::View<T,P...> & src) {
  return typename Impl::MirrorType<Space,T,P ...>::view_type(src.label(),src.layout());
}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror_view( const Kokkos::View<T,P...> & src
                  , typename std::enable_if<(
                      std::is_same< typename Kokkos::View<T,P...>::memory_space
                                  , typename Kokkos::View<T,P...>::HostMirror::memory_space
                                  >::value
                      &&
                      std::is_same< typename Kokkos::View<T,P...>::data_type
                                  , typename Kokkos::View<T,P...>::HostMirror::data_type
                                  >::value
                    )>::type * = 0
                  )
{
  return src ;
}

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror_view( const Kokkos::View<T,P...> & src
                  , typename std::enable_if< ! (
                      std::is_same< typename Kokkos::View<T,P...>::memory_space
                                  , typename Kokkos::View<T,P...>::HostMirror::memory_space
                                  >::value
                      &&
                      std::is_same< typename Kokkos::View<T,P...>::data_type
                                  , typename Kokkos::View<T,P...>::HostMirror::data_type
                                  >::value
                    )>::type * = 0
                  )
{
  return Kokkos::create_mirror( src );
}

// Create a mirror view in a new space (specialization for same space)
template<class Space, class T, class ... P>
typename Impl::MirrorViewType<Space,T,P ...>::view_type
create_mirror_view(const Space& , const Kokkos::View<T,P...> & src
  , typename std::enable_if<Impl::MirrorViewType<Space,T,P ...>::is_same_memspace>::type* = 0 ) {
  return src;
}

// Create a mirror view in a new space (specialization for different space)
template<class Space, class T, class ... P>
typename Impl::MirrorViewType<Space,T,P ...>::view_type
create_mirror_view(const Space& , const Kokkos::View<T,P...> & src
  , typename std::enable_if<!Impl::MirrorViewType<Space,T,P ...>::is_same_memspace>::type* = 0 ) {
  return typename Impl::MirrorViewType<Space,T,P ...>::view_type(src.label(),src.layout());
}

// Create a mirror view and deep_copy in a new space (specialization for same space)
template<class Space, class T, class ... P>
typename Impl::MirrorViewType<Space,T,P ...>::view_type
create_mirror_view_and_copy(const Space& , const Kokkos::View<T,P...> & src
  , std::string const& name = ""
  , typename std::enable_if<Impl::MirrorViewType<Space,T,P ...>::is_same_memspace>::type* = 0 ) {
  (void)name;
  return src;
}

// Create a mirror view and deep_copy in a new space (specialization for different space)
template<class Space, class T, class ... P>
typename Impl::MirrorViewType<Space,T,P ...>::view_type
create_mirror_view_and_copy(const Space& , const Kokkos::View<T,P...> & src
  , std::string const& name = ""
  , typename std::enable_if<!Impl::MirrorViewType<Space,T,P ...>::is_same_memspace>::type* = 0 ) {
  using Mirror = typename Impl::MirrorViewType<Space,T,P ...>::view_type;
  std::string label = name.empty() ? src.label() : name;
  auto mirror = Mirror(ViewAllocateWithoutInitializing(label), src.layout());
  deep_copy(mirror, src);
  return mirror;
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Resize a view with copying old data to new data at the corresponding indices. */
template< class T , class ... P >
inline
typename std::enable_if<
  std::is_same<typename Kokkos::View<T,P...>::array_layout,Kokkos::LayoutLeft>::value ||
  std::is_same<typename Kokkos::View<T,P...>::array_layout,Kokkos::LayoutRight>::value
>::type
resize( Kokkos::View<T,P...> & v ,
             const size_t n0 = 0 ,
             const size_t n1 = 0 ,
             const size_t n2 = 0 ,
             const size_t n3 = 0 ,
             const size_t n4 = 0 ,
             const size_t n5 = 0 ,
             const size_t n6 = 0 ,
             const size_t n7 = 0 )
{
  typedef Kokkos::View<T,P...>  view_type ;

  static_assert( Kokkos::ViewTraits<T,P...>::is_managed , "Can only resize managed views" );

  // Fix #904 by checking dimensions before actually resizing.
  //
  // Rank is known at compile time, so hopefully the compiler will
  // remove branches that are compile-time false.  The upcoming "if
  // constexpr" language feature would make this certain.
  if (view_type::Rank == 1 &&
      n0 == static_cast<size_t> (v.extent(0))) {
    return;
  }
  if (view_type::Rank == 2 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1))) {
    return;
  }
  if (view_type::Rank == 3 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2))) {
    return;
  }
  if (view_type::Rank == 4 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2)) &&
      n3 == static_cast<size_t> (v.extent(3))) {
    return;
  }
  if (view_type::Rank == 5 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2)) &&
      n3 == static_cast<size_t> (v.extent(3)) &&
      n4 == static_cast<size_t> (v.extent(4))) {
    return;
  }
  if (view_type::Rank == 6 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2)) &&
      n3 == static_cast<size_t> (v.extent(3)) &&
      n4 == static_cast<size_t> (v.extent(4)) &&
      n5 == static_cast<size_t> (v.extent(5))) {
    return;
  }
  if (view_type::Rank == 7 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2)) &&
      n3 == static_cast<size_t> (v.extent(3)) &&
      n4 == static_cast<size_t> (v.extent(4)) &&
      n5 == static_cast<size_t> (v.extent(5)) &&
      n6 == static_cast<size_t> (v.extent(6))) {
    return;
  }
  if (view_type::Rank == 8 &&
      n0 == static_cast<size_t> (v.extent(0)) &&
      n1 == static_cast<size_t> (v.extent(1)) &&
      n2 == static_cast<size_t> (v.extent(2)) &&
      n3 == static_cast<size_t> (v.extent(3)) &&
      n4 == static_cast<size_t> (v.extent(4)) &&
      n5 == static_cast<size_t> (v.extent(5)) &&
      n6 == static_cast<size_t> (v.extent(6)) &&
      n7 == static_cast<size_t> (v.extent(7))) {
    return;
  }
  // If Kokkos ever supports Views of rank > 8, the above code won't
  // be incorrect, because avoiding reallocation in resize() is just
  // an optimization.

  // TODO (mfh 27 Jun 2017) If the old View has enough space but just
  // different dimensions (e.g., if the product of the dimensions,
  // including extra space for alignment, will not change), then
  // consider just reusing storage.  For now, Kokkos always
  // reallocates if any of the dimensions change, even if the old View
  // has enough space.

  view_type v_resized( v.label(), n0, n1, n2, n3, n4, n5, n6, n7 );

  Kokkos::Impl::ViewRemap< view_type , view_type >( v_resized , v );

  v = v_resized ;
}

/** \brief  Resize a view with copying old data to new data at the corresponding indices. */
template< class T , class ... P >
inline
void resize(       Kokkos::View<T,P...> & v ,
    const typename Kokkos::View<T,P...>::array_layout & layout)
{
  typedef Kokkos::View<T,P...>  view_type ;

  static_assert( Kokkos::ViewTraits<T,P...>::is_managed , "Can only resize managed views" );

  view_type v_resized( v.label(), layout );

  Kokkos::Impl::ViewRemap< view_type , view_type >( v_resized , v );

  v = v_resized ;
}

/** \brief  Resize a view with discarding old data. */
template< class T , class ... P >
inline
typename std::enable_if<
  std::is_same<typename Kokkos::View<T,P...>::array_layout,Kokkos::LayoutLeft>::value ||
  std::is_same<typename Kokkos::View<T,P...>::array_layout,Kokkos::LayoutRight>::value
>::type
realloc( Kokkos::View<T,P...> & v ,
              const size_t n0 = 0 ,
              const size_t n1 = 0 ,
              const size_t n2 = 0 ,
              const size_t n3 = 0 ,
              const size_t n4 = 0 ,
              const size_t n5 = 0 ,
              const size_t n6 = 0 ,
              const size_t n7 = 0 )
{
  typedef Kokkos::View<T,P...>  view_type ;

  static_assert( Kokkos::ViewTraits<T,P...>::is_managed , "Can only realloc managed views" );

  const std::string label = v.label();

  v = view_type(); // Deallocate first, if the only view to allocation
  v = view_type( label, n0, n1, n2, n3, n4, n5, n6, n7 );
}

/** \brief  Resize a view with discarding old data. */
template< class T , class ... P >
inline
void realloc(      Kokkos::View<T,P...> & v ,
    const typename Kokkos::View<T,P...>::array_layout & layout)
{
  typedef Kokkos::View<T,P...>  view_type ;

  static_assert( Kokkos::ViewTraits<T,P...>::is_managed , "Can only realloc managed views" );

  const std::string label = v.label();

  v = view_type(); // Deallocate first, if the only view to allocation
  v = view_type( label, layout );
}
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos { namespace Impl {

template < class Specialize, typename A, typename B >
struct CommonViewValueType;

template < typename A, typename B >
struct CommonViewValueType< void, A, B >
{
  using value_type = typename std::common_type< A , B >::type;
};


template < class Specialize, class ValueType >
struct CommonViewAllocProp;

template < class ValueType >
struct CommonViewAllocProp< void, ValueType >
{
  using value_type = ValueType;
  using scalar_array_type = ValueType;

  template < class ... Views >
  KOKKOS_INLINE_FUNCTION
  CommonViewAllocProp( const Views & ... ) {}
};


template < class ... Views >
struct DeduceCommonViewAllocProp;

// Base case must provide types for:
// 1. specialize  2. value_type  3. is_view  4. prop_type
template < class FirstView >
struct DeduceCommonViewAllocProp< FirstView >
{
  using specialize = typename FirstView::traits::specialize;

  using value_type = typename FirstView::traits::value_type;

  enum : bool { is_view = is_view< FirstView >::value };

  using prop_type = CommonViewAllocProp< specialize, value_type >;
};


template < class FirstView, class ... NextViews >
struct DeduceCommonViewAllocProp< FirstView, NextViews... >
{
  using NextTraits = DeduceCommonViewAllocProp< NextViews... >;

  using first_specialize = typename FirstView::traits::specialize;
  using first_value_type = typename FirstView::traits::value_type;

  enum : bool { first_is_view = is_view< FirstView >::value };

  using next_specialize = typename NextTraits::specialize;
  using next_value_type = typename NextTraits::value_type;

  enum : bool { next_is_view = NextTraits::is_view };

  // common types

  // determine specialize type
  // if first and next specialize differ, but are not the same specialize, error out
  static_assert( !(!std::is_same< first_specialize, next_specialize >::value && !std::is_same< first_specialize, void>::value && !std::is_same< void, next_specialize >::value)  , "Kokkos DeduceCommonViewAllocProp ERROR: Only one non-void specialize trait allowed" );

  // otherwise choose non-void specialize if either/both are non-void
  using specialize = typename std::conditional< std::is_same< first_specialize, next_specialize >::value
                                              , first_specialize
                                              , typename std::conditional< ( std::is_same< first_specialize, void >::value
                                                                             && !std::is_same< next_specialize, void >::value)
                                                                           , next_specialize
                                                                           , first_specialize
                                                                         >::type
                                               >::type;

  using value_type = typename CommonViewValueType< specialize, first_value_type, next_value_type >::value_type;

  enum : bool { is_view = (first_is_view && next_is_view) };

  using prop_type = CommonViewAllocProp< specialize, value_type >;
};

} // end namespace Impl

template < class ... Views >
using DeducedCommonPropsType = typename Impl::DeduceCommonViewAllocProp<Views...>::prop_type ;

// User function
template < class ... Views >
KOKKOS_INLINE_FUNCTION
DeducedCommonPropsType<Views...> 
common_view_alloc_prop( Views const & ... views )
{
  return DeducedCommonPropsType<Views...>( views... );
}

} // namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// For backward compatibility:

namespace Kokkos {
namespace Experimental {

using Kokkos::ViewTraits ;
using Kokkos::View ;
using Kokkos::Subview ;
using Kokkos::is_view ;
using Kokkos::subview ;
using Kokkos::ALL ;
using Kokkos::WithoutInitializing ;
using Kokkos::AllowPadding ;
using Kokkos::view_alloc ;
using Kokkos::view_wrap ;
using Kokkos::deep_copy ;
using Kokkos::create_mirror ;
using Kokkos::create_mirror_view ;
using Kokkos::resize ;
using Kokkos::realloc ;

namespace Impl {

using Kokkos::Impl::ViewFill ;
using Kokkos::Impl::ViewRemap ;
using Kokkos::Impl::ViewCtorProp ;
using Kokkos::Impl::is_view_label ;
using Kokkos::Impl::WithoutInitializing_t ;
using Kokkos::Impl::AllowPadding_t ;
using Kokkos::Impl::SharedAllocationRecord ;
using Kokkos::Impl::SharedAllocationTracker ;
using Kokkos::Impl::ViewMapping ;
using Kokkos::Impl::ViewDataAnalysis ;


} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

namespace Kokkos {
namespace Impl {

using Kokkos::is_view ;

template< class SrcViewType
        , class Arg0Type
        , class Arg1Type
        , class Arg2Type
        , class Arg3Type
        , class Arg4Type
        , class Arg5Type
        , class Arg6Type
        , class Arg7Type
        >
struct ViewSubview /* { typedef ... type ; } */ ;

} /* namespace Impl */
} /* namespace Kokkos */

#include <impl/Kokkos_Atomic_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEW_HPP */

