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

/// \file Kokkos_DynRankView.hpp
/// \brief Declaration and definition of Kokkos::Experimental::DynRankView.
///
/// This header file declares and defines Kokkos::Experimental::DynRankView and its
/// related nonmember functions.
/*
 *   Changes from View
 *   1. The rank of the DynRankView is returned by the method rank()
 *   2. Max rank of a DynRankView is 7
 *   3. subview name is subdynrankview
 *   4. Every subdynrankview is returned with LayoutStride
 */

#ifndef KOKKOS_DYNRANKVIEW_HPP
#define KOKKOS_DYNRANKVIEW_HPP

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <type_traits>

namespace Kokkos {
namespace Experimental {

namespace Impl {

template <typename Specialize>
struct DynRankDimTraits {

  // Compute the rank of the view from the nonzero dimension arguments.
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const size_t N0
                           , const size_t N1
                           , const size_t N2
                           , const size_t N3
                           , const size_t N4
                           , const size_t N5
                           , const size_t N6
                           , const size_t N7 )
  {
    return
      (   (N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0 && N1 == 0 && N0 == 0) ? 0
      : ( (N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0 && N1 == 0) ? 1
      : ( (N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0) ? 2
      : ( (N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0) ? 3
      : ( (N6 == 0 && N5 == 0 && N4 == 0) ? 4
      : ( (N6 == 0 && N5 == 0) ? 5
      : ( (N6 == 0) ? 6
      : 7 ) ) ) ) ) ) );
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank( const Layout& layout )
  {
    return computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );
  }

  // Create the layout for the rank-7 view.
  // Non-strided Layout
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutRight>::value || std::is_same<Layout , Kokkos::LayoutLeft>::value) , Layout >::type createLayout( const Layout& layout )
  {
    return Layout( layout.dimension[0] != 0 ? layout.dimension[0] : 1
                 , layout.dimension[1] != 0 ? layout.dimension[1] : 1
                 , layout.dimension[2] != 0 ? layout.dimension[2] : 1
                 , layout.dimension[3] != 0 ? layout.dimension[3] : 1
                 , layout.dimension[4] != 0 ? layout.dimension[4] : 1
                 , layout.dimension[5] != 0 ? layout.dimension[5] : 1
                 , layout.dimension[6] != 0 ? layout.dimension[6] : 1
                 , layout.dimension[7] != 0 ? layout.dimension[7] : 1
                 );
  }

  // LayoutStride
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutStride>::value) , Layout>::type createLayout( const Layout& layout )
  {
    return Layout( layout.dimension[0] != 0 ? layout.dimension[0] : 1
                 , layout.stride[0] 
                 , layout.dimension[1] != 0 ? layout.dimension[1] : 1
                 , layout.stride[1] 
                 , layout.dimension[2] != 0 ? layout.dimension[2] : 1
                 , layout.stride[2] 
                 , layout.dimension[3] != 0 ? layout.dimension[3] : 1
                 , layout.stride[3] 
                 , layout.dimension[4] != 0 ? layout.dimension[4] : 1
                 , layout.stride[4] 
                 , layout.dimension[5] != 0 ? layout.dimension[5] : 1
                 , layout.stride[5] 
                 , layout.dimension[6] != 0 ? layout.dimension[6] : 1
                 , layout.stride[6] 
                 , layout.dimension[7] != 0 ? layout.dimension[7] : 1
                 , layout.stride[7] 
                 );
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  template <typename ViewType, typename ViewArg>
  static ViewType createView( const ViewArg& arg
                            , const size_t N0
                            , const size_t N1
                            , const size_t N2
                            , const size_t N3
                            , const size_t N4
                            , const size_t N5
                            , const size_t N6
                            , const size_t N7 )
  {
    return ViewType( arg
                   , N0 != 0 ? N0 : 1
                   , N1 != 0 ? N1 : 1
                   , N2 != 0 ? N2 : 1
                   , N3 != 0 ? N3 : 1
                   , N4 != 0 ? N4 : 1
                   , N5 != 0 ? N5 : 1
                   , N6 != 0 ? N6 : 1
                   , N7 != 0 ? N7 : 1 );
  }
};

} //end Impl

/* \class DynRankView
 * \brief Container that creates a Kokkos view with rank determined at runtime. 
 *   Essentially this is a rank 7 view that wraps the access operators
 *   to yield the functionality of a view 
 *
 *   Changes from View
 *   1. The rank of the DynRankView is returned by the method rank()
 *   2. Max rank of a DynRankView is 7
 *   3. subview name is subdynrankview
 *   4. Every subdynrankview is returned with LayoutStride
 *
 */

template< typename DataType , class ... Properties >
class DynRankView : private View< DataType*******, Properties... >
{
  static_assert( !std::is_array<DataType>::value && !std::is_pointer<DataType>::value , "Cannot template DynRankView with array or pointer datatype - must be pod" );

public: 
  using view_type = View< DataType******* , Properties...>;
  using reference_type = typename view_type::reference_type; 

private: 
  template < class , class ... > friend class DynRankView ;
  template< class , class ... > friend class Impl::ViewMapping ;
  unsigned m_rank;

public:
  KOKKOS_INLINE_FUNCTION
  view_type & DownCast() const { return static_cast< view_type & > (*this); }
  KOKKOS_INLINE_FUNCTION
  const view_type & ConstDownCast() const { return static_cast< const view_type & > (*this); }

  typedef ViewTraits< DataType , Properties ... > traits ;

  // Data type traits:
  typedef typename traits::data_type            data_type;
  typedef typename traits::const_data_type      const_data_type;
  typedef typename traits::non_const_data_type  non_const_data_type;

  // Compatible array of trivial type traits:
  typedef typename traits::scalar_array_type            scalar_array_type ;
  typedef typename traits::const_scalar_array_type      const_scalar_array_type ;
  typedef typename traits::non_const_scalar_array_type  non_const_scalar_array_type ;

  // Value type traits:
  typedef typename traits::value_type            value_type ;
  typedef typename traits::const_value_type      const_value_type ;
  typedef typename traits::non_const_value_type  non_const_value_type ;

  // Mapping traits:
  typedef typename traits::array_layout   array_layout ;
  typedef typename traits::specialize     specialize ;

  // Execution space, memory space, memory access traits, and host mirror space:
  typedef typename traits::execution_space    execution_space ;
  typedef typename traits::memory_space       memory_space ;
  typedef typename traits::device_type        device_type ;
  typedef typename traits::memory_traits      memory_traits ;
  typedef typename traits::host_mirror_space  host_mirror_space ;

  typedef typename traits::size_type size_type ;

  using view_type::is_hostspace ;
  using view_type::is_managed ;
  using view_type::is_random_access ;

  /** \brief  Compatible view of array of scalar types */
  typedef DynRankView< typename traits::scalar_array_type ,
                       typename traits::array_layout ,
                       typename traits::device_type ,
                       typename traits::memory_traits >
    array_type ;

  /** \brief  Compatible view of const data type */
  typedef DynRankView< typename traits::const_data_type ,
                       typename traits::array_layout ,
                       typename traits::device_type ,
                       typename traits::memory_traits >
    const_type ;

  /** \brief  Compatible view of non-const data type */
  typedef DynRankView< typename traits::non_const_data_type ,
                       typename traits::array_layout ,
                       typename traits::device_type ,
                       typename traits::memory_traits >
    non_const_type ;

  /** \brief  Compatible HostMirror view */
  typedef DynRankView< typename traits::non_const_data_type ,
                       typename traits::array_layout ,
                       typename traits::host_mirror_space >
    HostMirror ;

  //----------------------------------------
  // Domain rank and extents

  KOKKOS_INLINE_FUNCTION
  DynRankView() : view_type() , m_rank(0) {}

  KOKKOS_INLINE_FUNCTION
  constexpr unsigned rank() const { return m_rank; }

  using view_type::extent; 
  using view_type::extent_int; 
  using view_type::layout;
  using view_type::dimension;
  using view_type::size;
  using view_type::stride;

  using pointer_type = typename view_type::pointer_type;
  using view_type::reference_type_is_lvalue_reference;
  using view_type::span;
  using view_type::capacity;
  using view_type::span_is_contiguous;
  using view_type::data;
  using view_type::implementation_map;

  using view_type::is_contiguous;
  using view_type::ptr_on_device;

  //Deprecated, remove soon (add for test)
  using view_type::dimension_0;
  using view_type::dimension_1;
  using view_type::dimension_2;
  using view_type::dimension_3;
  using view_type::dimension_4;
  using view_type::dimension_5;
  using view_type::dimension_6;
  using view_type::dimension_7;
  using view_type::stride_0;
  using view_type::stride_1;
  using view_type::stride_2;
  using view_type::stride_3;
  using view_type::stride_4;
  using view_type::stride_5;
  using view_type::stride_6;
  using view_type::stride_7;

  //operators ()
  // Rank 0
  KOKKOS_INLINE_FUNCTION
  reference_type operator()() const
    { return view_type::operator()(0,0,0,0,0,0,0); }
  
  // Rank 1
  // This assumes a contiguous underlying memory (i.e. no padding, no striding...)
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if< std::is_same<value_type, scalar_array_type>::value && std::is_integral<iType>::value, reference_type>::type
  operator[](const iType & i0) const
    {
      return data()[i0];
    }

  // This assumes a contiguous underlying memory (i.e. no padding, no striding...
  // AND a Trilinos/Sacado scalar type )
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if< !std::is_same<value_type, scalar_array_type>::value && std::is_integral<iType>::value, reference_type>::type
  operator[](const iType & i0) const
    {
      auto map = implementation_map();

      const size_t dim_scalar = map.dimension_scalar();
      const size_t bytes = this->span() / dim_scalar;

      typedef Kokkos::View<DataType*, array_layout, device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged | memory_traits::RandomAccess | memory_traits::Atomic> > tmp_view_type;
      tmp_view_type rankone_view(this->data(), bytes, dim_scalar);
      return rankone_view(i0);
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType & i0 ) const 
    { return view_type::operator()(i0,0,0,0,0,0,0); }

  // Rank 2
  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 ) const 
    { return view_type::operator()(i0,i1,0,0,0,0,0); }

  // Rank 3
  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const 
    { return view_type::operator()(i0,i1,i2,0,0,0,0); }

  // Rank 4
  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const 
    { return view_type::operator()(i0,i1,i2,i3,0,0,0); }

  // Rank 5
  template< typename iType0 , typename iType1 , typename iType2 , typename iType3, typename iType4 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 ) const 
    { return view_type::operator()(i0,i1,i2,i3,i4,0,0); }

  // Rank 6
  template< typename iType0 , typename iType1 , typename iType2 , typename iType3, typename iType4 , typename iType5 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 , const iType5 & i5 ) const 
    { return view_type::operator()(i0,i1,i2,i3,i4,i5,0); }

  // Rank 7
  template< typename iType0 , typename iType1 , typename iType2 , typename iType3, typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()(const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 , const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const 
    { return view_type::operator()(i0,i1,i2,i3,i4,i5,i6); }

  //----------------------------------------
  // Standard constructor, destructor, and assignment operators... 

  KOKKOS_INLINE_FUNCTION
  ~DynRankView() {}

  KOKKOS_INLINE_FUNCTION
  DynRankView( const DynRankView & ) = default ;

  KOKKOS_INLINE_FUNCTION
  DynRankView( DynRankView && ) = default ;

  KOKKOS_INLINE_FUNCTION
  DynRankView & operator = ( const DynRankView & ) = default ;

  KOKKOS_INLINE_FUNCTION
  DynRankView & operator = ( DynRankView && ) = default ;

  //----------------------------------------
  // Compatible view copy constructor and assignment
  // may assign unmanaged from managed.

  template< class RT , class ... RP >
  KOKKOS_INLINE_FUNCTION
  DynRankView( const DynRankView<RT,RP...> & rhs )
    : view_type( rhs.ConstDownCast() )  
    , m_rank(rhs.m_rank)
    {}

  template< class RT , class ... RP >
  KOKKOS_INLINE_FUNCTION
  DynRankView & operator = (const DynRankView<RT,RP...> & rhs )
  {
    view_type::operator = ( rhs.ConstDownCast() );
    m_rank = rhs.rank();
    return *this;
  }

  //----------------------------------------
  // Allocation tracking properties

  using view_type::use_count;
  using view_type::label;

  //----------------------------------------
  // Allocation according to allocation properties and array layout

  template< class ... P >
  explicit inline
  DynRankView( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< ! Impl::ViewCtorProp< P... >::has_pointer
                               , typename traits::array_layout
                               >::type const & arg_layout
      )
      : view_type( arg_prop
                 , Impl::DynRankDimTraits<typename traits::specialize>::createLayout(arg_layout) )
      , m_rank( Impl::DynRankDimTraits<typename traits::specialize>::computeRank(arg_layout) )
    {}

//Wrappers
  template< class ... P >
  explicit KOKKOS_INLINE_FUNCTION
  DynRankView( const Impl::ViewCtorProp< P ... > & arg_prop
      , typename std::enable_if< Impl::ViewCtorProp< P... >::has_pointer
                               , typename traits::array_layout
                               >::type const & arg_layout
      )
      : view_type( arg_prop
                 , Impl::DynRankDimTraits<typename traits::specialize>::createLayout(arg_layout) )
      , m_rank( Impl::DynRankDimTraits<typename traits::specialize>::computeRank(arg_layout) )
    {}

  //----------------------------------------
  //Constructor(s)

  // Simple dimension-only layout
  template< class ... P >
  explicit inline
  DynRankView( const Impl::ViewCtorProp< P ... > & arg_prop
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
    : DynRankView( arg_prop
    , typename traits::array_layout
          ( arg_N0 , arg_N1 , arg_N2 , arg_N3 , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
      )
    {}

  template< class ... P >
  explicit KOKKOS_INLINE_FUNCTION
  DynRankView( const Impl::ViewCtorProp< P ... > & arg_prop
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
    : DynRankView( arg_prop
    , typename traits::array_layout
          ( arg_N0 , arg_N1 , arg_N2 , arg_N3 , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
      )
    {}

  // Allocate with label and layout
  template< typename Label >
  explicit inline
  DynRankView( const Label & arg_label
      , typename std::enable_if<
          Kokkos::Experimental::Impl::is_view_label<Label>::value ,
          typename traits::array_layout >::type const & arg_layout
      )
    : DynRankView( Impl::ViewCtorProp< std::string >( arg_label ) , arg_layout )
    {}

  // Allocate label and layout
  template< typename Label >
  explicit inline
  DynRankView( const Label & arg_label
      , typename std::enable_if<
          Kokkos::Experimental::Impl::is_view_label<Label>::value ,
        const size_t >::type arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : DynRankView( Impl::ViewCtorProp< std::string >( arg_label )
    , typename traits::array_layout
          ( arg_N0 , arg_N1 , arg_N2 , arg_N3 , arg_N4 , arg_N5 , arg_N6 , arg_N7 )
          )
    {}

  // For backward compatibility
/*
  explicit inline
  DynRankView( const ViewAllocateWithoutInitializing & arg_prop
      , const typename traits::array_layout & arg_layout
      )
    : view_type( Impl::ViewCtorProp< std::string , Kokkos::Experimental::Impl::WithoutInitializing_t >( arg_prop.label , Kokkos::Experimental::WithoutInitializing )
          , arg_layout
          )
    //, m_rank(arg_N0 == 0 ? 0 : ( arg_N1 == 0 ? 1 : ( arg_N2 == 0 ? 2 : ( arg_N3 == 0 ? 3 : ( arg_N4 == 0 ? 4 : ( arg_N5 == 0 ? 5 : ( arg_N6 == 0 ? 6 : ( arg_N7 == 0 ? 7 : 8 ) ) ) ) ) ) ) ) //how to extract rank?
    {}
*/

  explicit inline
  DynRankView( const ViewAllocateWithoutInitializing & arg_prop
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : DynRankView(Impl::ViewCtorProp< std::string , Kokkos::Experimental::Impl::WithoutInitializing_t >( arg_prop.label , Kokkos::Experimental::WithoutInitializing ), arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7 ) 
    {}

  using view_type::memory_span;

  explicit KOKKOS_INLINE_FUNCTION
  DynRankView( pointer_type arg_ptr
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0
      )
    : DynRankView( Impl::ViewCtorProp<pointer_type>(arg_ptr) , arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7 )
    {}

  explicit KOKKOS_INLINE_FUNCTION
  DynRankView( pointer_type arg_ptr
      , typename traits::array_layout & arg_layout
      )
    : DynRankView( Impl::ViewCtorProp<pointer_type>(arg_ptr) , arg_layout )
    {}


  //----------------------------------------
  // Shared scratch memory constructor

  using view_type::shmem_size; 

  explicit KOKKOS_INLINE_FUNCTION
  DynRankView( const typename traits::execution_space::scratch_memory_space & arg_space
      , const size_t arg_N0 = 0
      , const size_t arg_N1 = 0
      , const size_t arg_N2 = 0
      , const size_t arg_N3 = 0
      , const size_t arg_N4 = 0
      , const size_t arg_N5 = 0
      , const size_t arg_N6 = 0
      , const size_t arg_N7 = 0 )
    : view_type( Impl::DynRankDimTraits<typename traits::specialize>::template createView<view_type>( arg_space
                                                                                                    , arg_N0
                                                                                                    , arg_N1
                                                                                                    , arg_N2
                                                                                                    , arg_N3
                                                                                                    , arg_N4
                                                                                                    , arg_N5
                                                                                                    , arg_N6
                                                                                                    , arg_N7 ) )
    , m_rank( Impl::DynRankDimTraits<typename traits::specialize>::computeRank( arg_N0
                                                                              , arg_N1
                                                                              , arg_N2
                                                                              , arg_N3
                                                                              , arg_N4
                                                                              , arg_N5
                                                                              , arg_N6
                                                                              , arg_N7 ) )
    {}

};

//----------------------------------------------------------------------------
// Subview mapping.
// Deduce destination view type from source view traits and subview arguments

namespace Impl {

struct DynRankSubviewTag {};

template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize , void >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      ) 
    ), DynRankSubviewTag >::type
  , SrcTraits
  , Args ... >
{
private:

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Args...>::value)
    , R2 = bool(is_integral_extent<2,Args...>::value)
    , R3 = bool(is_integral_extent<3,Args...>::value)
    , R4 = bool(is_integral_extent<4,Args...>::value)
    , R5 = bool(is_integral_extent<5,Args...>::value)
    , R6 = bool(is_integral_extent<6,Args...>::value)
    };

  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  typedef Kokkos::LayoutStride array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

  typedef value_type******* data_type ; 

public:

  typedef Kokkos::Experimental::ViewTraits
    < data_type
    , array_layout 
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::Experimental::View
    < data_type
    , array_layout 
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  template< class MemoryTraits >
  struct apply {

    static_assert( Kokkos::Impl::is_memory_traits< MemoryTraits >::value , "" );

    typedef Kokkos::Experimental::ViewTraits
      < data_type 
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > traits_type ;

    typedef Kokkos::Experimental::View
      < data_type 
      , array_layout
      , typename SrcTraits::device_type
      , MemoryTraits > type ;
  }; 


  typedef typename SrcTraits::dimension dimension ;

  template < class Arg0 = int, class Arg1 = int, class Arg2 = int, class Arg3 = int, class Arg4 = int, class Arg5 = int, class Arg6 = int >
  struct ExtentGenerator {
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents< 7 , rank > generator ( const dimension & dim , Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(), Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(), Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6() )
    {
       return SubviewExtents< 7 , rank>( dim , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 );
    }
  };


  typedef DynRankView< value_type , array_layout , typename SrcTraits::device_type , typename SrcTraits::memory_traits >  ret_type;

  template < typename T , class ... P >
  KOKKOS_INLINE_FUNCTION
  static ret_type subview( const unsigned src_rank , Kokkos::Experimental::View< T******* , P...> const & src 
                    , Args ... args )
    {

      typedef ViewMapping< traits_type, void >  DstType ;

       typedef typename std::conditional< (rank==0) , ViewDimension<>
                                                    , typename std::conditional< (rank==1) , ViewDimension<0>
                                                    , typename std::conditional< (rank==2) , ViewDimension<0,0>
                                                    , typename std::conditional< (rank==3) , ViewDimension<0,0,0>
                                                    , typename std::conditional< (rank==4) , ViewDimension<0,0,0,0>
                                                    , typename std::conditional< (rank==5) , ViewDimension<0,0,0,0,0>
                                                    , typename std::conditional< (rank==6) , ViewDimension<0,0,0,0,0,0>
                                                                                           , ViewDimension<0,0,0,0,0,0,0>
                                                    >::type >::type >::type >::type >::type >::type >::type  DstDimType ;

      typedef ViewOffset< DstDimType , Kokkos::LayoutStride > dst_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      ret_type dst ;

      const SubviewExtents< 7 , rank > extents = 
        ExtentGenerator< Args ... >::generator( src.m_map.m_offset.m_dim , args... ) ; 

      dst_offset_type tempdst( src.m_map.m_offset , extents ) ;

      dst.m_track = src.m_track ;

      dst.m_map.m_offset.m_dim.N0 = tempdst.m_dim.N0 ;
      dst.m_map.m_offset.m_dim.N1 = tempdst.m_dim.N1 ;
      dst.m_map.m_offset.m_dim.N2 = tempdst.m_dim.N2 ;
      dst.m_map.m_offset.m_dim.N3 = tempdst.m_dim.N3 ;
      dst.m_map.m_offset.m_dim.N4 = tempdst.m_dim.N4 ;
      dst.m_map.m_offset.m_dim.N5 = tempdst.m_dim.N5 ;
      dst.m_map.m_offset.m_dim.N6 = tempdst.m_dim.N6 ;

      dst.m_map.m_offset.m_stride.S0 = tempdst.m_stride.S0 ;
      dst.m_map.m_offset.m_stride.S1 = tempdst.m_stride.S1 ;
      dst.m_map.m_offset.m_stride.S2 = tempdst.m_stride.S2 ;
      dst.m_map.m_offset.m_stride.S3 = tempdst.m_stride.S3 ;
      dst.m_map.m_offset.m_stride.S4 = tempdst.m_stride.S4 ;
      dst.m_map.m_offset.m_stride.S5 = tempdst.m_stride.S5 ;
      dst.m_map.m_offset.m_stride.S6 = tempdst.m_stride.S6 ;

      dst.m_map.m_handle = dst_handle_type( src.m_map.m_handle +
                                      src.m_map.m_offset( extents.domain_offset(0)
                                                  , extents.domain_offset(1)
                                                  , extents.domain_offset(2)
                                                  , extents.domain_offset(3)
                                                  , extents.domain_offset(4)
                                                  , extents.domain_offset(5)
                                                  , extents.domain_offset(6)
                                                  ) );

      dst.m_rank = ( src_rank > 0 ? unsigned(R0) : 0 )
                 + ( src_rank > 1 ? unsigned(R1) : 0 )
                 + ( src_rank > 2 ? unsigned(R2) : 0 )
                 + ( src_rank > 3 ? unsigned(R3) : 0 )
                 + ( src_rank > 4 ? unsigned(R4) : 0 )
                 + ( src_rank > 5 ? unsigned(R5) : 0 )
                 + ( src_rank > 6 ? unsigned(R6) : 0 ) ;

      return dst ;
    }
};

} // end Impl


template< class V , class ... Args >
using Subdynrankview = typename Kokkos::Experimental::Impl::ViewMapping< Kokkos::Experimental::Impl::DynRankSubviewTag , V , Args... >::ret_type ;

template< class D , class ... P , class ...Args >
KOKKOS_INLINE_FUNCTION
Subdynrankview< ViewTraits<D******* , P...> , Args... > 
subdynrankview( const Kokkos::Experimental::DynRankView< D , P... > &src , Args...args)
  {
    if ( src.rank() > sizeof...(Args) ) //allow sizeof...(Args) >= src.rank(), ignore the remaining args
      { Kokkos::abort("subdynrankview: num of args must be >= rank of the source DynRankView"); }
  
    typedef Kokkos::Experimental::Impl::ViewMapping< Kokkos::Experimental::Impl::DynRankSubviewTag , Kokkos::Experimental::ViewTraits< D*******, P... > , Args... > metafcn ;

    return metafcn::subview( src.rank() , src.ConstDownCast() , args... );
  }

} // namespace Experimental
} // namespace Kokkos


namespace Kokkos {
namespace Experimental {

// overload == and !=
template< class LT , class ... LP , class RT , class ... RP >
KOKKOS_INLINE_FUNCTION
bool operator == ( const DynRankView<LT,LP...> & lhs ,
                   const DynRankView<RT,RP...> & rhs )
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
    lhs.rank()       ==  rhs.rank() &&
    lhs.data()       == rhs.data() &&
    lhs.span()       == rhs.span() &&
    lhs.dimension(0) == rhs.dimension(0) &&
    lhs.dimension(1) == rhs.dimension(1) &&
    lhs.dimension(2) == rhs.dimension(2) &&
    lhs.dimension(3) == rhs.dimension(3) &&
    lhs.dimension(4) == rhs.dimension(4) &&
    lhs.dimension(5) == rhs.dimension(5) &&
    lhs.dimension(6) == rhs.dimension(6) &&
    lhs.dimension(7) == rhs.dimension(7);
}

template< class LT , class ... LP , class RT , class ... RP >
KOKKOS_INLINE_FUNCTION
bool operator != ( const DynRankView<LT,LP...> & lhs ,
                   const DynRankView<RT,RP...> & rhs )
{
  return ! ( operator==(lhs,rhs) );
}

} //end Experimental
} //end Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/** \brief  Deep copy a value from Host memory into a view.  */
template< class DT , class ... DP >
inline
void deep_copy
  ( const DynRankView<DT,DP...> & dst
  , typename ViewTraits<DT,DP...>::const_value_type & value )
{
  deep_copy( dst.ConstDownCast() , value );
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template< class ST , class ... SP >
inline
void deep_copy
  ( typename ViewTraits<ST,SP...>::non_const_value_type & dst
  , const DynRankView<ST,SP...> & src )
{
  deep_copy( dst , src.ConstDownCast() );
}


//----------------------------------------------------------------------------
/** \brief  A deep copy between views of compatible type */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy
  ( const DynRankView<DT,DP...> & dst
  , const DynRankView<ST,SP...> & src )
{
  deep_copy( dst.ConstDownCast() , src.ConstDownCast() );
}

} //end Experimental
} //end Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

template< class T , class ... P >
inline
typename DynRankView<T,P...>::HostMirror
create_mirror( const DynRankView<T,P...> & src
             , typename std::enable_if<
                 ! std::is_same< typename Kokkos::Experimental::ViewTraits<T,P...>::array_layout
                               , Kokkos::LayoutStride >::value
               >::type * = 0
             )
{
  typedef DynRankView<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  return dst_type( std::string( src.label() ).append("_mirror")
                 , src.dimension(0)
                 , src.dimension(1)
                 , src.dimension(2)
                 , src.dimension(3)
                 , src.dimension(4)
                 , src.dimension(5)
                 , src.dimension(6)
                 , src.dimension(7) );
}


template< class T , class ... P >
inline
typename DynRankView<T,P...>::HostMirror
create_mirror( const DynRankView<T,P...> & src
             , typename std::enable_if<
                 std::is_same< typename Kokkos::Experimental::ViewTraits<T,P...>::array_layout
                             , Kokkos::LayoutStride >::value
               >::type * = 0
             )
{
  typedef DynRankView<T,P...>                   src_type ;
  typedef typename src_type::HostMirror  dst_type ;

  Kokkos::LayoutStride layout ;

  layout.dimension[0] = src.dimension(0);
  layout.dimension[1] = src.dimension(1);
  layout.dimension[2] = src.dimension(2);
  layout.dimension[3] = src.dimension(3);
  layout.dimension[4] = src.dimension(4);
  layout.dimension[5] = src.dimension(5);
  layout.dimension[6] = src.dimension(6);
  layout.dimension[7] = src.dimension(7);

  layout.stride[0] = src.stride(0);
  layout.stride[1] = src.stride(1);
  layout.stride[2] = src.stride(2);
  layout.stride[3] = src.stride(3);
  layout.stride[4] = src.stride(4);
  layout.stride[5] = src.stride(5);
  layout.stride[6] = src.stride(6);
  layout.stride[7] = src.stride(7);

  return dst_type( std::string( src.label() ).append("_mirror") , layout );
}

template< class T , class ... P >
inline
typename DynRankView<T,P...>::HostMirror
create_mirror_view( const DynRankView<T,P...> & src
                  , typename std::enable_if<(
                      std::is_same< typename DynRankView<T,P...>::memory_space
                                  , typename DynRankView<T,P...>::HostMirror::memory_space
                                  >::value
                      &&
                      std::is_same< typename DynRankView<T,P...>::data_type
                                  , typename DynRankView<T,P...>::HostMirror::data_type
                                  >::value
                    )>::type * = 0
                  )
{
  return src ;
}

template< class T , class ... P >
inline
typename DynRankView<T,P...>::HostMirror
create_mirror_view( const DynRankView<T,P...> & src
                  , typename std::enable_if< ! (
                      std::is_same< typename DynRankView<T,P...>::memory_space
                                  , typename DynRankView<T,P...>::HostMirror::memory_space
                                  >::value
                      &&
                      std::is_same< typename DynRankView<T,P...>::data_type
                                  , typename DynRankView<T,P...>::HostMirror::data_type
                                  >::value
                    )>::type * = 0
                  )
{
  return Kokkos::Experimental::create_mirror( src ); 
}

} //end Experimental
} //end Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/** \brief  Resize a view with copying old data to new data at the corresponding indices. */
template< class T , class ... P >
inline
void resize( DynRankView<T,P...> & v ,
             const size_t n0 = 0 ,
             const size_t n1 = 0 ,
             const size_t n2 = 0 ,
             const size_t n3 = 0 ,
             const size_t n4 = 0 ,
             const size_t n5 = 0 ,
             const size_t n6 = 0 ,
             const size_t n7 = 0 )
{
  typedef DynRankView<T,P...>  drview_type ;

  static_assert( Kokkos::Experimental::ViewTraits<T,P...>::is_managed , "Can only resize managed views" );

  drview_type v_resized( v.label(), n0, n1, n2, n3, n4, n5, n6, n7 );

  Kokkos::Experimental::Impl::ViewRemap< drview_type , drview_type >( v_resized , v );

  v = v_resized ;
}

/** \brief  Resize a view with copying old data to new data at the corresponding indices. */
template< class T , class ... P >
inline
void realloc( DynRankView<T,P...> & v ,
              const size_t n0 = 0 ,
              const size_t n1 = 0 ,
              const size_t n2 = 0 ,
              const size_t n3 = 0 ,
              const size_t n4 = 0 ,
              const size_t n5 = 0 ,
              const size_t n6 = 0 ,
              const size_t n7 = 0 )
{
  typedef DynRankView<T,P...>  view_type ;

  static_assert( Kokkos::Experimental::ViewTraits<T,P...>::is_managed , "Can only realloc managed views" );

  const std::string label = v.label();

  v = view_type(); // Deallocate first, if the only view to allocation
  v = view_type( label, n0, n1, n2, n3, n4, n5, n6, n7 );
}

} //end Experimental

} //end Kokkos


namespace Kokkos {

template< typename D , class ... P >
using DynRankView = Kokkos::Experimental::DynRankView< D , P... > ;

using Kokkos::Experimental::deep_copy ;
using Kokkos::Experimental::create_mirror ;
using Kokkos::Experimental::create_mirror_view ;
using Kokkos::Experimental::subdynrankview ;
using Kokkos::Experimental::resize ;
using Kokkos::Experimental::realloc ;

} //end Kokkos
#endif
