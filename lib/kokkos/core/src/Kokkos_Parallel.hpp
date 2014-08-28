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

/// \file Kokkos_Parallel.hpp
/// \brief Declaration of parallel operators

#ifndef KOKKOS_PARALLEL_HPP
#define KOKKOS_PARALLEL_HPP

#include <cstddef>
#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Given a Functor and Execution Policy query an execution space.
 *
 *  if       the Policy has an execution space use that
 *  else if  the Functor has a device_type use that
 *  else     use the default
 */
template< class Functor
        , class Policy
        , class EnableFunctor = void
        , class EnablePolicy  = void
        >
struct FunctorPolicyExecutionSpace {
  typedef Kokkos::DefaultExecutionSpace execution_space ;
};

template< class Functor , class Policy >
struct FunctorPolicyExecutionSpace
  < Functor , Policy
  , typename enable_if_type< typename Functor::device_type     >::type
  , typename enable_if_type< typename Policy ::execution_space >::type
  >
{
  typedef typename Policy ::execution_space execution_space ;
};

template< class Functor , class Policy , class EnableFunctor >
struct FunctorPolicyExecutionSpace
  < Functor , Policy
  , EnableFunctor
  , typename enable_if_type< typename Policy::execution_space >::type
  >
{
  typedef typename Policy ::execution_space execution_space ;
};

template< class Functor , class Policy , class EnablePolicy >
struct FunctorPolicyExecutionSpace
  < Functor , Policy
  , typename enable_if_type< typename Functor::device_type >::type
  , EnablePolicy
  >
{
  typedef typename Functor::device_type execution_space ;
};

//----------------------------------------------------------------------------

template< class FunctorType , class Enable = void >
struct ReduceAdapterValueType ;

template< class FunctorType >
struct ReduceAdapterValueType< FunctorType , typename enable_if_type< typename FunctorType::value_type >::type >
{
  typedef typename FunctorType::value_type type ;
};

/// \class ReduceAdapter
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class ValueType = typename ReduceAdapterValueType< FunctorType >::type >
struct ReduceAdapter ;

//----------------------------------------------------------------------------
/// \class ParallelFor
/// \brief Implementation of the ParallelFor operator that has a
///   partial specialization for the device.
///
/// This is an implementation detail of parallel_for.  Users should
/// skip this and go directly to the nonmember function parallel_for.
template< class FunctorType 
        , class ExecPolicy 
        , class ExecSpace = typename FunctorPolicyExecutionSpace< FunctorType , ExecPolicy >::execution_space
        >
class ParallelFor ;

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType 
        , class ExecPolicy 
        , class ExecSpace = typename FunctorPolicyExecutionSpace< FunctorType , ExecPolicy >::execution_space
        >
class ParallelReduce ;

/// \class ParallelScan
/// \brief Implementation detail of parallel_scan.
///
/// This is an implementation detail of parallel_scan.  Users should
/// skip this and go directly to the documentation of the nonmember
/// template function Kokkos::parallel_scan.
template< class FunctorType 
        , class ExecPolicy 
        , class ExecSpace = typename FunctorPolicyExecutionSpace< FunctorType , ExecPolicy >::execution_space
        >
class ParallelScan ;

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief Execute \c functor in parallel according to the execution \c policy.
 *
 * A "functor" is a class containing the function to execute in parallel,
 * data needed for that execution, and an optional \c device_type
 * typedef.  Here is an example functor for parallel_for:
 *
 * \code
 *  class FunctorType {
 *  public:
 *    typedef  ...  device_type ;
 *    void operator() ( WorkType iwork ) const ;
 *  };
 * \endcode
 *
 * In the above example, \c WorkType is any integer type for which a
 * valid conversion from \c size_t to \c IntType exists.  Its
 * <tt>operator()</tt> method defines the operation to parallelize,
 * over the range of integer indices <tt>iwork=[0,work_count-1]</tt>.
 * This compares to a single iteration \c iwork of a \c for loop.
 * If \c device_type is not defined DefaultExecutionSpace will be used.
 */
template< class ExecPolicy , class FunctorType >
inline
void parallel_for( const ExecPolicy  & policy
                 , const FunctorType & functor
                 , typename Impl::enable_if< ! Impl::is_integral< ExecPolicy >::value >::type * = 0
                 )
{
  (void) Impl::ParallelFor< FunctorType , ExecPolicy >( functor , policy );
}

template< class FunctorType >
inline
void parallel_for( const size_t        work_count ,
                   const FunctorType & functor )
{
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;
  typedef RangePolicy< execution_space > policy ;
  (void) Impl::ParallelFor< FunctorType , policy >( functor , policy(0,work_count) );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Parallel reduction
 *
 * Example of a parallel_reduce functor for a POD (plain old data) value type:
 * \code
 *  class FunctorType { // For POD value type
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type ;
 *    void operator()( <intType> iwork , <podType> & update ) const ;
 *    void init( <podType> & update ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 * \endcode
 *
 * Example of a parallel_reduce functor for an array of POD (plain old data) values:
 * \code
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( <intType> , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */
template< class ExecPolicy , class FunctorType >
inline
void parallel_reduce( const ExecPolicy  & policy 
                    , const FunctorType & functor
                    , typename Impl::enable_if< ! Impl::is_integral< ExecPolicy >::value >::type * = 0
                    )
{
  (void) Impl::ParallelReduce< FunctorType , ExecPolicy >( functor , policy );
}

// integral range policy
template< class FunctorType >
inline
void parallel_reduce( const size_t        work_count
                    , const FunctorType & functor
                    )
{
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef RangePolicy< execution_space > policy ;

  typedef Kokkos::Impl::ReduceAdapter< FunctorType >  Reduce ;

  typedef typename Kokkos::Impl::if_c< Reduce::StaticValueSize
                                     , typename Reduce::scalar_type
                                     , typename Reduce::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , typename execution_space::host_mirror_device_type
              , Kokkos::MemoryUnmanaged
              >
    result_view ;

  (void) Impl::ParallelReduce< FunctorType , policy >( functor , policy(0,work_count) , result_view );
}

// general policy and view ouput
template< class ExecPolicy , class FunctorType , class ViewType >
inline
void parallel_reduce( const ExecPolicy  & policy 
                    , const FunctorType & functor 
                    , const ViewType    & result_view
                    , typename Impl::enable_if<
                      ( Impl::is_view<ViewType>::value && ! Impl::is_integral< ExecPolicy >::value
                      )>::type * = 0 )
{
  (void) Impl::ParallelReduce< FunctorType, ExecPolicy >( functor , policy , result_view );
}

// general policy and pod or array of pod output
template< class ExecPolicy , class FunctorType >
inline
void parallel_reduce( const ExecPolicy  & policy 
                    , const FunctorType & functor 
                    , typename Impl::enable_if<
                      ( ! Impl::is_integral< ExecPolicy >::value )
                      , typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type
                      >::type result_ref )
{
  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace< FunctorType , ExecPolicy >::execution_space
      execution_space ;

  typedef Kokkos::Impl::ReduceAdapter< FunctorType >  Reduce ;

  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< Reduce::StaticValueSize
                                     , typename Reduce::scalar_type
                                     , typename Reduce::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , typename execution_space::host_mirror_device_type
              , Kokkos::MemoryUnmanaged
              >
    result_view( Reduce::pointer( result_ref )
               , Reduce::value_count( functor )
               );

  (void) Impl::ParallelReduce< FunctorType, ExecPolicy >( functor , policy , result_view );
}

// integral range policy and view ouput
template< class FunctorType , class ViewType >
inline
void parallel_reduce( const size_t        work_count
                    , const FunctorType & functor 
                    , const ViewType    & result_view
                    , typename Impl::enable_if<( Impl::is_view<ViewType>::value )>::type * = 0 )
{
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef RangePolicy< execution_space > ExecPolicy ;

  (void) Impl::ParallelReduce< FunctorType, ExecPolicy >( functor , ExecPolicy(0,work_count) , result_view );
}

// integral range policy and pod or array of pod output
template< class FunctorType >
inline
void parallel_reduce( const size_t        work_count ,
                      const FunctorType & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef Kokkos::RangePolicy< execution_space > policy ;

  typedef Kokkos::Impl::ReduceAdapter< FunctorType >  Reduce ;

  // Wrap the result output request in a view to inform the implementation
  // of the type and memory space.

  typedef typename Kokkos::Impl::if_c< Reduce::StaticValueSize
                                     , typename Reduce::scalar_type
                                     , typename Reduce::pointer_type
                                     >::type value_type ;

  Kokkos::View< value_type
              , typename execution_space::host_mirror_device_type
              , Kokkos::MemoryUnmanaged
              >
    result_view( Reduce::pointer( result )
               , Reduce::value_count( functor )
               );
  
  (void) Impl::ParallelReduce< FunctorType , policy >( functor , policy(0,work_count) , result_view );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/// \fn parallel_scan
/// \tparam ExecutionPolicy The execution policy type.
/// \tparam FunctorType     The scan functor type.
///
/// \param policy  [in] The execution policy.
/// \param functor [in] The scan functor.
///
/// This function implements a parallel scan pattern.  The scan can
/// be either inclusive or exclusive, depending on how you implement
/// the scan functor.
///
/// A scan functor looks almost exactly like a reduce functor, except
/// that its operator() takes a third \c bool argument, \c final_pass,
/// which indicates whether this is the last pass of the scan
/// operation.  We will show below how to use the \c final_pass
/// argument to control whether the scan is inclusive or exclusive.
///
/// Here is the minimum required interface of a scan functor for a POD
/// (plain old data) value type \c PodType.  That is, the result is a
/// View of zero or more PodType.  It is also possible for the result
/// to be an array of (same-sized) arrays of PodType, but we do not
/// show the required interface for that here.
/// \code
/// template< class ExecPolicy , class FunctorType >
/// class ScanFunctor {
/// public:
///   // The Kokkos device type
///   typedef ... device_type;
///   // Type of an entry of the array containing the result;
///   // also the type of each of the entries combined using
///   // operator() or join().
///   typedef PodType value_type;
///
///   void operator () (const ExecPolicy::member_type & i, value_type& update, const bool final_pass) const;
///   void init (value_type& update) const;
///   void join (volatile value_type& update, volatile const value_type& input) const
/// };
/// \endcode
///
/// Here is an example of a functor which computes an inclusive plus-scan
/// of an array of \c int, in place.  If given an array [1, 2, 3, 4], this
/// scan will overwrite that array with [1, 3, 6, 10].
///
/// \code
/// template<class DeviceType>
/// class InclScanFunctor {
/// public:
///   typedef DeviceType device_type;
///   typedef int value_type;
///   typedef typename DeviceType::size_type size_type;
///
///   InclScanFunctor (Kokkos::View<value_type*, device_type> x) : x_ (x) {}
///
///   void operator () (const size_type i, value_type& update, const bool final_pass) const {
///     update += x_(i);
///     if (final_pass) {
///       x_(i) = update;
///     }
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input) const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, device_type> x_;
/// };
/// \endcode
///
/// Here is an example of a functor which computes an <i>exclusive</i>
/// scan of an array of \c int, in place.  In operator(), note both
/// that the final_pass test and the update have switched places, and
/// the use of a temporary.  If given an array [1, 2, 3, 4], this scan
/// will overwrite that array with [0, 1, 3, 6].
///
/// \code
/// template<class DeviceType>
/// class ExclScanFunctor {
/// public:
///   typedef DeviceType device_type;
///   typedef int value_type;
///   typedef typename DeviceType::size_type size_type;
///
///   ExclScanFunctor (Kokkos::View<value_type*, device_type> x) : x_ (x) {}
///
///   void operator () (const size_type i, value_type& update, const bool final_pass) const {
///     const value_type x_i = x_(i);
///     if (final_pass) {
///       x_(i) = update;
///     }
///     update += x_i;
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input) const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, device_type> x_;
/// };
/// \endcode
///
/// Here is an example of a functor which builds on the above
/// exclusive scan example, to compute an offsets array from a
/// population count array, in place.  We assume that the pop count
/// array has an extra entry at the end to store the final count.  If
/// given an array [1, 2, 3, 4, 0], this scan will overwrite that
/// array with [0, 1, 3, 6, 10].
///
/// \code
/// template<class DeviceType>
/// class OffsetScanFunctor {
/// public:
///   typedef DeviceType device_type;
///   typedef int value_type;
///   typedef typename DeviceType::size_type size_type;
///
///   // lastIndex_ is the last valid index (zero-based) of x.
///   // If x has length zero, then lastIndex_ won't be used anyway.
///   OffsetScanFunctor (Kokkos::View<value_type*, device_type> x) :
///     x_ (x), last_index_ (x.dimension_0 () == 0 ? 0 : x.dimension_0 () - 1)
///   {}
///
///   void operator () (const size_type i, int& update, const bool final_pass) const {
///     const value_type x_i = x_(i);
///     if (final_pass) {
///       x_(i) = update;
///     }
///     update += x_i;
///     // The last entry of x_ gets the final sum.
///     if (final_pass && i == last_index_) {
///       x_(i) = update;
///     }
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input) const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, device_type> x_;
///   const size_type last_index_;
/// };
/// \endcode
///
template< class ExecutionPolicy , class FunctorType >
inline
void parallel_scan( const ExecutionPolicy & policy 
                  , const FunctorType     & functor 
                  , typename Impl::enable_if< ! Impl::is_integral< ExecutionPolicy >::value >::type * = 0
                  )
{
  Impl::ParallelScan< FunctorType , ExecutionPolicy > scan( functor , policy );
}

template< class FunctorType >
inline
void parallel_scan( const size_t        work_count ,
                    const FunctorType & functor )
{
  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef Kokkos::RangePolicy< execution_space > policy ;

  (void) Impl::ParallelScan< FunctorType , policy >( functor , policy(0,work_count) );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Enable = void >
struct FunctorHasInit : public false_type {};

template< class FunctorType >
struct FunctorHasInit< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::init ) >::type >
  : public true_type {};

template< class FunctorType , class Enable = void >
struct FunctorHasJoin : public false_type {};

template< class FunctorType >
struct FunctorHasJoin< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::join ) >::type >
  : public true_type {};

template< class FunctorType , class Enable = void >
struct FunctorHasFinal : public false_type {};

template< class FunctorType >
struct FunctorHasFinal< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::final ) >::type >
  : public true_type {};

template< class FunctorType , class Enable = void >
struct FunctorTeamShmemSize
{
  static inline size_t value( const FunctorType & , int ) { return 0 ; }
};

template< class FunctorType >
struct FunctorTeamShmemSize< FunctorType , typename enable_if< sizeof( & FunctorType::team_shmem_size ) >::type >
{
  static inline size_t value( const FunctorType & f , int team_size ) { return f.team_shmem_size( team_size ) ; }
};

template< class FunctorType >
struct FunctorTeamShmemSize< FunctorType , typename enable_if< sizeof( & FunctorType::shmem_size ) >::type >
{
  static inline size_t value( const FunctorType & f , int team_size ) { return f.shmem_size( team_size ) ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#if defined( KOKKOS_HAVE_CXX11 )

template< class FunctionPtr >
struct ReduceAdapterFunctorOperatorArgType ;

template< class Functor , class Arg0 , class Arg1 >
struct ReduceAdapterFunctorOperatorArgType< void ( Functor::*)( Arg0 , Arg1 & ) const > {
  typedef Arg1 type ;
};

// Functor does not have a 'typedef ... value_type' and C++11 is enabled.
// Deduce the value type from the functor's argument list.
template< class FunctorType , class Enable >
struct ReduceAdapterValueType {
private:
  typedef decltype( & FunctorType::operator() ) function_pointer_type ;
public:
  typedef typename ReduceAdapterFunctorOperatorArgType< function_pointer_type >::type type ;
};

#endif

template< class FunctorType , class ScalarType >
struct ReduceAdapter
{
  enum { StaticValueSize = sizeof(ScalarType) };

  typedef ScalarType & reference_type  ;
  typedef ScalarType * pointer_type  ;
  typedef ScalarType   scalar_type  ;

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p ) { return *((ScalarType*) p); }

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p , unsigned i ) { return ((ScalarType*) p)[i]; }

  KOKKOS_INLINE_FUNCTION static
  pointer_type pointer( reference_type p ) { return & p ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_count( const FunctorType & ) { return 1 ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_size( const FunctorType & ) { return sizeof(ScalarType); }

  KOKKOS_INLINE_FUNCTION static
  void copy( const FunctorType & , void * const dst , const void * const src )
    { *((scalar_type*)dst) = *((const scalar_type*)src); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void join( const F & f
           , volatile void * update
           , typename enable_if< is_same<F,FunctorType>::value &&
                                 FunctorHasJoin<F>::value
                               , volatile const void *
                               >::type input )
    { f.join( *((volatile ScalarType*)update) , *((volatile const ScalarType*)input) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void join( const F & f
           , volatile void * update
           , typename enable_if< is_same<F,FunctorType>::value &&
                                 ! FunctorHasJoin<F>::value
                               , volatile const void *
                               >::type input )
    { *((volatile ScalarType*)update) += *((volatile const ScalarType*)input); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  reference_type
  init( const F & f ,
        typename enable_if< ( is_same<F,FunctorType>::value &&
                              FunctorHasInit<F>::value )
                          >::type * p )
    { f.init( *((ScalarType *) p ) );  return *((ScalarType *) p ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  reference_type
  init( const F & ,
        typename enable_if< ( is_same<F,FunctorType>::value &&
                              ! FunctorHasInit<F>::value )
                          >::type * p )
    { return *( new(p) ScalarType() ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & f ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    FunctorHasFinal<F>::value )
                                >::type * p )
    { f.final( *((ScalarType *) p ) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    ! FunctorHasFinal<F>::value )
                                >::type * )
    {}
};

template< class FunctorType , class ScalarType >
struct ReduceAdapter< FunctorType , ScalarType[] >
{
  enum { StaticValueSize = 0 };

  typedef ScalarType * reference_type  ;
  typedef ScalarType * pointer_type  ;
  typedef ScalarType   scalar_type  ;

  KOKKOS_INLINE_FUNCTION static
  ScalarType * reference( void * p ) { return (ScalarType*) p ; }

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p , unsigned i ) { return ((ScalarType*) p)+i; }

  KOKKOS_INLINE_FUNCTION static
  pointer_type pointer( reference_type p ) { return p ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_count( const FunctorType & f ) { return f.value_count ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_size( const FunctorType & f ) { return f.value_count * sizeof(ScalarType); }

  KOKKOS_INLINE_FUNCTION static
  void copy( const FunctorType & f , void * const dst , const void * const src )
    {
      for ( int i = 0 ; i < int(f.value_count) ; ++i ) {
        ((scalar_type*)dst)[i] = ((const scalar_type*)src)[i];
      }
    }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void join( const F & f
           , volatile void * const update
           , typename enable_if< is_same<F,FunctorType>::value &&
                                 FunctorHasJoin<F>::value
                               , volatile const void * const
                               >::type input )
    { f.join( ((volatile ScalarType*)update) , ((volatile const ScalarType*)input) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void join( const F & f
           , volatile void * const update
           , typename enable_if< is_same<F,FunctorType>::value &&
                                 ! FunctorHasJoin<F>::value
                               , volatile const void * const
                               >::type input )
    {
      for ( int i = 0 ; i < int(f.value_count) ; ++i ) {
        ((volatile ScalarType*)update)[i] += ((volatile const ScalarType*)input)[i] ;
      }
    }


  template< class F >
  KOKKOS_INLINE_FUNCTION static
  reference_type
  init( const F & f ,
        typename enable_if< ( is_same<F,FunctorType>::value &&
                              FunctorHasInit<F>::value )
                          >::type * p )
    { f.init( ((ScalarType *) p ) ); return (ScalarType*) p ; }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  reference_type
  init( const F & f ,
        typename enable_if< ( is_same<F,FunctorType>::value &&
                              ! FunctorHasInit<F>::value )
                          >::type * p )
    {
      for ( int i = 0 ; i < int(f.value_count) ; ++i ) {
        new(((ScalarType*)p)+i) ScalarType();
      }
      return (ScalarType*)p ;
    }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & f ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    FunctorHasFinal<F>::value )
                                >::type * p )
    { f.final( ((ScalarType *) p ) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    ! FunctorHasFinal<F>::value )
                                >::type * )
    {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_PARALLEL_HPP */

