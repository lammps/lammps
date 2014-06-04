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
#include <impl/Kokkos_Traits.hpp>

namespace Kokkos {
#if   defined ( KOKKOS_HAVE_CUDA )
class Cuda ;
#endif
#if   defined ( KOKKOS_HAVE_OPENMP )
class OpenMP ;
#endif
#if   defined ( KOKKOS_HAVE_PTHREAD )
class Threads ;
#endif
#if   defined ( KOKKOS_HAVE_SERIAL )
class Serial ;
#endif
} // namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
  #if   defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef Cuda DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef OpenMP DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
    typedef Threads DefaultDeviceType;
  #elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) && \
       !defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
    typedef Serial DefaultDeviceType;
  #else
    #if   defined ( KOKKOS_HAVE_CUDA )
      typedef Kokkos::Cuda DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA
    #elif defined ( KOKKOS_HAVE_OPENMP )
      typedef OpenMP DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP
    #elif defined ( KOKKOS_HAVE_PTHREAD )
      typedef Threads DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS
    #else
      typedef Serial DefaultDeviceType;
      #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL
    #endif
  #endif
}
}

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Enable = void >
struct FunctorHasDeviceType : public false_type {};

template< class FunctorType >
struct FunctorHasDeviceType< FunctorType , typename
   enable_if< ! is_same<typename FunctorType::device_type,int>::value >::type >
  : public true_type {};
}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelFor
/// \brief Implementation of the ParallelFor operator that has a
///   partial specialization for the device.
///
/// This is an implementation detail of parallel_for.  Users should
/// skip this and go directly to the nonmember function parallel_for.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelFor ;

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

/// \class VectorParallel
/// \brief Request for parallel_for to attempt thread+vector parallelism.
struct VectorParallel
{
  const size_t nwork ;
  VectorParallel( const size_t n ) : nwork(n) {}
  operator size_t () const { return nwork ; }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief Execute \c functor \c work_count times in parallel.
 *
 * A "functor" is a class containing the function to execute in
 * parallel, any data needed for that execution, and a \c device_type
 * typedef.  Here is an example functor for parallel_for:
 *
 * \code
 *  class FunctorType {
 *  public:
 *    typedef  ...  device_type ;
 *    void operator() (IntType iwork) const ;
 *  };
 * \endcode
 *
 * In the above example, \c IntType is any integer type for which a
 * valid conversion from \c size_t to \c IntType exists.  Its
 * <tt>operator()</tt> method defines the operation to parallelize,
 * over the range of integer indices <tt>iwork=[0,work_count-1]</tt>.
 * This compares to a single iteration \c iwork of a \c for loop.
 */
template< class FunctorType >
inline
void parallel_for( const size_t        work_count ,
                   const FunctorType & functor ,
     typename Impl::enable_if<Impl::FunctorHasDeviceType<FunctorType>::value,int>::type = 0 )
{
  Impl::ParallelFor< FunctorType , size_t > tmp( functor , work_count );
}

template< class FunctorType >
inline
void parallel_for( const size_t        work_count ,
                   const FunctorType & functor ,
                   typename Impl::enable_if<!Impl::FunctorHasDeviceType<FunctorType>::value,int>::type = 0 )
{
  Impl::ParallelFor< FunctorType , size_t, Impl::DefaultDeviceType >
    tmp( functor , work_count );
}

/** \brief Execute \c functor \c work_count times in parallel, with vectorization.
 *
 * This is like parallel_for, except that it <i>mandates</i>
 * vectorization as well as parallelization of the given functor.  We
 * emphasize "mandates": this means that the user asserts that
 * vectorization is correct, and insists that the compiler vectorize.
 * Mandating vectorization is not always desirable, for example if the
 * body of the functor is complicated.  In some cases, users might
 * want to parallelize over threads, and use vectorization inside the
 * parallel operation.  Furthermore, the compiler might still be able
 * to vectorize through a parallel_for.  Thus, users should take care
 * not to use this execution option arbitrarily.
 */
template< class FunctorType >
inline
void vector_parallel_for( const size_t        work_count ,
                          const FunctorType & functor )
{
  Impl::ParallelFor< FunctorType , VectorParallel > tmp( functor , work_count );
}

template< class DeviceType >
class MultiFunctorParallelFor ;

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelReduce ;

/// \class ReduceAdapter
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class ValueType = typename FunctorType::value_type >
struct ReduceAdapter ;

} // namespace Impl
} // namespace Kokkos


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
template< class FunctorType >
inline
void parallel_reduce( const size_t        work_count ,
                      const FunctorType & functor )
{
  Impl::ParallelReduce< FunctorType , size_t > reduce( functor , work_count );
}

/** \brief  Parallel reduction and output to host.
 *
 *  If FunctorType::value_type is
 *    - \c PodType,  then \c reference_type is <tt>PodType & </tt>.
 *    - <tt>PodType[]</tt>, then \c reference_type is <tt>PodType * </tt>.
 */
template< class FunctorType >
inline
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType, size_t >
    reduce( functor , work_count , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait();
}

template< class FunctorType >
inline
void parallel_reduce( const VectorParallel & work_count ,
                      const FunctorType & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType, VectorParallel >
    reduce( functor , work_count , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait();
}

template< class DeviceType >
class MultiFunctorParallelReduce ;

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelScan
/// \brief Implementation detail of parallel_scan.
///
/// This is an implementation detail of parallel_scan.  Users should
/// skip this and go directly to the documentation of the nonmember
/// template function Kokkos::parallel_scan.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelScan ;

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

/// \fn parallel_scan
/// \tparam FunctorType Type of the scan functor.
///
/// \param work_count [in] Number of work items.
/// \param functor [in] The scan functor.
///
/// This function implements a parallel scan operation.  The scan can
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
/// class ScanFunctor {
/// public:
///   // The Kokkos device type
///   typedef ... device_type;
///   // Type of an entry of the array containing the result;
///   // also the type of each of the entries combined using
///   // operator() or join().
///   typedef PodType value_type;
///   typedef typename DeviceType::size_type size_type;
///
///   void operator () (const size_type i, value_type& update, const bool final_pass) const;
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
///   ExclScanFunctor (Kokkos::View<value_type*, device_type> x) :
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
template< class FunctorType >
inline
void parallel_scan( const size_t        work_count ,
                    const FunctorType & functor )
{
  Impl::ParallelScan< FunctorType , size_t > scan( functor , work_count );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Parallel work request for shared memory, league size, and team size.
 *
 *  If the shared size is too large then slow (global) memory will be used.
 *  If the league or team size are too large then they will be reduced.
 */
struct ParallelWorkRequest {
  size_t  league_size ; ///<  Size of league (number of teams in a league)
  size_t  team_size ;   ///<  Size of team (number of threads in a team)

  KOKKOS_INLINE_FUNCTION
  ParallelWorkRequest() : league_size(0), team_size(0) {}

  KOKKOS_INLINE_FUNCTION
  ParallelWorkRequest( size_t s0 , size_t s1 ) : league_size(s0), team_size(s1) {}
};

/** \brief  Execute functor in parallel with work request,
 *          the actual league_size and team_size may be smaller.
 *
 *  class FunctorType {
 *  public:
 *    typedef  ...  device_type ;
 *    void operator()( device_type ) const ;
 *  };
 */
template< class FunctorType >
inline
void parallel_for( const ParallelWorkRequest & request ,
                   const FunctorType         & functor )
{
  Kokkos::Impl::ParallelFor< FunctorType , ParallelWorkRequest >( functor , request );
}

} // namespace Kokkos

namespace Kokkos {

/** \brief  Parallel reduction.
 *
 *  class FunctorType {
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type ; // POD type
 *    void operator()( device_type , <podType> & ) const ;
 *    void init( <podType> & ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 *
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( device_type , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 */
template< class FunctorType >
inline
void parallel_reduce( const Kokkos::ParallelWorkRequest  & request ,
                      const FunctorType          & functor )
{
  Impl::ParallelReduce< FunctorType , Kokkos::ParallelWorkRequest > reduce( functor , request );
}

template< class FunctorType >
inline
void parallel_reduce( const Kokkos::ParallelWorkRequest  & request ,
                      const FunctorType          & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType , Kokkos::ParallelWorkRequest >
    reduce( functor , request , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait(); // Wait for reduce to complete and output result
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

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
struct FunctorShmemSize
{
  static inline size_t value( const FunctorType & ) { return 0 ; }
};

template< class FunctorType >
struct FunctorShmemSize< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::shmem_size ) >::type >
{
  static inline size_t value( const FunctorType & f ) { return f.shmem_size() ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

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

  KOKKOS_INLINE_FUNCTION static
  void join( const FunctorType & f , volatile void * update , volatile const void * input )
    { f.join( *((volatile ScalarType*)update) , *((volatile const ScalarType*)input) ); }

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

  KOKKOS_INLINE_FUNCTION static
  void join( const FunctorType & f , volatile void * update , volatile const void * input )
    { f.join( ((volatile ScalarType*)update) , ((volatile const ScalarType*)input) ); }

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

