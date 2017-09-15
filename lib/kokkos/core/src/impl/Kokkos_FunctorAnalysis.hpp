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

#ifndef KOKKOS_FUNCTORANALYSIS_HPP
#define KOKKOS_FUNCTORANALYSIS_HPP

#include <cstddef>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct FunctorPatternInterface {
  struct FOR {};
  struct REDUCE {};
  struct SCAN {};
};

/** \brief  Query Functor and execution policy argument tag for value type.
 *
 *  If 'value_type' is not explicitly declared in the functor
 *  then attempt to deduce the type from FunctorType::operator()
 *  interface used by the pattern and policy.
 *
 *  For the REDUCE pattern generate a Reducer and finalization function
 *  derived from what is available within the functor.
 */
template< typename PatternInterface , class Policy , class Functor >
struct FunctorAnalysis {
private:

  using FOR    = FunctorPatternInterface::FOR ;
  using REDUCE = FunctorPatternInterface::REDUCE ;
  using SCAN   = FunctorPatternInterface::SCAN ;

  //----------------------------------------

  struct VOID {};

  template< typename P = Policy , typename = std::false_type >
  struct has_work_tag
    {
      using type = void ;
      using wtag = VOID ;
    };

  template< typename P >
  struct has_work_tag
    < P , typename std::is_same< typename P::work_tag , void >::type >
    {
      using type = typename P::work_tag ;
      using wtag = typename P::work_tag ;
    };

  using Tag  = typename has_work_tag<>::type ;
  using WTag = typename has_work_tag<>::wtag ;

  //----------------------------------------
  // Check for T::execution_space

  template< typename T , typename = std::false_type >
  struct has_execution_space { using type = void ; enum { value = false }; };

  template< typename T >
  struct has_execution_space
    < T , typename std::is_same< typename T::execution_space , void >::type >
  {
    using type = typename T::execution_space ;
    enum { value = true };
  };

  using policy_has_space  = has_execution_space< Policy > ;
  using functor_has_space = has_execution_space< Functor > ;

  static_assert( ! policy_has_space::value ||
                 ! functor_has_space::value ||
                 std::is_same< typename policy_has_space::type
                             , typename functor_has_space::type >::value
               , "Execution Policy and Functor execution space must match" );

  //----------------------------------------
  // Check for Functor::value_type, which is either a simple type T or T[]

  template< typename F , typename = std::false_type >
  struct has_value_type { using type = void ; };

  template< typename F >
  struct has_value_type
    < F , typename std::is_same< typename F::value_type , void >::type >
  {
    using type = typename F::value_type ;

    static_assert( ! std::is_reference< type >::value &&
                   std::rank< type >::value <= 1 &&
                   std::extent< type >::value == 0
                 , "Kokkos Functor::value_type is T or T[]" );
  };

  //----------------------------------------
  // If Functor::value_type does not exist then evaluate operator(),
  // depending upon the pattern and whether the policy has a work tag,
  // to determine the reduction or scan value_type.

  template< typename F
          , typename P = PatternInterface
          , typename V = typename has_value_type<F>::type
          , bool     T = std::is_same< Tag , void >::value
          >
  struct deduce_value_type { using type = V ; };

  template< typename F >
  struct deduce_value_type< F , REDUCE , void , true > {

    template< typename M , typename A >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( M , A & ) const );

    using type = decltype( deduce( & F::operator() ) );
  };

  template< typename F >
  struct deduce_value_type< F , REDUCE , void , false > {

    template< typename M , typename A >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( WTag , M , A & ) const );

    template< typename M , typename A >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( WTag const & , M , A & ) const );

    using type = decltype( deduce( & F::operator() ) );
  };

  template< typename F >
  struct deduce_value_type< F , SCAN , void , true > {

    template< typename M , typename A , typename I >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( M , A & , I ) const );

    using type = decltype( deduce( & F::operator() ) );
  };

  template< typename F >
  struct deduce_value_type< F , SCAN , void , false > {

    template< typename M , typename A , typename I >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( WTag , M , A & , I ) const );

    template< typename M , typename A , typename I >
    KOKKOS_INLINE_FUNCTION static
    A deduce( void (Functor::*)( WTag const & , M , A & , I ) const );

    using type = decltype( deduce( & F::operator() ) );
  };

  //----------------------------------------

  using candidate_type = typename deduce_value_type< Functor >::type ;

  enum { candidate_is_void  = std::is_same< candidate_type , void >::value
       , candidate_is_array = std::rank< candidate_type >::value == 1 };

  //----------------------------------------

public:

  using execution_space = typename std::conditional
    < functor_has_space::value
    , typename functor_has_space::type
    , typename std::conditional
      < policy_has_space::value
      , typename policy_has_space::type
      , Kokkos::DefaultExecutionSpace
      >::type
    >::type ;

  using value_type = typename std::remove_extent< candidate_type >::type ;

  static_assert( ! std::is_const< value_type >::value
               , "Kokkos functor operator reduce argument cannot be const" );

private:

  // Stub to avoid defining a type 'void &'
  using ValueType = typename
    std::conditional< candidate_is_void , VOID , value_type >::type ;

public:

  using pointer_type = typename
    std::conditional< candidate_is_void , void , ValueType * >::type ;

  using reference_type = typename
    std::conditional< candidate_is_array  , ValueType * , typename
    std::conditional< ! candidate_is_void , ValueType & , void >
    ::type >::type ;

private:

  template< bool IsArray , class FF >
  KOKKOS_INLINE_FUNCTION static constexpr
  typename std::enable_if< IsArray , unsigned >::type
  get_length( FF const & f ) { return f.value_count ; }

  template< bool IsArray , class FF >
  KOKKOS_INLINE_FUNCTION static constexpr
  typename std::enable_if< ! IsArray , unsigned >::type
  get_length( FF const & ) { return candidate_is_void ? 0 : 1 ; }

public:

  enum { StaticValueSize = ! candidate_is_void &&
                           ! candidate_is_array
                         ? sizeof(ValueType) : 0 };

  KOKKOS_FORCEINLINE_FUNCTION static constexpr
  unsigned value_count( const Functor & f )
    { return FunctorAnalysis::template get_length< candidate_is_array >(f); }

  KOKKOS_FORCEINLINE_FUNCTION static constexpr
  unsigned value_size( const Functor & f )
    { return FunctorAnalysis::template get_length< candidate_is_array >(f) * sizeof(ValueType); }

  //----------------------------------------

  template< class Unknown >
  KOKKOS_FORCEINLINE_FUNCTION static constexpr
  unsigned value_count( const Unknown & )
    { return candidate_is_void ? 0 : 1  ; }

  template< class Unknown >
  KOKKOS_FORCEINLINE_FUNCTION static constexpr
  unsigned value_size( const Unknown & )
    { return candidate_is_void ? 0 : sizeof(ValueType); }

private:

  enum INTERFACE : int
    { DISABLE           = 0
    , NO_TAG_NOT_ARRAY  = 1
    , NO_TAG_IS_ARRAY   = 2
    , HAS_TAG_NOT_ARRAY = 3
    , HAS_TAG_IS_ARRAY  = 4
    , DEDUCED =
       ! std::is_same< PatternInterface , REDUCE >::value ? DISABLE : (
       std::is_same<Tag,void>::value
         ? (candidate_is_array ? NO_TAG_IS_ARRAY  : NO_TAG_NOT_ARRAY)
         : (candidate_is_array ? HAS_TAG_IS_ARRAY : HAS_TAG_NOT_ARRAY) )
    };

  //----------------------------------------
  // parallel_reduce join operator

  template< class F , INTERFACE >
  struct has_join_function ;

  template< class F >
  struct has_join_function< F , NO_TAG_NOT_ARRAY >
    {
      typedef volatile       ValueType & vref_type ;
      typedef volatile const ValueType & cvref_type ;

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const f
               , ValueType volatile * dst
               , ValueType volatile const * src )
        { f->join( *dst , *src ); }
    };

  template< class F >
  struct has_join_function< F , NO_TAG_IS_ARRAY >
    {
      typedef volatile       ValueType * vref_type ;
      typedef volatile const ValueType * cvref_type ;

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const f
               , ValueType volatile * dst
               , ValueType volatile const * src )
        { f->join( dst , src ); }
    };

  template< class F >
  struct has_join_function< F , HAS_TAG_NOT_ARRAY >
    {
      typedef volatile       ValueType & vref_type ;
      typedef volatile const ValueType & cvref_type ;

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const f
               , ValueType volatile * dst
               , ValueType volatile const * src )
        { f->join( WTag() , *dst , *src ); }
    };

  template< class F >
  struct has_join_function< F , HAS_TAG_IS_ARRAY >
    {
      typedef volatile       ValueType * vref_type ;
      typedef volatile const ValueType * cvref_type ;

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , vref_type , cvref_type ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , vref_type , cvref_type ) );

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const f
               , ValueType volatile * dst
               , ValueType volatile const * src )
        { f->join( WTag() , dst , src ); }
    };


  template< class F   = Functor
          , INTERFACE = DEDUCED
          , typename  = void >
  struct DeduceJoin
    {
      enum { value = false };

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const f
               , ValueType volatile * dst
               , ValueType volatile const * src )
       {
         const int n = FunctorAnalysis::value_count( *f );
         for ( int i = 0 ; i < n ; ++i ) dst[i] += src[i];
       }
    };

  template< class F >
  struct DeduceJoin< F , DISABLE , void >
    {
      enum { value = false };

      KOKKOS_INLINE_FUNCTION static
      void join( F const * const
               , ValueType volatile *
               , ValueType volatile const * ) {}
    };

  template< class F , INTERFACE I >
  struct DeduceJoin< F , I ,
    decltype( has_join_function<F,I>::enable_if( & F::join ) ) >
    : public has_join_function<F,I>
    { enum { value = true }; };

  //----------------------------------------

  template< class , INTERFACE >
  struct has_init_function ;

  template< class F >
  struct has_init_function< F , NO_TAG_NOT_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const f , ValueType * dst )
        { f->init( *dst ); }
    };

  template< class F >
  struct has_init_function< F , NO_TAG_IS_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const f , ValueType * dst )
        { f->init( dst ); }
    };

  template< class F >
  struct has_init_function< F , HAS_TAG_NOT_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const f , ValueType * dst )
        { f->init( WTag(), *dst ); }
    };

  template< class F >
  struct has_init_function< F , HAS_TAG_IS_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const f , ValueType * dst )
        { f->init( WTag(), dst ); }
    };

  template< class F   = Functor
          , INTERFACE = DEDUCED
          , typename  = void >
  struct DeduceInit
    {
      enum { value = false };

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const , ValueType * dst ) { new(dst) ValueType(); }
    };

  template< class F >
  struct DeduceInit< F , DISABLE , void >
    {
      enum { value = false };

      KOKKOS_INLINE_FUNCTION static
      void init( F const * const , ValueType * ) {}
    };

  template< class F , INTERFACE I >
  struct DeduceInit< F , I ,
    decltype( has_init_function<F,I>::enable_if( & F::init ) ) >
    : public has_init_function<F,I>
    { enum { value = true }; };

  //----------------------------------------

  template< class , INTERFACE >
  struct has_final_function ;

  // No tag, not array
  template< class F >
  struct has_final_function< F , NO_TAG_NOT_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void final( F const * const f , ValueType * dst )
        { f->final( *dst ); }
    };

  // No tag, is array
  template< class F >
  struct has_final_function< F , NO_TAG_IS_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void final( F const * const f , ValueType * dst )
        { f->final( dst ); }
    };

  // Has tag, not array
  template< class F >
  struct has_final_function< F , HAS_TAG_NOT_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , ValueType & ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , ValueType & ) );

      KOKKOS_INLINE_FUNCTION static
      void final( F const * const f , ValueType * dst )
        { f->final( WTag(), *dst ); }
    };

  // Has tag, is array
  template< class F >
  struct has_final_function< F , HAS_TAG_IS_ARRAY >
    {
      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag , ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (F::*)( WTag const & , ValueType * ) const );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag , ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void enable_if( void (*)( WTag const & , ValueType * ) );

      KOKKOS_INLINE_FUNCTION static
      void final( F const * const f , ValueType * dst )
        { f->final( WTag(), dst ); }
    };

  template< class F   = Functor
          , INTERFACE = DEDUCED
          , typename  = void >
  struct DeduceFinal
    {
      enum { value = false };

      KOKKOS_INLINE_FUNCTION
      static void final( F const * const , ValueType * ) {}
    };

  template< class F , INTERFACE I >
  struct DeduceFinal< F , I ,
    decltype( has_final_function<F,I>::enable_if( & F::final ) ) >
    : public has_final_function<F,I>
    { enum { value = true }; };

  //----------------------------------------

  template< class F = Functor , typename = void >
  struct DeduceTeamShmem
    {
      enum { value = false };

      static size_t team_shmem_size( F const & , int ) { return 0 ; }
    };

  template< class F >
  struct DeduceTeamShmem< F , typename std::enable_if< 0 < sizeof( & F::team_shmem_size ) >::type >
    {
      enum { value = true };

      static size_t team_shmem_size( F const * const f , int team_size )
        { return f->team_shmem_size( team_size ); }
    };

  template< class F >
  struct DeduceTeamShmem< F , typename std::enable_if< 0 < sizeof( & F::shmem_size ) >::type >
    {
      enum { value = true };

      static size_t team_shmem_size( F const * const f , int team_size )
        { return f->shmem_size( team_size ); }
    };

  //----------------------------------------

public:

  inline static
  size_t team_shmem_size( Functor const & f )
    { return DeduceTeamShmem<>::team_shmem_size( f ); }

  //----------------------------------------

  enum { has_join_member_function  = DeduceJoin<>::value };
  enum { has_init_member_function  = DeduceInit<>::value };
  enum { has_final_member_function = DeduceFinal<>::value };


  template< class MemorySpace = typename execution_space::memory_space >
  struct Reducer
  {
  private:

    Functor const * const m_functor ;
    ValueType     * const m_result ;

    template< bool IsArray >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< IsArray , FunctorAnalysis::ValueType * >::type
    ref() const noexcept { return m_result ; }

    template< bool IsArray >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< ! IsArray , FunctorAnalysis::ValueType & >::type
    ref() const noexcept { return *m_result ; }

    template< bool IsArray >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< IsArray , int >::type
    len() const noexcept { return m_functor->value_count ; }

    template< bool IsArray >
    KOKKOS_INLINE_FUNCTION constexpr
    typename std::enable_if< ! IsArray , int >::type
    len() const noexcept { return candidate_is_void ? 0 : 1 ; }

  public:

    using reducer        = Reducer ;
    using value_type     = FunctorAnalysis::value_type ;
    using memory_space   = MemorySpace ;
    using reference_type = FunctorAnalysis::reference_type ;
    using functor_type   = Functor ; // Adapts a functor

    KOKKOS_INLINE_FUNCTION constexpr
    value_type * data() const noexcept { return m_result ; }

    KOKKOS_INLINE_FUNCTION constexpr
    reference_type reference() const noexcept
      { return Reducer::template ref< candidate_is_array >(); }

    KOKKOS_INLINE_FUNCTION constexpr
    int length() const noexcept
      { return Reducer::template len< candidate_is_array >(); }

    KOKKOS_INLINE_FUNCTION
    void copy( ValueType * const dst
             , ValueType const * const src ) const noexcept
      { for ( int i = 0 ; i < Reducer::template len< candidate_is_array >() ; ++i ) dst[i] = src[i] ; }

    KOKKOS_INLINE_FUNCTION
    void join( ValueType volatile * dst
             , ValueType volatile const * src ) const noexcept
      { DeduceJoin<>::join( m_functor , dst , src ); }

    KOKKOS_INLINE_FUNCTION 
    void init( ValueType * dst ) const noexcept
      { DeduceInit<>::init( m_functor , dst ); }

    KOKKOS_INLINE_FUNCTION
    void final( ValueType * dst ) const noexcept
      { DeduceFinal<>::final( m_functor , dst ); }

    Reducer( Reducer const & ) = default ;
    Reducer( Reducer && ) = default ;
    Reducer & operator = ( Reducer const & ) = delete ;
    Reducer & operator = ( Reducer && ) = delete ;

    template< class S >
    using rebind = Reducer< S > ;

    KOKKOS_INLINE_FUNCTION explicit constexpr
    Reducer( Functor const * arg_functor = 0
           , ValueType * arg_value = 0 ) noexcept
      : m_functor(arg_functor), m_result(arg_value) {}
  };
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_FUNCTORANALYSIS_HPP */

