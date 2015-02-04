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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_QTHREAD_TASKPOLICY_HPP
#define KOKKOS_QTHREAD_TASKPOLICY_HPP

#include <string>
#include <typeinfo>
#include <stdexcept>

#include <qthread.h>

#include <Kokkos_Qthread.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <Kokkos_View.hpp>

#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class TaskMember< Kokkos::Qthread , void , void >
{
public:

  typedef void         (* function_apply_type)  ( TaskMember * );
  typedef void         (* function_dealloc_type)( TaskMember * );
  typedef TaskMember * (* function_verify_type) ( TaskMember * );

private:

  const function_dealloc_type  m_dealloc ; ///< Deallocation
  const function_verify_type   m_verify ;  ///< Result type verification
  const function_apply_type    m_apply ;   ///< Apply function
  TaskMember ** const          m_dep ;     ///< Dependences
  const int                    m_dep_capacity ; ///< Capacity of dependences
  int                          m_dep_size ;     ///< Actual count of dependences
  int                          m_ref_count ;    ///< Reference count
  int                          m_state ;        ///< State of the task
  aligned_t                    m_qfeb ;         ///< Qthread full/empty bit

  TaskMember() /* = delete */ ;
  TaskMember( const TaskMember & ) /* = delete */ ;
  TaskMember & operator = ( const TaskMember & ) /* = delete */ ;

  static aligned_t qthread_func( void * arg );

  static void * allocate( const unsigned arg_sizeof_derived , const unsigned arg_dependence_capacity );
  static void   deallocate( void * );

  void throw_error_add_dependence() const ;
  static void throw_error_verify_type();

  template < class DerivedTaskType >
  static
  void deallocate( TaskMember * t )
    {
      DerivedTaskType * ptr = static_cast< DerivedTaskType * >(t);
      ptr->~DerivedTaskType();
      deallocate( (void *) ptr );
    }

protected :

  ~TaskMember();

  // Used by TaskMember< Qthread , ResultType , void >
  TaskMember( const function_verify_type   arg_verify
            , const function_dealloc_type  arg_dealloc
            , const function_apply_type    arg_apply
            , const unsigned               arg_sizeof_derived
            , const unsigned               arg_dependence_capacity
            );

  // Used for TaskMember< Qthread , void , void >
  TaskMember( const function_dealloc_type  arg_dealloc
            , const function_apply_type    arg_apply
            , const unsigned               arg_sizeof_derived
            , const unsigned               arg_dependence_capacity
            );

public:

  template< typename ResultType >
  KOKKOS_FUNCTION static
  TaskMember * verify_type( TaskMember * t )
    {
      enum { check_type = ! Impl::is_same< ResultType , void >::value };

      if ( check_type && t != 0 ) {

        // Verify that t->m_verify is this function
        const function_verify_type self = & TaskMember::template verify_type< ResultType > ;

        if ( t->m_verify != self ) {
          t = 0 ;
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
          throw_error_verify_type();
#endif
        }
      }
      return t ;
    }

  //----------------------------------------
  /*  Inheritence Requirements on task types:
   *    typedef  FunctorType::value_type  value_type ;
   *    class DerivedTaskType
   *      : public TaskMember< Qthread , value_type , FunctorType >
   *      { ... };
   *    class TaskMember< Qthread , value_type , FunctorType >
   *      : public TaskMember< Qthread , value_type , void >
   *      , public Functor
   *      { ... };
   *  If value_type != void
   *    class TaskMember< Qthread , value_type , void >
   *      : public TaskMember< Qthread , void , void >
   *
   *  Allocate space for DerivedTaskType followed by TaskMember*[ dependence_capacity ]
   *
   */

  /** \brief  Allocate and construct a single-thread task */
  template< class DerivedTaskType >
  static
  TaskMember * create( const typename DerivedTaskType::functor_type &  arg_functor
                     , const unsigned                                  arg_dependence_capacity )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;
      typedef typename functor_type::value_type       value_type ;

      DerivedTaskType * const task =
        new( allocate( sizeof(DerivedTaskType) , arg_dependence_capacity ) )
          DerivedTaskType( & TaskMember::template deallocate< DerivedTaskType >
                         , & TaskMember::template apply_single< functor_type , value_type >
                         , sizeof(DerivedTaskType)
                         , arg_dependence_capacity
                         , arg_functor );

      return static_cast< TaskMember * >( task );
    }

  /** \brief  Allocate and construct a data parallel task */
  template< class DerivedTaskType >
  static
  TaskMember * create( const typename DerivedTaskType::policy_type &   arg_policy
                     , const typename DerivedTaskType::functor_type &  arg_functor
                     , const unsigned                                  arg_dependence_capacity )
    {
      DerivedTaskType * const task =
        new( allocate( sizeof(DerivedTaskType) , arg_dependence_capacity ) )
          DerivedTaskType( & TaskMember::template deallocate< DerivedTaskType >
                         , sizeof(DerivedTaskType)
                         , arg_dependence_capacity
                         , arg_policy
                         , arg_functor
                         );

      return static_cast< TaskMember * >( task );
    }

  void schedule();
  static void wait( const Future< void , Kokkos::Qthread > & );

  //----------------------------------------

  typedef FutureValueTypeIsVoidError get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return get_result_type() ; }

  KOKKOS_INLINE_FUNCTION
  Kokkos::TaskState get_state() const { return Kokkos::TaskState( m_state ); }

  //----------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  static
  void assign( TaskMember ** const lhs , TaskMember * const rhs , const bool no_throw = false );
#else
  KOKKOS_INLINE_FUNCTION static
  void assign( TaskMember ** const lhs , TaskMember * const rhs , const bool no_throw = false ) {}
#endif

  KOKKOS_INLINE_FUNCTION
  TaskMember * get_dependence( int i ) const
    { return ( Kokkos::TASK_STATE_EXECUTING == m_state && 0 <= i && i < m_dep_size ) ? m_dep[i] : (TaskMember*) 0 ; }

  KOKKOS_INLINE_FUNCTION
  int get_dependence() const
    { return m_dep_size ; }

  KOKKOS_INLINE_FUNCTION
  void clear_dependence()
    {
      for ( int i = 0 ; i < m_dep_size ; ++i ) assign( m_dep + i , 0 );
      m_dep_size = 0 ;
    }

  KOKKOS_INLINE_FUNCTION
  void add_dependence( TaskMember * before )
    {
      if ( ( Kokkos::TASK_STATE_CONSTRUCTING == m_state ||
             Kokkos::TASK_STATE_EXECUTING    == m_state ) &&
           m_dep_size < m_dep_capacity ) {
        assign( m_dep + m_dep_size , before );
        ++m_dep_size ;
      }
      else {
        throw_error_add_dependence();
      }
    }

  //----------------------------------------

  template< class FunctorType , class ResultType >
  KOKKOS_INLINE_FUNCTION static
  void apply_single( typename Impl::enable_if< ! Impl::is_same< ResultType , void >::value , TaskMember * >::type t )
    {
      typedef TaskMember< Kokkos::Qthread , ResultType , FunctorType > derived_type ;

      // TaskMember< Kokkos::Qthread , ResultType , FunctorType >
      //   : public TaskMember< Kokkos::Qthread , ResultType , void >
      //   , public FunctorType
      //   { ... };

      derived_type & m = * static_cast< derived_type * >( t );

      Impl::FunctorApply< FunctorType , void , ResultType & >::apply( (FunctorType &) m , & m.m_result );
    }

  template< class FunctorType , class ResultType >
  KOKKOS_INLINE_FUNCTION static
  void apply_single( typename Impl::enable_if< Impl::is_same< ResultType , void >::value , TaskMember * >::type t )
    {
      typedef TaskMember< Kokkos::Qthread , ResultType , FunctorType > derived_type ;

      // TaskMember< Kokkos::Qthread , ResultType , FunctorType >
      //   : public TaskMember< Kokkos::Qthread , ResultType , void >
      //   , public FunctorType
      //   { ... };

      derived_type & m = * static_cast< derived_type * >( t );

      Impl::FunctorApply< FunctorType , void , void >::apply( (FunctorType &) m );
    }
};

//----------------------------------------------------------------------------
/** \brief  Base class for tasks with a result value in the Qthread execution space.
 *
 *  The FunctorType must be void because this class is accessed by the
 *  Future class for the task and result value.
 *
 *  Must be derived from TaskMember<S,void,void> 'root class' so the Future class
 *  can correctly static_cast from the 'root class' to this class.
 */
template < class ResultType >
class TaskMember< Kokkos::Qthread , ResultType , void >
  : public TaskMember< Kokkos::Qthread , void , void >
{
public:

  ResultType  m_result ;

  typedef const ResultType & get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return m_result ; }

protected:

  typedef TaskMember< Kokkos::Qthread , void , void >  task_root_type ;
  typedef task_root_type::function_dealloc_type       function_dealloc_type ;
  typedef task_root_type::function_apply_type         function_apply_type ;

  inline
  TaskMember( const function_dealloc_type  arg_dealloc
            , const function_apply_type    arg_apply
            , const unsigned               arg_sizeof_derived
            , const unsigned               arg_dependence_capacity
            )
    : task_root_type( & task_root_type::template verify_type< ResultType >
                    , arg_dealloc
                    , arg_apply
                    , arg_sizeof_derived
                    , arg_dependence_capacity )
    , m_result()
    {}

};

template< class ResultType , class FunctorType >
class TaskMember< Kokkos::Qthread , ResultType , FunctorType >
  : public TaskMember< Kokkos::Qthread , ResultType , void >
  , public FunctorType
{
public:

  typedef FunctorType  functor_type ;

  typedef TaskMember< Kokkos::Qthread , void , void >        task_root_type ;
  typedef TaskMember< Kokkos::Qthread , ResultType , void >  task_base_type ;
  typedef task_root_type::function_dealloc_type             function_dealloc_type ;
  typedef task_root_type::function_apply_type               function_apply_type ;

  inline
  TaskMember( const function_dealloc_type  arg_dealloc
            , const function_apply_type    arg_apply
            , const unsigned               arg_sizeof_derived
            , const unsigned               arg_dependence_capacity
            , const functor_type &         arg_functor
            )
    : task_base_type( arg_dealloc , arg_apply , arg_sizeof_derived , arg_dependence_capacity )
    , functor_type( arg_functor )
    {}
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  ForEach task in the Qthread execution space
 *
 *  Derived from TaskMember< Kokkos::Qthread , ResultType , FunctorType >
 *  so that Functor can be cast to task root type without knowing policy.
 */
template< class Arg0 , class Arg1 , class Arg2 , class ResultType , class FunctorType >
class TaskForEach< Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >
                 , ResultType
                 , FunctorType >
  : TaskMember< Kokkos::Qthread , ResultType , FunctorType >
{
public:

  typedef FunctorType                                              functor_type ;
  typedef RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >       policy_type ;

private:

  friend class Kokkos::TaskPolicy< Kokkos::Qthread > ;
  friend class Kokkos::Impl::TaskMember< Kokkos::Qthread , void , void > ;

  typedef TaskMember< Kokkos::Qthread , void , void >               task_root_type ;
  typedef TaskMember< Kokkos::Qthread , ResultType , FunctorType >  task_base_type ;
  typedef task_root_type::function_dealloc_type                    function_dealloc_type ;

  policy_type  m_policy ;

  template< class Tag >
  inline
  typename Impl::enable_if< Impl::is_same<Tag,void>::value >::type
    apply_policy() const
    {
      const typename policy_type::member_type e = m_policy.end();
      for ( typename policy_type::member_type i = m_policy.begin() ; i < e ; ++i ) {
        functor_type::operator()(i);
      }
    }

  template< class Tag >
  inline
  typename Impl::enable_if< ! Impl::is_same<Tag,void>::value >::type
    apply_policy() const
    {
      const Tag tag ;
      const typename policy_type::member_type e = m_policy.end();
      for ( typename policy_type::member_type i = m_policy.begin() ; i < e ; ++i ) {
        functor_type::operator()(tag,i);
      }
    }

  static
  void apply_parallel( task_root_type * t )
    {
      static_cast<TaskForEach*>(t)->template apply_policy< typename policy_type::work_tag >();

      task_root_type::template apply_single< functor_type , ResultType >( t );
    }

  TaskForEach( const function_dealloc_type  arg_dealloc
             , const int                    arg_sizeof_derived
             , const int                    arg_dependence_capacity
             , const policy_type &          arg_policy
             , const functor_type &         arg_functor
             )
    : task_base_type( arg_dealloc
                    , & apply_parallel
                    , arg_sizeof_derived
                    , arg_dependence_capacity
                    , arg_functor )
    , m_policy( arg_policy )
    {}

  TaskForEach() /* = delete */ ;
  TaskForEach( const TaskForEach & ) /* = delete */ ;
  TaskForEach & operator = ( const TaskForEach & ) /* = delete */ ;
};

//----------------------------------------------------------------------------
/** \brief  Reduce task in the Qthread execution space
 *
 *  Derived from TaskMember< Kokkos::Qthread , ResultType , FunctorType >
 *  so that Functor can be cast to task root type without knowing policy.
 */
template< class Arg0 , class Arg1 , class Arg2 , class ResultType , class FunctorType >
class TaskReduce< Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >
                , ResultType
                , FunctorType >
  : TaskMember< Kokkos::Qthread , ResultType , FunctorType >
{
public:

  typedef FunctorType                                              functor_type ;
  typedef RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Qthread >       policy_type ;

private:

  friend class Kokkos::TaskPolicy< Kokkos::Qthread > ;
  friend class Kokkos::Impl::TaskMember< Kokkos::Qthread , void , void > ;

  typedef TaskMember< Kokkos::Qthread , void , void >               task_root_type ;
  typedef TaskMember< Kokkos::Qthread , ResultType , FunctorType >  task_base_type ;
  typedef task_root_type::function_dealloc_type                    function_dealloc_type ;

  policy_type  m_policy ;

  template< class Tag >
  inline
  void apply_policy( typename Impl::enable_if< Impl::is_same<Tag,void>::value , ResultType & >::type result ) const
    {
      Impl::FunctorValueInit< functor_type , Tag >::init( *this , & result );
      const typename policy_type::member_type e = m_policy.end();
      for ( typename policy_type::member_type i = m_policy.begin() ; i < e ; ++i ) {
        functor_type::operator()( i, result );
      }
    }

  template< class Tag >
  inline
  void apply_policy( typename Impl::enable_if< ! Impl::is_same<Tag,void>::value , ResultType & >::type result ) const
    {
      Impl::FunctorValueInit< functor_type , Tag >::init( *this , & result );
      const Tag tag ;
      const typename policy_type::member_type e = m_policy.end();
      for ( typename policy_type::member_type i = m_policy.begin() ; i < e ; ++i ) {
        functor_type::operator()( tag, i, result );
      }
    }

  static
  void apply_parallel( task_root_type * t )
    {
      TaskReduce * const task = static_cast<TaskReduce*>(t);

      task->template apply_policy< typename policy_type::work_tag >( task->task_base_type::m_result );

      task_root_type::template apply_single< functor_type , ResultType >( t );
    }

  TaskReduce( const function_dealloc_type  arg_dealloc
            , const int                    arg_sizeof_derived
            , const int                    arg_dependence_capacity
            , const policy_type &          arg_policy
            , const functor_type &         arg_functor
            )
    : task_base_type( arg_dealloc
                    , & apply_parallel
                    , arg_sizeof_derived
                    , arg_dependence_capacity
                    , arg_functor )
    , m_policy( arg_policy )
    {}

  TaskReduce() /* = delete */ ;
  TaskReduce( const TaskReduce & ) /* = delete */ ;
  TaskReduce & operator = ( const TaskReduce & ) /* = delete */ ;
};


} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template<>
class TaskPolicy< Kokkos::Qthread >
{
public:

  typedef Kokkos::Qthread execution_space ;

private:

  typedef Impl::TaskMember< execution_space , void , void > task_root_type ;

  TaskPolicy & operator = ( const TaskPolicy & ) /* = delete */ ;

  template< class FunctorType >
  static inline
  const task_root_type * get_task_root( const FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< const task_root_type * >( static_cast< const task_type * >(f) );
    }

  template< class FunctorType >
  static inline
  task_root_type * get_task_root( FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< task_root_type * >( static_cast< task_type * >(f) );
    }

  const unsigned m_default_dependence_capacity ;

public:

  KOKKOS_INLINE_FUNCTION
  TaskPolicy() : m_default_dependence_capacity(4) {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicy( const TaskPolicy & rhs ) : m_default_dependence_capacity( rhs.m_default_dependence_capacity ) {}

  KOKKOS_INLINE_FUNCTION
  explicit
  TaskPolicy( const unsigned arg_default_dependence_capacity )
    : m_default_dependence_capacity( arg_default_dependence_capacity ) {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicy( const TaskPolicy &
            , const unsigned arg_default_dependence_capacity )
    : m_default_dependence_capacity( arg_default_dependence_capacity ) {}

  //----------------------------------------

  template< class ValueType >
  const Future< ValueType , execution_space > &
    spawn( const Future< ValueType , execution_space > & f ) const
      {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        f.m_task->schedule();
#endif
        return f ;
      }

  // Create single-thread task

  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create( const FunctorType & functor
        , const unsigned dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type value_type ;
      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;
      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type >(
          functor , ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity ) )
#endif
        );
    }

  // Create parallel foreach task

  template< class PolicyType , class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create_foreach( const PolicyType  & policy
                , const FunctorType & functor
                , const unsigned      dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type value_type ;
      typedef Impl::TaskForEach< PolicyType , value_type , FunctorType > task_type ;
      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type >( policy , functor ,
          ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity ) )
#endif
        );
    }

  // Create parallel reduce task

  template< class PolicyType , class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create_reduce( const PolicyType  & policy
               , const FunctorType & functor
               , const unsigned      dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type value_type ;
      typedef Impl::TaskReduce< PolicyType , value_type , FunctorType > task_type ;
      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type >( policy , functor ,
          ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity ) )
#endif
        );
    }

  // Add dependence
  template< class A1 , class A2 , class A3 , class A4 >
  void add_dependence( const Future<A1,A2> & after
                     , const Future<A3,A4> & before
                     , typename Impl::enable_if
                        < Impl::is_same< typename Future<A1,A2>::execution_space , execution_space >::value
                          &&
                          Impl::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      )
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      after.m_task->add_dependence( before.m_task );
#endif
    }

  //----------------------------------------
  // Functions for an executing task functor to query dependences,
  // set new dependences, and respawn itself.

  template< class FunctorType >
  Future< void , execution_space >
  get_dependence( const FunctorType * task_functor , int i ) const
    {
      return Future<void,execution_space>(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        get_task_root(task_functor)->get_dependence(i)
#endif
        );
    }

  template< class FunctorType >
  int get_dependence( const FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return get_task_root(task_functor)->get_dependence(); }
#else
    { return 0 ; }
#endif

  template< class FunctorType >
  void clear_dependence( FunctorType * task_functor ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      get_task_root(task_functor)->clear_dependence();
#endif
    }

  template< class FunctorType , class A3 , class A4 >
  void add_dependence( FunctorType * task_functor
                     , const Future<A3,A4> & before
                     , typename Impl::enable_if
                        < Impl::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      )
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      get_task_root(task_functor)->add_dependence( before.m_task );
#endif
    }

  template< class FunctorType >
  void respawn( FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { get_task_root(task_functor)->schedule(); }
#else
    {}
#endif

};

inline
void wait( TaskPolicy< Kokkos::Qthread > & );

inline
void wait( const Future< void , Kokkos::Qthread > & future )
{ Impl::TaskMember< Kokkos::Qthread , void , void >::wait( future ); }

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_QTHREAD_TASK_HPP */

