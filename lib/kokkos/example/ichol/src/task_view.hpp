#pragma once
#ifndef __TASK_VIEW_HPP__
#define __TASK_VIEW_HPP__

/// \file task_view.hpp
/// \brief Task view is inherited from matrix view and have a member for the task handler.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  using namespace std;

  template<typename MatrixViewType,
           typename TaskFactoryType>
  class TaskView : public MatrixViewType {
  public:
    typedef          MatrixViewType                matrix_type ;
    typedef typename MatrixViewType::value_type    value_type;
    typedef typename MatrixViewType::ordinal_type  ordinal_type;

    typedef TaskFactoryType task_factory_type;
    typedef typename task_factory_type::policy_type policy_type;
    typedef typename task_factory_type::future_type future_type;

  private:
    future_type _f;

  public:
    KOKKOS_INLINE_FUNCTION
    void setFuture(const future_type &f)
      { _f = f; }

    KOKKOS_INLINE_FUNCTION
    future_type Future() const { return _f; }

    KOKKOS_INLINE_FUNCTION
    ~TaskView() = default ;

    KOKKOS_INLINE_FUNCTION
    TaskView() 
      : MatrixViewType(), _f()
    { } 

    TaskView(const TaskView &b) = delete ;

    KOKKOS_INLINE_FUNCTION
    TaskView(typename MatrixViewType::mat_base_type const & b) 
      : MatrixViewType(b), _f() 
    { }

    KOKKOS_INLINE_FUNCTION
    TaskView(typename MatrixViewType::mat_base_type const & b,
             const ordinal_type offm, const ordinal_type m,
             const ordinal_type offn, const ordinal_type n) 
      : MatrixViewType(b, offm, m, offn, n), _f() 
    { }

  };
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if ! KOKKOS_USING_EXP_VIEW

namespace Kokkos {
  namespace Impl {

    //  The Kokkos::View allocation will by default assign each allocated datum to zero.
    //  This is not the required initialization behavior when
    //  non-trivial objects are used within a Kokkos::View.
    //  Create a partial specialization of the Kokkos::Impl::AViewDefaultConstruct
    //  to replace the assignment initialization with placement new initialization.
    //
    //  This work-around is necessary until a TBD design refactorization of Kokkos::View.

    template< class ExecSpace , typename T1, typename T2 >
    struct ViewDefaultConstruct< ExecSpace , Tacho::TaskView<T1,T2> , true >
    {
      typedef Tacho::TaskView<T1,T2> type ;
      type * const m_ptr ;

      KOKKOS_FORCEINLINE_FUNCTION
      void operator()( const typename ExecSpace::size_type& i ) const
      { new(m_ptr+i) type(); }

      ViewDefaultConstruct( type * pointer , size_t capacity )
        : m_ptr( pointer )
      {
        Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
        parallel_for( range , *this );
        ExecSpace::fence();
      }
    };

  } // namespace Impl
} // namespace Kokkos

#endif /* #if ! KOKKOS_USING_EXP_VIEW */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
