#pragma once
#ifndef __CHOL_U_RIGHT_LOOK_BY_BLOCKS_HPP__
#define __CHOL_U_RIGHT_LOOK_BY_BLOCKS_HPP__

/// \file chol_u_right_look_by_blocks.hpp
/// \brief Cholesky factorization by-blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

/// The Partitioned-Block Matrix (PBM) is sparse and a block itself is a view of a sparse matrix. 
/// The algorithm generates tasks with a given sparse block matrix structure.

// basic utils
#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho { 
  
  using namespace std;

  template< typename CrsTaskViewType >
  KOKKOS_INLINE_FUNCTION
  int releaseFutures( typename CrsTaskViewType::matrix_type & A )
    {
      typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
      typedef typename CrsTaskViewType::row_view_type     row_view_type;
      typedef typename CrsTaskViewType::future_type       future_type;
      
      row_view_type a(A,0);
      
      const ordinal_type nnz = a.NumNonZeros();

      for (ordinal_type j=0;j<nnz;++j) {
        a.Value(j).setFuture( future_type() );
      }

      return nnz ;
    }
  
  // ========================================
  // detailed workflow of by-blocks algorithm
  // ========================================
  template<int ArgVariant, 
           template<int,int> class ControlType,
           typename CrsTaskViewType>
  class CholUpperRightLookByBlocks {
  public:
    KOKKOS_INLINE_FUNCTION
    static int genScalarTask(typename CrsTaskViewType::policy_type &policy,
                             typename CrsTaskViewType::matrix_type &A) {
      typedef typename CrsTaskViewType::value_type        value_type;
      typedef typename CrsTaskViewType::row_view_type     row_view_type;
      
      typedef typename CrsTaskViewType::future_type       future_type;
      typedef typename CrsTaskViewType::task_factory_type task_factory_type;
      
      row_view_type a(A, 0); 
      value_type &aa = a.Value(0);
      
      // construct a task
      future_type f = task_factory_type::create(policy,
                                                typename Chol<Uplo::Upper,
                                                CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Chol)>
                                                ::template TaskFunctor<value_type>(policy,aa));
      

if ( false ) {
 printf("Chol [%d +%d)x[%d +%d) spawn depend %d\n"
       , aa.OffsetRows()
       , aa.NumRows()
       , aa.OffsetCols()
       , aa.NumCols()
       , int( ! aa.Future().is_null() )
       );
}

      // manage dependence
      task_factory_type::addDependence(policy, f, aa.Future());
      aa.setFuture(f);

      // spawn a task
      task_factory_type::spawn(policy, f, true /* high priority */ );
      
      return 1;
    }
    
    KOKKOS_INLINE_FUNCTION
    static int genTrsmTasks(typename CrsTaskViewType::policy_type &policy,
                            typename CrsTaskViewType::matrix_type &A,
                            typename CrsTaskViewType::matrix_type &B) {
      typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
      typedef typename CrsTaskViewType::row_view_type     row_view_type;
      typedef typename CrsTaskViewType::value_type        value_type;

      typedef typename CrsTaskViewType::future_type       future_type;
      typedef typename CrsTaskViewType::task_factory_type task_factory_type;
      
      row_view_type a(A,0), b(B,0); 
      value_type &aa = a.Value(0);

if ( false ) {
  printf("genTrsmTasks after aa.Future().reference_count = %d\n"
        , aa.Future().reference_count());
}
      const ordinal_type nnz = b.NumNonZeros();
      for (ordinal_type j=0;j<nnz;++j) {
        typedef typename
           Trsm< Side::Left,Uplo::Upper,Trans::ConjTranspose,
                 CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Trsm)>
           ::template TaskFunctor<double,value_type,value_type>
             FunctorType ;

        value_type &bb = b.Value(j);
        
        future_type f = task_factory_type
          ::create(policy, FunctorType(policy,Diag::NonUnit, 1.0, aa, bb));
        
if ( false ) {
 printf("Trsm [%d +%d)x[%d +%d) spawn depend %d %d\n"
       , bb.OffsetRows()
       , bb.NumRows()
       , bb.OffsetCols()
       , bb.NumCols()
       , int( ! aa.Future().is_null() )
       , int( ! bb.Future().is_null() )
       );
}

        // trsm dependence
        task_factory_type::addDependence(policy, f, aa.Future());
        
        // self
        task_factory_type::addDependence(policy, f, bb.Future());
        
        // place task signature on b
        bb.setFuture(f);
        
        // spawn a task
        task_factory_type::spawn(policy, f, true /* high priority */);              
      }

if ( false ) {
  printf("genTrsmTasks after aa.Future().reference_count = %d\n"
        , aa.Future().reference_count());
}
      
      return nnz ;
    }
    
    KOKKOS_INLINE_FUNCTION
    static int genHerkTasks(typename CrsTaskViewType::policy_type &policy,
                            typename CrsTaskViewType::matrix_type &A,
                            typename CrsTaskViewType::matrix_type &C) {
      typedef typename CrsTaskViewType::ordinal_type      ordinal_type;
      typedef typename CrsTaskViewType::value_type        value_type;
      typedef typename CrsTaskViewType::row_view_type     row_view_type;
      
      typedef typename CrsTaskViewType::future_type       future_type;
      typedef typename CrsTaskViewType::task_factory_type task_factory_type;
      
      // case that X.transpose, A.no_transpose, Y.no_transpose
      
      row_view_type a(A,0), c; 
      
      const ordinal_type nnz = a.NumNonZeros();
      ordinal_type herk_count = 0 ; 
      ordinal_type gemm_count = 0 ; 

      // update herk
      for (ordinal_type i=0;i<nnz;++i) {
        const ordinal_type row_at_i = a.Col(i);
        value_type &aa = a.Value(i);
        
        c.setView(C, row_at_i);
        
        ordinal_type idx = 0;
        for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
          const ordinal_type col_at_j = a.Col(j);
          value_type &bb = a.Value(j);
          
          if (row_at_i == col_at_j) {
            idx = c.Index(row_at_i, idx);
            if (idx >= 0) {
              ++herk_count ;
              value_type &cc = c.Value(idx);
              future_type f = task_factory_type
                ::create(policy, 
                         typename Herk<Uplo::Upper,Trans::ConjTranspose,
                         CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Herk)>
                         ::template TaskFunctor<double,value_type,value_type>(policy,-1.0, aa, 1.0, cc));
            

if ( false ) {
 printf("Herk [%d +%d)x[%d +%d) spawn %d %d\n"
       , cc.OffsetRows()
       , cc.NumRows()
       , cc.OffsetCols()
       , cc.NumCols()
       , int( ! aa.Future().is_null() )
       , int( ! cc.Future().is_null() )
       );
}

              // dependence
              task_factory_type::addDependence(policy, f, aa.Future());              
            
              // self
              task_factory_type::addDependence(policy, f, cc.Future());
            
              // place task signature on y
              cc.setFuture(f);

              // spawn a task
              task_factory_type::spawn(policy, f);
            }
          } else {
            idx = c.Index(col_at_j, idx);
            if (idx >= 0) {
              ++gemm_count ;
              value_type &cc = c.Value(idx);
              future_type f = task_factory_type
                ::create(policy, 
                         typename Gemm<Trans::ConjTranspose,Trans::NoTranspose,
                         CtrlDetail(ControlType,AlgoChol::ByBlocks,ArgVariant,Gemm)>
                         ::template TaskFunctor<double,value_type,value_type,value_type>(policy,-1.0, aa, bb, 1.0, cc));
            

if ( false ) {
 printf("Gemm [%d +%d)x[%d +%d) spawn %d %d %d\n"
       , cc.OffsetRows()
       , cc.NumRows()
       , cc.OffsetCols()
       , cc.NumCols()
       , int( ! aa.Future().is_null() )
       , int( ! bb.Future().is_null() )
       , int( ! cc.Future().is_null() )
       );
}
 
              // dependence
              task_factory_type::addDependence(policy, f, aa.Future());
              task_factory_type::addDependence(policy, f, bb.Future());
            
              // self
              task_factory_type::addDependence(policy, f, cc.Future());
            
              // place task signature on y
              cc.setFuture(f);
            
              // spawn a task
              task_factory_type::spawn(policy, f);
            }
          }
        }
      }

if ( false ) {
printf("genHerkTask Herk(%ld) Gemm(%ld)\n",(long)herk_count,(long)gemm_count);
}
    
      return herk_count + gemm_count ;
    }
    
  };
  
  // specialization for different task generation in right looking by-blocks algorithm
  // =================================================================================
  template<int ArgVariant, template<int,int> class ControlType>
  class Chol<Uplo::Upper,AlgoChol::RightLookByBlocks,ArgVariant,ControlType> {
  public:

    // function interface
    // ==================
    template<typename ExecViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename ExecViewType::policy_type &policy, 
                      const typename ExecViewType::policy_type::member_type &member, 
                      typename ExecViewType::matrix_type & A,
                      int checkpoint )
      {
        typedef typename ExecViewType::row_view_type  row_view_type ;

        enum { CYCLE = 2 };

        typename ExecViewType::matrix_type
          ATL, ATR,      A00, A01, A02,
          ABL, ABR,      A10, A11, A12,
                         A20, A21, A22;

        Part_2x2(A,  ATL, ATR,
                 /**/ABL, ABR,
                 checkpoint, checkpoint, Partition::TopLeft);

        int tasks_spawned = 0 ;
        int futures_released = 0 ;

        for ( int i = 0 ; i < CYCLE && ATL.NumRows() < A.NumRows() ; ++i ) {
          Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                          /*******/ /**/  A10, A11, A12,
                          ABL, ABR, /**/  A20, A21, A22,
                          1, 1, Partition::BottomRight);
          // -----------------------------------------------------
          // Spawning tasks:

          // A11 = chol(A11) : #task = 1
          tasks_spawned +=
          CholUpperRightLookByBlocks<ArgVariant,ControlType,ExecViewType>
            ::genScalarTask(policy, A11);
          
          // A12 = inv(triu(A11)') * A12 : #tasks = non-zero row blocks
          tasks_spawned +=
          CholUpperRightLookByBlocks<ArgVariant,ControlType,ExecViewType>
            ::genTrsmTasks(policy, A11, A12);

          // A22 = A22 - A12' * A12 : #tasks = highly variable
          tasks_spawned +=
          CholUpperRightLookByBlocks<ArgVariant,ControlType,ExecViewType>
            ::genHerkTasks(policy, A12, A22);

          // -----------------------------------------------------
          // Can release futures of A11 and A12 

          futures_released += releaseFutures<ExecViewType>( A11 );
          futures_released += releaseFutures<ExecViewType>( A12 );

if ( false ) {
  printf("Chol iteration(%d) task_count(%d) cumulative: spawn(%d) release(%d)\n"
        , int(ATL.NumRows())
        , policy.allocated_task_count()
        , tasks_spawned , futures_released
        );
}

          // -----------------------------------------------------
          Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                           A10, A11, A12, /**/ /******/
                           A20, A21, A22, /**/ ABL, ABR,
                           Partition::TopLeft);

        }
      
      return ATL.NumRows();
    }
    
    // task-data parallel interface
    // ============================
    template<typename ExecViewType>
    class TaskFunctor {
    public:
      typedef typename ExecViewType::policy_type  policy_type;
      typedef typename ExecViewType::future_type  future_type;
      typedef typename policy_type::member_type   member_type;
      typedef int value_type;
      
    private:
      typename ExecViewType::matrix_type _A;
      
      policy_type _policy;
      int         _checkpoint ;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const policy_type & P ,
                  const typename ExecViewType::matrix_type & A)
        : _A(A),
          _policy(P),
          _checkpoint(0)
      { } 
      
      string Label() const { return "Chol"; }
      
      // task-data execution
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val)
      {
        if (member.team_rank() == 0) {
          // Clear out previous dependence
          _policy.clear_dependence( this );

          _checkpoint = Chol::invoke<ExecViewType>(_policy, member, _A, _checkpoint);

          if ( _checkpoint < _A.NumRows() ) _policy.respawn_needing_memory(this);

          r_val = 0 ;
        }
        return ;
      }

    };

  };
}

#endif
