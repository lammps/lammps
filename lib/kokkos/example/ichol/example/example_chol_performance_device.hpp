#pragma once
#ifndef __EXAMPLE_CHOL_PERFORMANCE_DEVICE_HPP__
#define __EXAMPLE_CHOL_PERFORMANCE_DEVICE_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "graph_helper_scotch.hpp"
#include "symbolic_factor_helper.hpp"
#include "crs_matrix_helper.hpp"

#include "task_view.hpp"

#include "task_factory.hpp"

#include "chol.hpp"

namespace Tacho {

  using namespace std;

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void>
  int exampleCholPerformanceDevice(const string file_input,
                                   const int treecut,
                                   const int prunecut,
                                   const int seed,
                                   const int nthreads,
                                   const int max_task_dependence,
                                   const int max_concurrency,
                                   const int team_size,
                                   const int fill_level,
                                   const int league_size,
                                   const bool skip_serial,
                                   const bool verbose) {
    typedef ValueType   value_type;
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;
    typedef typename
       Kokkos::Impl::is_space< SpaceType >::host_mirror_space::execution_space
         HostSpaceType ;

    typedef TaskFactory<Kokkos::Experimental::TaskPolicy<SpaceType>,
      Kokkos::Experimental::Future<int,SpaceType> > TaskFactoryType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType>
      CrsMatrixBaseType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>
      CrsMatrixBaseHostType;

    typedef Kokkos::MemoryUnmanaged MemoryUnmanaged ;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,SpaceType,MemoryUnmanaged >
      CrsMatrixNestedType;


    typedef GraphHelper_Scotch<CrsMatrixBaseHostType> GraphHelperType;
    typedef SymbolicFactorHelper<CrsMatrixBaseHostType> SymbolicFactorHelperType;

    typedef CrsMatrixView<CrsMatrixNestedType> CrsMatrixViewType;
    typedef TaskView<CrsMatrixViewType,TaskFactoryType> CrsTaskViewType;

    typedef CrsMatrixBase<CrsTaskViewType,ordinal_type,size_type,SpaceType> CrsHierMatrixBaseType;

    typedef CrsMatrixView<CrsHierMatrixBaseType> CrsHierMatrixViewType;
    typedef TaskView<CrsHierMatrixViewType,TaskFactoryType> CrsHierTaskViewType;

    int r_val = 0;

    Kokkos::Timer timer;
    double
      t_import = 0.0,
      t_reorder = 0.0,
      t_symbolic = 0.0,
      t_flat2hier = 0.0,
      t_factor_task = 0.0;

    cout << "CholPerformanceDevice:: import input file = " << file_input << endl;
    CrsMatrixBaseHostType AA("AA");
    {
      timer.reset();

      ifstream in;
      in.open(file_input);
      if (!in.good()) {
        cout << "Failed in open the file: " << file_input << endl;
        return ++r_val;
      }
      AA.importMatrixMarket(in);

      t_import = timer.seconds();

      if (verbose) {
        AA.showMe( std::cout );
        std::cout << endl;
      }
    }
    cout << "CholPerformanceDevice:: import input file::time = " << t_import << endl;

    cout << "CholPerformanceDevice:: reorder the matrix" << endl;
    CrsMatrixBaseHostType PA("Permuted AA");

    // '*_UU' is the permuted base upper triangular matrix
    CrsMatrixBaseHostType host_UU("host_UU");
    CrsMatrixBaseType     device_UU("UU");
    CrsHierMatrixBaseType device_HU("HU");;

    // typename CrsMatrixBaseHostType host_UU("host_UU");

    {
      typename GraphHelperType::size_type_array rptr("Graph::RowPtrArray", AA.NumRows() + 1);
      typename GraphHelperType::ordinal_type_array cidx("Graph::ColIndexArray", AA.NumNonZeros());

      AA.convertGraph(rptr, cidx);
      GraphHelperType S("ScotchHelper",
                        AA.NumRows(),
                        rptr,
                        cidx,
                        seed);
      {
        timer.reset();

        S.computeOrdering(treecut, 0);
        S.pruneTree(prunecut);

        PA.copy(S.PermVector(), S.InvPermVector(), AA);

        t_reorder = timer.seconds();

        if (verbose) {
          S.showMe( std::cout );
          std::cout << std::endl ;
          PA.showMe( std::cout );
          std::cout << std::endl ;
        }
      }

      // Symbolic factorization adds non-zero entries
      // for factorization levels.
      // Runs on the host process and currently requires std::sort.

      cout << "CholPerformanceDevice:: reorder the matrix::time = " << t_reorder << endl;
      {
        SymbolicFactorHelperType F(PA, league_size);
        timer.reset();
        F.createNonZeroPattern(fill_level, Uplo::Upper, host_UU);
        t_symbolic = timer.seconds();
        cout << "CholPerformanceDevice:: AA (nnz) = " << AA.NumNonZeros() << ", host_UU (nnz) = " << host_UU.NumNonZeros() << endl;

        if (verbose) {
          F.showMe( std::cout );
          std::cout << std::endl ;
          host_UU.showMe( std::cout );
          std::cout << std::endl ;
        }
      }
      cout << "CholPerformanceDevice:: symbolic factorization::time = " << t_symbolic << endl;

    //----------------------------------------------------------------------
    // Allocate device_UU conformal to host_UU 
    // and deep_copy host_UU arrays to device_UU arrays.
    // Set up device_HU referencing blocks of device_UU

      {
        timer.reset();

        device_UU.copy( host_UU );

        CrsMatrixHelper::flat2hier(Uplo::Upper, device_UU, device_HU,
                                   S.NumBlocks(),
                                   S.RangeVector(),
                                   S.TreeVector());

        // Filling non-zero block matrixes' row ranges within block view.
        // This is performed entirely in the 'device_HU' space.

        CrsMatrixHelper::fillRowViewArray( device_HU );

        t_flat2hier = timer.seconds();

        cout << "CholPerformanceDevice:: Hier (dof, nnz) = " << device_HU.NumRows() << ", " << device_HU.NumNonZeros() << endl;
      }
      cout << "CholPerformanceDevice:: copy base matrix and construct hierarchical matrix::time = " << t_flat2hier << endl;
    }

    cout << "CholPerformanceDevice:: max concurrency = " << max_concurrency << endl;

    const size_t max_task_size = 4*sizeof(CrsTaskViewType)+128;
    cout << "CholPerformanceDevice:: max task size   = " << max_task_size << endl;

    //----------------------------------------------------------------------
    // From here onward all work is on the device.
    //----------------------------------------------------------------------

    {
      typename TaskFactoryType::policy_type policy(max_concurrency,
                                                   max_task_size,
                                                   max_task_dependence,
                                                   team_size);

      cout << "CholPerformanceDevice:: ByBlocks factorize the matrix:: team_size = " << team_size << endl;
      CrsHierTaskViewType H( device_HU );
      {
        timer.reset();
        {
          // auto future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks>::
          auto future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::Two>::
                                                TaskFunctor<CrsHierTaskViewType>(policy,H), 0);
          policy.spawn(future);
          Kokkos::Experimental::wait(policy);
        }
        t_factor_task += timer.seconds();

        cout << "CholPerformanceDevice:: policy.allocated_task_count = "
             << policy.allocated_task_count()
             << endl ;

        if (verbose) {
          host_UU.copy( device_UU );
          host_UU.showMe( std::cout );
          std::cout << endl;
        }
      }
      cout << "CholPerformanceDevice:: ByBlocks factorize the matrix::time = " << t_factor_task << endl;
    }

    return r_val;
  }
}

#endif
