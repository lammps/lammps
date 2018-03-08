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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP
#define KOKKOS_EXAMPLE_FENL_IMPL_HPP

#include <cmath>

// Kokkos libraries' headers:

#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <impl/Kokkos_Timer.hpp>

// Examples headers:

#include <BoxElemFixture.hpp>
#include <VectorImport.hpp>
#include <CGSolve.hpp>

#include <fenl.hpp>
#include <fenl_functors.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

inline
double maximum( MPI_Comm comm , double local )
{
  double global = local ;
#if defined( KOKKOS_ENABLE_MPI )
  MPI_Allreduce( & local , & global , 1 , MPI_DOUBLE , MPI_MAX , comm );
#endif
  return global ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

class ManufacturedSolution {
public:

  // Manufactured solution for one dimensional nonlinear PDE
  //
  //  -K T_zz + T^2 = 0 ; T(zmin) = T_zmin ; T(zmax) = T_zmax
  //
  //  Has an analytic solution of the form:
  //
  //    T(z) = ( a ( z - zmin ) + b )^(-2) where K = 1 / ( 6 a^2 )
  //
  //  Given T_0 and T_L compute K for this analytic solution.
  //
  //  Two analytic solutions:
  //
  //    Solution with singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) + 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( -1.0 / sqrt(T_zmin) )
  //
  //    Solution without singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) - 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( 1.0 / sqrt(T_zmin) )

  const double zmin ;
  const double zmax ;
  const double T_zmin ;
  const double T_zmax ;
  const double a ;
  const double b ;
  const double K ;

  ManufacturedSolution( const double arg_zmin ,
                        const double arg_zmax ,
                        const double arg_T_zmin ,
                        const double arg_T_zmax )
    : zmin( arg_zmin )
    , zmax( arg_zmax )
    , T_zmin( arg_T_zmin )
    , T_zmax( arg_T_zmax )
    , a( ( 1.0 / std::sqrt(T_zmax) - 1.0 / std::sqrt(T_zmin) ) / ( zmax - zmin ) )
    , b( 1.0 / std::sqrt(T_zmin) )
    , K( 1.0 / ( 6.0 * a * a ) )
    {}

  double operator()( const double z ) const
  {
    const double tmp = a * ( z - zmin ) + b ;
    return 1.0 / ( tmp * tmp );
  }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template < class Space , BoxElemPart::ElemOrder ElemOrder >
Perf fenl(
  MPI_Comm comm ,
  const int use_print ,
  const int use_trials ,
  const int use_atomic ,
  const int use_elems[] )
{
  typedef Kokkos::Example::BoxElemFixture< Space , ElemOrder > FixtureType ;

  typedef Kokkos::Example::CrsMatrix< double , Space >
    SparseMatrixType ;

  typedef typename SparseMatrixType::StaticCrsGraphType
    SparseGraphType ;

  typedef Kokkos::Example::FENL::NodeNodeGraph< typename FixtureType::elem_node_type , SparseGraphType , FixtureType::ElemNode >
     NodeNodeGraphType ;

  typedef Kokkos::Example::FENL::ElementComputation< FixtureType , SparseMatrixType >
    ElementComputationType ;

  typedef Kokkos::Example::FENL::DirichletComputation< FixtureType , SparseMatrixType >
    DirichletComputationType ;

  typedef NodeElemGatherFill< ElementComputationType >
    NodeElemGatherFillType ;

  typedef typename ElementComputationType::vector_type VectorType ;

  typedef Kokkos::Example::VectorImport<
     typename FixtureType::comm_list_type ,
     typename FixtureType::send_nodeid_type ,
     VectorType > ImportType ;

  //------------------------------------

  const unsigned newton_iteration_limit     = 10 ;
  const double   newton_iteration_tolerance = 1e-7 ;
  const unsigned cg_iteration_limit         = 200 ;
  const double   cg_iteration_tolerance     = 1e-7 ;

  //------------------------------------

  const int print_flag = use_print && std::is_same< Kokkos::HostSpace , typename Space::memory_space >::value ;

  int comm_rank ;
  int comm_size ;

  MPI_Comm_rank( comm , & comm_rank );
  MPI_Comm_size( comm , & comm_size );

  // Decompose by node to avoid mpi-communication for assembly

  const float bubble_x = 1.0 ;
  const float bubble_y = 1.0 ;
  const float bubble_z = 1.0 ;

  const FixtureType fixture( BoxElemPart::DecomposeNode , comm_size , comm_rank ,
                             use_elems[0] , use_elems[1] , use_elems[2] ,
                             bubble_x , bubble_y , bubble_z );


  {
    int global_error = ! fixture.ok();

#if defined( KOKKOS_ENABLE_MPI )
    int local_error = global_error ;
    global_error = 0 ;
    MPI_Allreduce( & local_error , & global_error , 1 , MPI_INT , MPI_SUM , comm );
#endif

    if ( global_error ) {
      throw std::runtime_error(std::string("Error generating finite element fixture"));
    }
  }

  //------------------------------------

  const ImportType comm_nodal_import(
    comm ,
    fixture.recv_node() ,
    fixture.send_node() ,
    fixture.send_nodeid() ,
    fixture.node_count_owned() ,
    fixture.node_count() - fixture.node_count_owned() );

  //------------------------------------

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;

  const Kokkos::Example::FENL::ManufacturedSolution
    manufactured_solution( 0 , 1 , bc_lower_value , bc_upper_value  );

  //------------------------------------

  for ( int k = 0 ; k < comm_size && use_print ; ++k ) {
    if ( k == comm_rank ) {
      typename FixtureType::node_grid_type::HostMirror
        h_node_grid = Kokkos::create_mirror_view( fixture.node_grid() );

      typename FixtureType::node_coord_type::HostMirror
        h_node_coord = Kokkos::create_mirror_view( fixture.node_coord() );

      typename FixtureType::elem_node_type::HostMirror
        h_elem_node = Kokkos::create_mirror_view( fixture.elem_node() );

      Kokkos::deep_copy( h_node_grid , fixture.node_grid() );
      Kokkos::deep_copy( h_node_coord , fixture.node_coord() );
      Kokkos::deep_copy( h_elem_node , fixture.elem_node() );

      std::cout << "MPI[" << comm_rank << "]" << std::endl ;
      std::cout << "Node grid {" ;
      for ( unsigned inode = 0 ; inode < fixture.node_count() ; ++inode ) {
        std::cout << " (" << h_node_grid(inode,0)
                  << "," << h_node_grid(inode,1)
                  << "," << h_node_grid(inode,2)
                  << ")" ;
      }
      std::cout << " }" << std::endl ;

      std::cout << "Node coord {" ;
      for ( unsigned inode = 0 ; inode < fixture.node_count() ; ++inode ) {
        std::cout << " (" << h_node_coord(inode,0)
                  << "," << h_node_coord(inode,1)
                  << "," << h_node_coord(inode,2)
                  << ")" ;
      }
      std::cout << " }" << std::endl ;

      std::cout << "Manufactured solution"
                << " a[" << manufactured_solution.a << "]"
                << " b[" << manufactured_solution.b << "]"
                << " K[" << manufactured_solution.K << "]"
                << " {" ;
      for ( unsigned inode = 0 ; inode < fixture.node_count() ; ++inode ) {
        std::cout << " " << manufactured_solution( h_node_coord( inode , 2 ) );
      }
      std::cout << " }" << std::endl ;

      std::cout << "ElemNode {" << std::endl ;
      for ( unsigned ielem = 0 ; ielem < fixture.elem_count() ; ++ielem ) {
        std::cout << "  elem[" << ielem << "]{" ;
        for ( unsigned inode = 0 ; inode < FixtureType::ElemNode ; ++inode ) {
          std::cout << " " << h_elem_node(ielem,inode);
        }
        std::cout << " }{" ;
        for ( unsigned inode = 0 ; inode < FixtureType::ElemNode ; ++inode ) {
          std::cout << " (" << h_node_grid(h_elem_node(ielem,inode),0)
                    << "," << h_node_grid(h_elem_node(ielem,inode),1)
                    << "," << h_node_grid(h_elem_node(ielem,inode),2)
                    << ")" ;
        }
        std::cout << " }" << std::endl ;
      }
      std::cout << "}" << std::endl ;
    }
    std::cout.flush();
    MPI_Barrier( comm );
  }

  //------------------------------------

  Kokkos::Timer wall_clock ;

  Perf perf_stats = Perf() ;

  for ( int itrial = 0 ; itrial < use_trials ; ++itrial ) {

    Perf perf = Perf() ;

    perf.global_elem_count = fixture.elem_count_global();
    perf.global_node_count = fixture.node_count_global();

    //----------------------------------
    // Create the sparse matrix graph and element-to-graph map
    // from the element->to->node identifier array.
    // The graph only has rows for the owned nodes.

    typename NodeNodeGraphType::Times graph_times;

    const NodeNodeGraphType
      mesh_to_graph( fixture.elem_node() , fixture.node_count_owned(), graph_times );

    perf.map_ratio          = maximum(comm, graph_times.ratio);
    perf.fill_node_set      = maximum(comm, graph_times.fill_node_set);
    perf.scan_node_count    = maximum(comm, graph_times.scan_node_count);
    perf.fill_graph_entries = maximum(comm, graph_times.fill_graph_entries);
    perf.sort_graph_entries = maximum(comm, graph_times.sort_graph_entries);
    perf.fill_element_graph = maximum(comm, graph_times.fill_element_graph);

    wall_clock.reset();
    // Create the sparse matrix from the graph:

    SparseMatrixType jacobian( mesh_to_graph.graph );

    Space::fence();

    perf.create_sparse_matrix = maximum( comm , wall_clock.seconds() );

    //----------------------------------

    for ( int k = 0 ; k < comm_size && print_flag ; ++k ) {
      if ( k == comm_rank ) {
        const unsigned nrow = jacobian.graph.numRows();
        std::cout << "MPI[" << comm_rank << "]" << std::endl ;
        std::cout << "JacobianGraph {" << std::endl ;
        for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
          std::cout << "  row[" << irow << "]{" ;
          const unsigned entry_end = jacobian.graph.row_map(irow+1);
          for ( unsigned entry = jacobian.graph.row_map(irow) ; entry < entry_end ; ++entry ) {
            std::cout << " " << jacobian.graph.entries(entry);
          }
          std::cout << " }" << std::endl ;
        }
        std::cout << "}" << std::endl ;

        std::cout << "ElemGraph {" << std::endl ;
        for ( unsigned ielem = 0 ; ielem < mesh_to_graph.elem_graph.extent(0) ; ++ielem ) {
          std::cout << "  elem[" << ielem << "]{" ;
          for ( unsigned irow = 0 ; irow < mesh_to_graph.elem_graph.extent(1) ; ++irow ) {
            std::cout << " {" ;
            for ( unsigned icol = 0 ; icol < mesh_to_graph.elem_graph.extent(2) ; ++icol ) {
              std::cout << " " << mesh_to_graph.elem_graph(ielem,irow,icol);
            }
            std::cout << " }" ;
          }
          std::cout << " }" << std::endl ;
        }
        std::cout << "}" << std::endl ;
      }
      std::cout.flush();
      MPI_Barrier( comm );
    }

    //----------------------------------

    // Allocate solution vector for each node in the mesh and residual vector for each owned node
    const VectorType nodal_solution( "nodal_solution" , fixture.node_count() );
    const VectorType nodal_residual( "nodal_residual" , fixture.node_count_owned() );
    const VectorType nodal_delta(    "nodal_delta" ,    fixture.node_count_owned() );

    // Create element computation functor
    const ElementComputationType elemcomp(
      use_atomic ? ElementComputationType( fixture , manufactured_solution.K , nodal_solution ,
                                           mesh_to_graph.elem_graph , jacobian , nodal_residual )
                 : ElementComputationType( fixture , manufactured_solution.K , nodal_solution ) );

    const NodeElemGatherFillType gatherfill(
      use_atomic ? NodeElemGatherFillType()
                 : NodeElemGatherFillType( fixture.elem_node() ,
                                           mesh_to_graph.elem_graph ,
                                           nodal_residual ,
                                           jacobian ,
                                           elemcomp.elem_residuals ,
                                           elemcomp.elem_jacobians ) );

    // Create boundary condition functor
    const DirichletComputationType dirichlet(
      fixture , nodal_solution , jacobian , nodal_residual ,
      2 /* apply at 'z' ends */ ,
      manufactured_solution.T_zmin ,
      manufactured_solution.T_zmax );

    //----------------------------------
    // Nonlinear Newton iteration:

    double residual_norm_init = 0 ;

    for ( perf.newton_iter_count = 0 ;
          perf.newton_iter_count < newton_iteration_limit ;
          ++perf.newton_iter_count ) {

      //--------------------------------

      comm_nodal_import( nodal_solution );

      //--------------------------------
      // Element contributions to residual and jacobian

      wall_clock.reset();

      Kokkos::deep_copy( nodal_residual , double(0) );
      Kokkos::deep_copy( jacobian.coeff , double(0) );

      elemcomp.apply();

      if ( ! use_atomic ) {
        gatherfill.apply();
      }

      Space::fence();
      perf.fill_time = maximum( comm , wall_clock.seconds() );

      //--------------------------------
      // Apply boundary conditions

      wall_clock.reset();

      dirichlet.apply();

      Space::fence();
      perf.bc_time = maximum( comm , wall_clock.seconds() );

      //--------------------------------
      // Evaluate convergence

      const double residual_norm =
        std::sqrt(
          Kokkos::Example::all_reduce(
            Kokkos::Example::dot( fixture.node_count_owned() , nodal_residual, nodal_residual ) , comm ) );

      perf.newton_residual = residual_norm ;

      if ( 0 == perf.newton_iter_count ) { residual_norm_init = residual_norm ; }

      if ( residual_norm < residual_norm_init * newton_iteration_tolerance ) { break ; }

      //--------------------------------
      // Solve for nonlinear update

      CGSolveResult cg_result ;

      Kokkos::Example::cgsolve( comm_nodal_import
                              , jacobian
                              , nodal_residual
                              , nodal_delta
                              , cg_iteration_limit
                              , cg_iteration_tolerance
                              , & cg_result
                              );

      // Update solution vector

      Kokkos::Example::waxpby( fixture.node_count_owned() , nodal_solution , -1.0 , nodal_delta , 1.0 , nodal_solution );

      perf.cg_iter_count += cg_result.iteration ;
      perf.matvec_time   += cg_result.matvec_time ;
      perf.cg_time       += cg_result.iter_time ;

      //--------------------------------

      if ( print_flag ) {
        const double delta_norm =
          std::sqrt(
            Kokkos::Example::all_reduce(
              Kokkos::Example::dot( fixture.node_count_owned() , nodal_delta, nodal_delta ) , comm ) );

        if ( 0 == comm_rank ) {
          std::cout << "Newton iteration[" << perf.newton_iter_count << "]"
                    << " residual[" << perf.newton_residual << "]"
                    << " update[" << delta_norm << "]"
                    << " cg_iteration[" << cg_result.iteration << "]"
                    << " cg_residual[" << cg_result.norm_res << "]"
                    << std::endl ;
        }

        for ( int k = 0 ; k < comm_size ; ++k ) {
          if ( k == comm_rank ) {
            const unsigned nrow = jacobian.graph.numRows();

            std::cout << "MPI[" << comm_rank << "]" << std::endl ;
            std::cout << "Residual {" ;
            for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
              std::cout << " " << nodal_residual(irow);
            }
            std::cout << " }" << std::endl ;

            std::cout << "Delta {" ;
            for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
              std::cout << " " << nodal_delta(irow);
            }
            std::cout << " }" << std::endl ;

            std::cout << "Solution {" ;
            for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
              std::cout << " " << nodal_solution(irow);
            }
            std::cout << " }" << std::endl ;

            std::cout << "Jacobian[ "
                      << jacobian.graph.numRows() << " x " << Kokkos::maximum_entry( jacobian.graph )
                      << " ] {" << std::endl ;
            for ( unsigned irow = 0 ; irow < nrow ; ++irow ) {
              std::cout << "  {" ;
              const unsigned entry_end = jacobian.graph.row_map(irow+1);
              for ( unsigned entry = jacobian.graph.row_map(irow) ; entry < entry_end ; ++entry ) {
                std::cout << " (" << jacobian.graph.entries(entry)
                          << "," << jacobian.coeff(entry)
                          << ")" ;
              }
              std::cout << " }" << std::endl ;
            }
            std::cout << "}" << std::endl ;
          }
          std::cout.flush();
          MPI_Barrier( comm );
        }
      }
      //--------------------------------
    }

    // Evaluate solution error

    if ( 0 == itrial ) {
      const typename FixtureType::node_coord_type::HostMirror
        h_node_coord = Kokkos::create_mirror_view( fixture.node_coord() );

      const typename VectorType::HostMirror
        h_nodal_solution = Kokkos::create_mirror_view( nodal_solution );

      Kokkos::deep_copy( h_node_coord , fixture.node_coord() );
      Kokkos::deep_copy( h_nodal_solution , nodal_solution );

      double error_max = 0 ;
      for ( unsigned inode = 0 ; inode < fixture.node_count_owned() ; ++inode ) {
        const double answer = manufactured_solution( h_node_coord( inode , 2 ) );
        const double error = ( h_nodal_solution(inode) - answer ) / answer ;
        if ( error_max < fabs( error ) ) { error_max = fabs( error ); }
      }

      perf.error_max = std::sqrt( Kokkos::Example::all_reduce_max( error_max , comm ) );

      perf_stats = perf ;
    }
    else {
      perf_stats.fill_node_set = std::min( perf_stats.fill_node_set , perf.fill_node_set );
      perf_stats.scan_node_count = std::min( perf_stats.scan_node_count , perf.scan_node_count );
      perf_stats.fill_graph_entries = std::min( perf_stats.fill_graph_entries , perf.fill_graph_entries );
      perf_stats.sort_graph_entries = std::min( perf_stats.sort_graph_entries , perf.sort_graph_entries );
      perf_stats.fill_element_graph = std::min( perf_stats.fill_element_graph , perf.fill_element_graph );
      perf_stats.create_sparse_matrix = std::min( perf_stats.create_sparse_matrix , perf.create_sparse_matrix );
      perf_stats.fill_time = std::min( perf_stats.fill_time , perf.fill_time );
      perf_stats.bc_time = std::min( perf_stats.bc_time , perf.bc_time );
      perf_stats.cg_time = std::min( perf_stats.cg_time , perf.cg_time );
    }
  }

  return perf_stats ;
}

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_IMPL_HPP */

