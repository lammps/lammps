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

#include <utility>
#include <iostream>

#include <Kokkos_Core.hpp>

#include <BoxElemPart.hpp>

namespace Kokkos {
namespace Example {
template< class > void test_fixture();
}
}

int test_box( const size_t global_size
            , const size_t global_box[][2]
            , const bool print_verbose )
{
  size_t global_count = 0 ;
  size_t global_max = 0 ;
  size_t global_min = Kokkos::Example::box_count( global_box );
  size_t global_box_max[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };
  size_t global_box_min[3][2] = { { 0 , global_box[0][1] } , { 0 , global_box[1][1] } , { 0 , global_box[2][1] } };
  size_t intersect_error = 0 ;
  size_t neighbor_max = 0 ;

  for ( size_t global_rank = 0 ; global_rank < global_size ; ++global_rank ) {
    size_t box[3][2] = { { 0 , global_box[0][1] } , { 0 , global_box[1][1] } , { 0 , global_box[2][1] } };
    size_t ghost_box[3][2] ;
    size_t neighbor_count = 0 ;

    Kokkos::Example::box_partition( global_size , global_rank , global_box , box );

    Kokkos::Example::box_ghost_layer( global_box , box , 1 , ghost_box );

    {
      const size_t n = Kokkos::Example::box_count( box );

      for ( int i = 0 ; i < 3 ; ++i ) {
        if ( ( box[i][1] - box[i][0] ) < ( global_box_min[i][1] - global_box_min[i][0] ) ) {
          global_box_min[i][0] = box[i][0] ;
          global_box_min[i][1] = box[i][1] ;
        }
        if ( ( box[i][1] - box[i][0] ) > ( global_box_max[i][1] - global_box_max[i][0] ) ) {
          global_box_max[i][0] = box[i][0] ;
          global_box_max[i][1] = box[i][1] ;
        }
      }

      global_max = std::max( global_max , n );
      global_min = std::min( global_min , n );
      global_count += n ;
    }

    for ( size_t other_rank = 0 ; other_rank  < global_size ; ++other_rank ) {

      if ( other_rank == global_rank ) continue ;

      size_t other_box[3][2] = { { 0 , global_box[0][1] } , { 0 , global_box[1][1] } , { 0 , global_box[2][1] } };
      size_t intersect_box[3][2] ;

      Kokkos::Example::box_partition( global_size , other_rank , global_box , other_box );

      Kokkos::Example::box_intersect( intersect_box , box , other_box );

      const size_t n = Kokkos::Example::box_count( intersect_box );

      intersect_error += n ;

      Kokkos::Example::box_intersect( intersect_box , ghost_box , other_box );

      neighbor_count += Kokkos::Example::box_count( intersect_box ) ? 1 : 0 ;

      if ( n ) {
        std::cout << "box partition intersection error" << std::endl ;
        std::cout << "box = {"
                  << " [ " << box[0][0] << " , " << box[0][1] << " )"
                  << " [ " << box[1][0] << " , " << box[1][1] << " )"
                  << " [ " << box[2][0] << " , " << box[2][1] << " )"
                  << " }" << std::endl ;
        std::cout << "other_box = {"
                  << " [ " << other_box[0][0] << " , " << other_box[0][1] << " )"
                  << " [ " << other_box[1][0] << " , " << other_box[1][1] << " )"
                  << " [ " << other_box[2][0] << " , " << other_box[2][1] << " )"
                  << " }" << std::endl ;
        return 0 ;
      }
    }

    neighbor_max = std::max( neighbor_max , neighbor_count );
  }

  if ( print_verbose ) {

    std::cout << "global_part = " << global_size << std::endl ;
    std::cout << "global_box  = { "
              << " [ " << global_box[0][0] << " .. " << global_box[0][1] << " ) X"
              << " [ " << global_box[1][0] << " .. " << global_box[1][1] << " ) X"
              << " [ " << global_box[2][0] << " .. " << global_box[2][1] << " )"
              << " }" << std::endl ;
    std::cout << "count( global_box ) = " << Kokkos::Example::box_count( global_box ) << std::endl ;
    std::cout << "sum partition( global_box ) = " << global_count << std::endl ;
    std::cout << "avg partition( global_box ) = " << size_t( double(global_count) / double(global_size)) << std::endl ;
    std::cout << "min partition( global_box ) = " << global_min << std::endl ;
    std::cout << "min part X   ( global_box ) = [ " << global_box_min[0][0] << " .. " << global_box_min[0][1] << " )" << std::endl ;
    std::cout << "min part Y   ( global_box ) = [ " << global_box_min[1][0] << " .. " << global_box_min[1][1] << " )" << std::endl ;
    std::cout << "min part Z   ( global_box ) = [ " << global_box_min[2][0] << " .. " << global_box_min[2][1] << " )" << std::endl ;
    std::cout << "max partition( global_box ) = " << global_max << std::endl ;
    std::cout << "max part X   ( global_box ) = [ " << global_box_max[0][0] << " .. " << global_box_max[0][1] << " )" << std::endl ;
    std::cout << "max part Y   ( global_box ) = [ " << global_box_max[1][0] << " .. " << global_box_max[1][1] << " )" << std::endl ;
    std::cout << "max part Z   ( global_box ) = [ " << global_box_max[2][0] << " .. " << global_box_max[2][1] << " )" << std::endl ;
    std::cout << "sum intersect( global_box ) = " << intersect_error << std::endl ;
    std::cout << "max neighbor = " << neighbor_max << std::endl ;
  }

  return neighbor_max ;
}

void test_elem()
{
  const Kokkos::Example::BoxElemPart::Decompose
    decompose = Kokkos::Example::BoxElemPart:: DecomposeElem ; // DecomposeElem | DecomposeNode ;
  const size_t global_size = 256 ;
  const size_t global_nx = 100 ;
  const size_t global_ny = 120 ;
  const size_t global_nz = 140 ;

  double node_count_avg = 0 ;
  size_t node_count_max = 0 ;
  size_t node_count_min = ( global_nx + 1 ) * ( global_ny + 1 ) * ( global_nz + 1 );
  double elem_count_avg = 0 ;
  size_t elem_count_max = 0 ;
  size_t elem_count_min = global_nx * global_ny * global_nz ;
  double recv_count_avg = 0 ;
  size_t recv_count_max = 0 ;
  size_t recv_count_min = global_size ;
  double send_count_avg = 0 ;
  size_t send_count_max = 0 ;
  size_t send_count_min = global_size ;

  for ( size_t r = 0 ; r < global_size ; ++r ) {
    const Kokkos::Example::BoxElemPart
       fixture( Kokkos::Example::BoxElemPart::ElemLinear ,
                decompose , global_size , r , global_nx , global_ny , global_nz );

    // Print a sample:

    // if ( r == global_size * 2 / 3 ) fixture.print( std::cout );

    // Verify recv/send alignment:

    {
      size_t recv_lid = fixture.owns_node_count();

      for ( size_t i = 0 ; i < fixture.recv_node_msg_count() ; ++i ) {
        const size_t recv_rank  = fixture.recv_node_rank( i );
        const size_t recv_count = fixture.recv_node_count( i );

        const Kokkos::Example::BoxElemPart other_fixture(
           Kokkos::Example::BoxElemPart::ElemLinear ,
           decompose , global_size , recv_rank , global_nx , global_ny , global_nz );

        size_t send_item = 0 ;

        size_t j = 0 ;
        while ( j < other_fixture.send_node_msg_count() && other_fixture.send_node_rank(j) != r ) {
          send_item += other_fixture.send_node_count( j );
           ++j ;
        }

        if ( recv_count != other_fixture.send_node_count(j) ) {
          std::cout << "Error P[" << r << "].recv(" << recv_count << ") != "
                    << "P[" << recv_rank << "].send(" << other_fixture.send_node_count(j) << ")"
                    << std::endl ;
        }
        else {

          for ( size_t k = 0 ; k < recv_count ; ++k , ++send_item , ++recv_lid ) {

            const size_t send_lid = other_fixture.send_node_id( send_item );

            size_t recv_coord[3] , send_coord[3] ;

            fixture.local_node_coord( recv_lid , recv_coord );

            other_fixture.local_node_coord( send_lid , send_coord );

            if ( recv_coord[0] != send_coord[0] ||
                 recv_coord[1] != send_coord[1] ||
                 recv_coord[2] != send_coord[2] ) {
              std::cout << "Error P[" << r << "].recv[" << recv_lid << "]{ "
                        << recv_coord[0] << " , "
                        << recv_coord[1] << " , "
                        << recv_coord[2] << " } != "
                        << "P[" << recv_rank << "].send[" << send_lid << "]{ "
                        << send_coord[0] << " , "
                        << send_coord[1] << " , "
                        << send_coord[2] << " }"
                        << std::endl ;
            }
          }
        }
      }
    }

    node_count_avg += fixture.owns_node_count();
    elem_count_avg += fixture.uses_elem_count();
    recv_count_avg += fixture.recv_node_msg_count();
    send_count_avg += fixture.send_node_msg_count();

    elem_count_min = std::min( (size_t) fixture.uses_elem_count() , elem_count_min );
    elem_count_max = std::max( (size_t) fixture.uses_elem_count() , elem_count_max );
    node_count_min = std::min( (size_t) fixture.owns_node_count() , node_count_min );
    node_count_max = std::max( (size_t) fixture.owns_node_count() , node_count_max );

    recv_count_max = std::max( (size_t) fixture.recv_node_msg_count() , recv_count_max );
    recv_count_min = std::min( (size_t) fixture.recv_node_msg_count() , recv_count_min );
    send_count_max = std::max( (size_t) fixture.send_node_msg_count() , send_count_max );
    send_count_min = std::min( (size_t) fixture.send_node_msg_count() , send_count_min );
  }

  node_count_avg /= double(global_size);
  elem_count_avg /= double(global_size);
  recv_count_avg /= double(global_size);
  send_count_avg /= double(global_size);

  std::cout << "Elem min(" << elem_count_min << ") avg(" << elem_count_avg << ") max(" << elem_count_max << ") " << std::endl
            << "Node min(" << node_count_min << ") avg(" << node_count_avg << ") max(" << node_count_max << ") " << std::endl
            << "Recv min(" << recv_count_min << ") avg(" << recv_count_avg << ") max(" << recv_count_max << ") " << std::endl
            << "Send min(" << send_count_min << ") avg(" << send_count_avg << ") max(" << send_count_max << ") " << std::endl
            ;
}

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc,argv);
  for ( int i = 1 ; i <= 32 ; ++i ) {
    const size_t global_size = 16 * i ;
    const size_t global_box[3][2] = { { 0 , 65 } , { 0 , 65 } , { 0 , 65 } };
    if ( 30 < test_box( global_size , global_box , false ) ) {
      test_box( global_size , global_box , true );
    }
  }

//  test_elem();

  {
    std::cout << "test_fixture< Host >" << std::endl ;
    Kokkos::Example::test_fixture< Kokkos::DefaultHostExecutionSpace >();
  }

#if defined( KOKKOS_ENABLE_CUDA )
  {
    std::cout << "test_fixture< Cuda >" << std::endl ;
    Kokkos::Example::test_fixture< Kokkos::Cuda >();
  }
#endif

#if defined( KOKKOS_ENABLE_ROCM )
  {
    std::cout << "test_fixture< ROCm >" << std::endl ;
    Kokkos::Example::test_fixture< Kokkos::Experimental::ROCm >();
  }
#endif
  Kokkos::finalize();
}

