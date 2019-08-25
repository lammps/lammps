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

#ifndef TESTFEMESHBOXFIXTURE_HPP
#define TESTFEMESHBOXFIXTURE_HPP

#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshFixture.hpp>

#include <ParallelComm.hpp>

//----------------------------------------------------------------------------

namespace TestFEMesh {

template< class ViewType >
struct VerifyUnpack  ;

template< typename DeviceType, typename T >
struct VerifyUnpack< Kokkos::View< T*[3] , DeviceType > >
{
  typedef DeviceType     execution_space ;
  typedef typename execution_space::size_type  size_type ;
  typedef size_type               value_type ;

  typedef Kokkos::View< T* ,    execution_space > buffer_type ;
  typedef Kokkos::View< T*[3] , execution_space > array_type ;

private:

  array_type  node_coords ;
  buffer_type buffer ;
  size_type   node_begin ;

public:

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i , value_type & update ) const
  {
    const size_type node_id = i + node_begin ;
    const size_type k = i * 3 ;

    const long xb = buffer[k];
    const long yb = buffer[k+1];
    const long zb = buffer[k+2];
    const long xn = node_coords(node_id,0);
    const long yn = node_coords(node_id,1);
    const long zn = node_coords(node_id,2);

    if ( xb != xn || yb != yn || zb != zn ) {
      printf("TestFEMesh::VerifyUnpack failed at %d : node %d : { %ld %ld %ld } != { %ld %ld %ld }\n",
             (int)i,(int)node_id, xb,yb,zb, xn, yn, zn );
      ++update ;
    }
  }

  static inline
  size_type unpack( const array_type  & arg_node_coords ,
                    const size_type     arg_node_begin ,
                    const size_type     arg_node_count ,
                    const buffer_type & arg_buffer )
  {
    VerifyUnpack op ;
    op.node_coords = arg_node_coords ;
    op.buffer      = arg_buffer ;
    op.node_begin  = arg_node_begin ;
    size_type count = 0 ;
    Kokkos::parallel_reduce( arg_node_count , op , count );
    return count ;
  }
};

}

//----------------------------------------------------------------------------

#ifdef KOKKOS_ENABLE_MPI

namespace TestFEMesh {

template< typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  const HybridFEM::FEMesh< coordinate_scalar_type ,
                           ElemNodeCount ,
                           Device > & mesh )
{
  typedef HybridFEM::FEMesh< coordinate_scalar_type, ElemNodeCount, Device > femesh_type ;
  typedef typename femesh_type::node_coords_type node_coords_type ;

  comm::Machine machine = mesh.parallel_data_map.machine ;

  // Communicate node coordinates to verify communication and setup.

  const size_t chunk_size = 3 ;

  Kokkos::AsyncExchange< coordinate_scalar_type, Device, Kokkos::ParallelDataMap >
    exchange( mesh.parallel_data_map , chunk_size );

  const size_t send_begin = mesh.parallel_data_map.count_interior ;
  const size_t send_count = mesh.parallel_data_map.count_send ;

  const size_t recv_begin = mesh.parallel_data_map.count_owned ;
  const size_t recv_count = mesh.parallel_data_map.count_receive ;

  typedef Kokkos::PackArray< node_coords_type > pack_type ;

  pack_type::pack( exchange.buffer(), send_begin, send_count, mesh.node_coords );

  exchange.setup();

  // Launch local-action device kernels

  exchange.send_receive();

  unsigned long local[3] ;
  local[0] = mesh.parallel_data_map.count_owned ;
  local[1] = mesh.parallel_data_map.count_receive ;
  local[2] = TestFEMesh::VerifyUnpack< node_coords_type >::unpack( mesh.node_coords, recv_begin, recv_count, exchange.buffer() );

  unsigned long global[3] = { 0 , 0 , 0 };

  MPI_Allreduce( local , global ,
                 3 , MPI_UNSIGNED_LONG , MPI_SUM , machine.mpi_comm );

  if ( 0 == comm::rank( machine ) ) {
    std::cout << ( global[2] ? "FAILED" : "PASSED" )
              << ": TestFEMesh::verify_parallel "
              << "NP(" << comm::size( machine )
              << ") total_node(" << global[0]
              << ") verified_nodes(" << global[1]
              << ") failed_nodes(" << global[2]
              << ")" << std::endl ;
  }
}

} // namespace TestFEMesh

#else /* ! #ifdef KOKKOS_ENABLE_MPI */

namespace TestFEMesh {

template< typename coordinate_scalar_type ,
          unsigned ElemNodeCount ,
          class Device >
void verify_parallel(
  const HybridFEM::FEMesh< coordinate_scalar_type ,
                           ElemNodeCount ,
                           Device > & )
{}

} // namespace TestFEMesh

#endif /* ! #ifdef KOKKOS_ENABLE_MPI */

//----------------------------------------------------------------------------

template< class Device >
void test_box_fixture( comm::Machine machine ,
                       const size_t gang_count ,
                       const size_t nodes_nx ,
                       const size_t nodes_ny ,
                       const size_t nodes_nz )
{
  typedef long                coordinate_scalar_type ;
  typedef FixtureElementHex8  fixture_element_type ;

  typedef BoxMeshFixture< coordinate_scalar_type ,
                          Device ,
                          fixture_element_type > fixture_type ;

  typedef typename fixture_type::FEMeshType  mesh_type ;

  const size_t proc_count = comm::size( machine );
  const size_t proc_local = comm::rank( machine ) ;

  mesh_type mesh =
    fixture_type::create( proc_count , proc_local , gang_count ,
                          nodes_nx - 1 , nodes_ny - 1 , nodes_nz - 1 );

  mesh.parallel_data_map.machine = machine ;

  TestFEMesh::verify_parallel( mesh );
}

#endif /* #ifndef TESTFEMESHBOXFIXTURE_HPP */


