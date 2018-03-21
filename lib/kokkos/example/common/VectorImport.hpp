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

#ifndef KOKKOS_VECTORIMPORT_HPP
#define KOKKOS_VECTORIMPORT_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include <WrapMPI.hpp>

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
struct VectorImport ;

} // namespace Example
} // namespace Kokkos

#if ! defined( KOKKOS_ENABLE_MPI )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
struct VectorImport {

  const MPI_Comm comm ;
  const unsigned count_owned ;
  const unsigned count_receive ;

  VectorImport( MPI_Comm arg_comm ,
                const CommMessageType & ,
                const CommMessageType & ,
                const CommIdentType   & ,
                const unsigned arg_count_owned ,
                const unsigned arg_count_receive )
    : comm( arg_comm )
    , count_owned( arg_count_owned )
    , count_receive( arg_count_receive )
    {}

  inline
  void operator()( const VectorType & ) const {}
};


} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else /* defined( KOKKOS_ENABLE_MPI ) */

namespace Kokkos {
namespace Example {

template< class CommMessageType , class CommIdentType , class VectorType >
class VectorImport {
private:

  // rank == 1 or array_layout == LayoutRight
  enum { OK = Kokkos::Impl::StaticAssert<
           ( VectorType::rank == 1 ) ||
           std::is_same< typename VectorType::array_layout , Kokkos::LayoutRight >::value
         >::value };

  typedef typename VectorType::HostMirror HostVectorType ;

  enum { ReceiveInPlace =
    std::is_same< typename VectorType::memory_space ,
                           typename HostVectorType::memory_space >::value };

  const CommMessageType  recv_msg ;
  const CommMessageType  send_msg ;
  const CommIdentType    send_nodeid ;
  VectorType             send_buffer ;
  HostVectorType         host_send_buffer ;
  HostVectorType         host_recv_buffer ;
  unsigned               chunk ;

public:

  const MPI_Comm         comm ;
  const unsigned         count_owned ;
  const unsigned         count_receive ;

  struct Pack {
    typedef typename VectorType::execution_space execution_space ;
    const CommIdentType  index ;
    const VectorType     source ;
    const VectorType     buffer ;

    KOKKOS_INLINE_FUNCTION
    void operator()( const unsigned i ) const
      { buffer( i ) = source( index(i) ); }

    Pack( const CommIdentType  & arg_index ,
          const VectorType     & arg_source ,
          const VectorType     & arg_buffer )
      : index( arg_index )
      , source( arg_source )
      , buffer( arg_buffer )
    {
      Kokkos::parallel_for( index.dimension_0() , *this );
      execution_space::fence();
    }
  };

  VectorImport( MPI_Comm arg_comm ,
                const CommMessageType & arg_recv_msg ,
                const CommMessageType & arg_send_msg ,
                const CommIdentType   & arg_send_nodeid ,
                const unsigned arg_count_owned ,
                const unsigned arg_count_receive )
    : recv_msg( arg_recv_msg )
    , send_msg( arg_send_msg )
    , send_nodeid( arg_send_nodeid )
    , send_buffer()
    , host_send_buffer()
    , host_recv_buffer()
    , comm( arg_comm )
    , count_owned( arg_count_owned )
    , count_receive( arg_count_receive )
    {
      if ( ! ReceiveInPlace ) {
        host_recv_buffer = HostVectorType("recv_buffer",count_receive);
      }

      unsigned send_count = 0 ;
      for ( unsigned i = 0 ; i < send_msg.dimension_0() ; ++i ) { send_count += send_msg(i,1); }
      send_buffer      = VectorType("send_buffer",send_count);
      host_send_buffer = Kokkos::create_mirror_view( send_buffer );
    }

  inline
  void operator()( const VectorType & v ) const
  {
    typedef typename VectorType::value_type  scalar_type ;

    const int mpi_tag = 42 ;
    const unsigned chunk = v.dimension_1();

    // Subvector for receives
    const std::pair<unsigned,unsigned> recv_range( count_owned , count_owned + count_receive );
    const VectorType recv_vector = Kokkos::subview( v , recv_range );

    std::vector< MPI_Request > recv_request( recv_msg.dimension_0() , MPI_REQUEST_NULL );

    { // Post receives
      scalar_type * ptr =
        ReceiveInPlace ? recv_vector.ptr_on_device() : host_recv_buffer.ptr_on_device();

      for ( size_t i = 0 ; i < recv_msg.dimension_0() ; ++i ) {
        const int proc  = recv_msg(i,0);
        const int count = recv_msg(i,1) * chunk ;

        MPI_Irecv( ptr , count * sizeof(scalar_type) , MPI_BYTE ,
                   proc , mpi_tag , comm , & recv_request[i] );

        ptr += count ;
      }
    }

    MPI_Barrier( comm );

    { // Pack and send 
      const Pack pack( send_nodeid , v , send_buffer );

      Kokkos::deep_copy( host_send_buffer , send_buffer );

      scalar_type * ptr = host_send_buffer.ptr_on_device();

      for ( size_t i = 0 ; i < send_msg.dimension_0() ; ++i ) {
        const int proc  = send_msg(i,0);
        const int count = send_msg(i,1) * chunk ;

        // MPI_Ssend blocks until
        // (1) a receive is matched for the message and
        // (2) the send buffer can be re-used.
        //
        // It is suggested that MPI_Ssend will have the best performance:
        // http://www.mcs.anl.gov/research/projects/mpi/sendmode.html .

        MPI_Ssend( ptr ,
                   count * sizeof(scalar_type) , MPI_BYTE ,
                   proc , mpi_tag , comm );

        ptr += count ;
      }
    }

    // Wait for receives and verify:

    for ( size_t i = 0 ; i < recv_msg.dimension_0() ; ++i ) {
      MPI_Status recv_status ;
      int recv_which = 0 ;
      int recv_size  = 0 ;

      MPI_Waitany( recv_msg.dimension_0() , & recv_request[0] , & recv_which , & recv_status );

      const int recv_proc = recv_status.MPI_SOURCE ;

      MPI_Get_count( & recv_status , MPI_BYTE , & recv_size );

      // Verify message properly received:

      const int  expected_proc = recv_msg(recv_which,0);
      const int  expected_size = recv_msg(recv_which,1) * chunk * sizeof(scalar_type);

      if ( ( expected_proc != recv_proc ) ||
           ( expected_size != recv_size ) ) {

        int local_rank  = 0 ;

        MPI_Comm_rank( comm , & local_rank );

        std::ostringstream msg ;
        msg << "VectorImport error:"
            << " P" << local_rank
            << " received from P" << recv_proc
            << " size "     << recv_size
            << " expected " << expected_size
            << " from P"    << expected_proc ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Copy received data to device memory.

    if ( ! ReceiveInPlace ) { Kokkos::deep_copy( recv_vector , host_recv_buffer ); }
  }
};

} // namespace Example
} // namespace Kokkos

#endif

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VECTORIMPORT_HPP */


