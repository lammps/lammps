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

#ifndef KOKKOS_PARALLELDATAMAP_HPP
#define KOKKOS_PARALLELDATAMAP_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <ParallelComm.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Parallel distributed data mapping
 *
 *  ordering { interior : { owned items not sent elsewhere }
 *             send     : { owned items sent }
 *             receive  : { not-owned items received } }
 *
 *  recv { { N ghosted items from process P : ( P , N ) } }
 *
 *  send { { N send items to process P : ( P , N ) } }
 *
 *  send_item { send item offsets within 'send' range }
 */
struct ParallelDataMap {
  typedef View< unsigned*[2], HostSpace >  host_recv_type ;
  typedef View< unsigned*[2], HostSpace >  host_send_type ;
  typedef View< unsigned* ,   HostSpace >  host_send_item_type ;

  comm::Machine        machine ;
  host_recv_type       host_recv ;
  host_send_type       host_send ;
  host_send_item_type  host_send_item ;
  unsigned             count_interior ;
  unsigned             count_send ;
  unsigned             count_owned ; // = count_interior + count_send
  unsigned             count_receive ;

  void assign( const unsigned arg_count_interior ,
               const unsigned arg_count_owned ,
               const unsigned arg_count_total ,
               const unsigned arg_recv_msg ,
               const unsigned arg_send_msg ,
               const unsigned arg_send_count )
  {
    const std::string label("Kokkos::ParallelDataMap buffer");

    count_interior = arg_count_interior ;
    count_owned    = arg_count_owned ;
    count_send     = arg_count_owned - arg_count_interior ;
    count_receive  = arg_count_total - arg_count_owned ;

    host_recv = host_recv_type( label , arg_recv_msg );
    host_send = host_send_type( label , arg_send_msg );
    host_send_item = host_send_item_type( label , arg_send_count );
  }
};

//----------------------------------------------------------------------------
//PackArray
//----------------------------------------------------------------------------
template< class ArrayType , class Rank = void >
struct PackArray ;

template< typename DeviceType, typename ValueType >
struct PackArray< View< ValueType* , DeviceType > , void >
{
  typedef DeviceType                         execution_space ;
  typedef typename DeviceType::size_type     size_type ;
  typedef View< ValueType* , execution_space >  array_type ;
  typedef View< ValueType* , execution_space >  buffer_type ;

private:

  buffer_type  output ;
  array_type   input ;
  size_type    base ;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  { output[i] = input(base+i); }

  inline
  static
  void pack( const buffer_type & arg_output ,
             const size_type     arg_begin ,
             const size_type     arg_count ,
             const array_type  & arg_input )
  {
    PackArray op ;
    op.output = arg_output ;
    op.input  = arg_input ;
    op.base   = arg_begin ;
    parallel_for( arg_count , op );
  }
};

template< typename DeviceType, typename ValueType , unsigned N1 >
struct PackArray< View< ValueType*[N1] , DeviceType > , void >
{
  typedef DeviceType                                  execution_space ;
  typedef typename DeviceType::size_type              size_type ;
  typedef View< ValueType*[N1] , execution_space >       array_type ;
  typedef View< ValueType* , execution_space >           buffer_type ;

private:

  buffer_type  output ;
  array_type   input ;
  size_type    base ;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
    for ( size_type j = 0 , k = i * N1 ; j < N1 ; ++j , ++k ) {
      output[k] = input(base+i,j);
    }
  }

  inline static
  void pack( const buffer_type & arg_output ,
             const size_type     arg_begin ,
             const size_type     arg_count ,
             const array_type  & arg_input )
  {
    if ( arg_count ) {
      PackArray op ;
      op.output = arg_output ;
      op.input  = arg_input ;
      op.base   = arg_begin ;
      parallel_for( arg_count , op );
    }
  }
};

//----------------------------------------------------------------------------
//UnpackArray
//----------------------------------------------------------------------------
template< class ArrayType , class Rank = void > struct UnpackArray ;

template< typename DeviceType, typename ValueType >
struct UnpackArray< View< ValueType* , DeviceType > , void >
{
  typedef DeviceType                         execution_space ;
  typedef typename DeviceType::size_type     size_type ;
  typedef View< ValueType* , execution_space >  array_type ;
  typedef View< ValueType* , execution_space >  buffer_type ;

private:

  array_type   output ;
  buffer_type  input ;
  size_type    base ;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  { output(base+i) = input[i]; }

  inline
  static
  void unpack( const array_type  & arg_output ,
               const buffer_type & arg_input ,
               const size_type     arg_begin ,
               const size_type     arg_count )
  {
    UnpackArray op ;
    op.output = arg_output ;
    op.input  = arg_input ;
    op.base   = arg_begin ;
    parallel_for( arg_count , op );
  }
};

template< typename DeviceType, typename ValueType , unsigned N1 >
struct UnpackArray< View< ValueType*[N1] , DeviceType > , void >
{
  typedef DeviceType                                  execution_space ;
  typedef typename DeviceType::size_type              size_type ;
  typedef View< ValueType* , execution_space >           buffer_type ;
  typedef View< ValueType*[N1] , execution_space >       array_type ;

private:

  array_type   output ;
  buffer_type  input ;
  size_type    base ;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
    for ( size_type j = 0 , k = i * N1 ; j < N1 ; ++j , ++k ) {
      output(base+i,j) = input(k);
    }
  }

  inline
  static
  void unpack( const array_type  & arg_output ,
               const buffer_type & arg_input ,
               const size_type     arg_begin ,
               const size_type     arg_count )
  {
    if ( arg_count ) {
      UnpackArray op ;
      op.output = arg_output ;
      op.input  = arg_input ;
      op.base   = arg_begin ;
      parallel_for( arg_count , op );
    }
  }
};
//----------------------------------------------------------------------------
template< class ValueType , class Device , class DataMap >
class AsyncExchange ;

} // namespace Kokkos

//----------------------------------------------------------------------------
// Application call procedure:
//
// construct: AsyncExchange object
// * pack send buffer on device
// initiate: copy send buffer from device to host
// * dispatch asynchronous local work
// complete: send/receive on host, copy receive buffer to device
// * unpack receive buffer on device
// destroy: AsyncExchange object
//
//----------------------------------------------------------------------------

#ifdef KOKKOS_ENABLE_MPI

namespace Kokkos {

template< class ValueType , class Device >
class AsyncExchange< ValueType, Device , Kokkos::ParallelDataMap > {
public:

  typedef Device                                    execution_space ;
  typedef Kokkos::ParallelDataMap                   data_map_type ;
  typedef Kokkos::View< ValueType* , execution_space >  buffer_dev_type ;
  typedef typename buffer_dev_type::HostMirror      buffer_host_type ;

private:

  static const int mpi_tag = 11 ;

  const data_map_type  data_map ;
  unsigned             chunk_size ;
  unsigned             send_count_max ;
  buffer_host_type     host_recv_buffer ;
  buffer_host_type     host_send_buffer ;
  buffer_host_type     send_msg_buffer ;
  buffer_dev_type      dev_buffer ;
  buffer_dev_type      dev_send_buffer ; // Subview for send
  buffer_dev_type      dev_recv_buffer ; // Subview for receive
  std::vector< MPI_Request > recv_request ;

public:

  const buffer_dev_type & buffer() const { return dev_buffer ; }

  AsyncExchange( const data_map_type & arg_data_map ,
                 const size_t          arg_chunk )
  : data_map( arg_data_map )
  , chunk_size( arg_chunk )
  , send_count_max( 0 )
  , host_recv_buffer()
  , host_send_buffer()
  , send_msg_buffer()
  , dev_buffer()
  , dev_send_buffer()
  , dev_recv_buffer()
  , recv_request()
  {
    const size_t send_msg_count = arg_data_map.host_send.dimension_0();
    const size_t recv_msg_count = arg_data_map.host_recv.dimension_0();

    const size_t send_msg_length = arg_chunk * arg_data_map.count_send ;
    const size_t recv_msg_length = arg_chunk * arg_data_map.count_receive ;

    for ( size_t i = 0 ; i < send_msg_count ; ++i ) {
      send_count_max = std::max( send_count_max ,
                                 (unsigned) arg_data_map.host_send(i,1) );
    }

    // A single shared buffer on the device can be used for
    // send and receive message buffers.
    dev_buffer = buffer_dev_type(
                     std::string("AsyncExchange dev_buffer") ,
                     std::max( send_msg_length , recv_msg_length ) );

    // Total send subview of the device buffer
    dev_send_buffer =
      Kokkos::subview( dev_buffer , std::pair<size_t,size_t>( 0 , send_msg_length ) );

    // Total receive subview of the device buffer
    dev_recv_buffer =
      Kokkos::subview( dev_buffer , std::pair<size_t,size_t>( 0 , recv_msg_length ) );

    // Total receive message buffer on the host:
    host_recv_buffer = buffer_host_type(
                           std::string("AsyncExchange host_recv_buffer") ,
                           recv_msg_length );

    // Total send message buffer on the host:
    host_send_buffer = buffer_host_type(
                           std::string("AsyncExchange host_send_buffer") ,
                           send_msg_length );

    // Individual send message buffer on the host:
    send_msg_buffer = buffer_host_type(
                          std::string("AsyncExchange send_msg_buffer") ,
                          arg_chunk * send_count_max );

    // MPI asynchronous receive request handles:
    recv_request.assign( recv_msg_count , MPI_REQUEST_NULL );
  }

  //------------------------------------------------------------------------

  void setup()
  {
    { // Post receives:
      const size_t recv_msg_count = data_map.host_recv.dimension_0();

      ValueType * ptr = host_recv_buffer.ptr_on_device();

      for ( size_t i = 0 ; i < recv_msg_count ; ++i ) {
        const int proc  = data_map.host_recv(i,0);
        const int count = data_map.host_recv(i,1) * chunk_size ;

        MPI_Irecv( ptr , count * sizeof(ValueType) , MPI_BYTE ,
                   proc , mpi_tag , data_map.machine.mpi_comm ,
                   & recv_request[i] );

        ptr += count ;
      }
    }

    // Copy send buffer from the device to host memory for sending

    Kokkos::deep_copy( host_send_buffer , dev_send_buffer );

    // Done with the device until communication is complete.
    // Application can dispatch asynchronous work on the device.
  }

  // Application can dispatch local work to device ...
  // No communication progress until main thread calls 'send_receive'

  void send_receive()
  {
    const size_t recv_msg_count = data_map.host_recv.dimension_0();
    const size_t send_msg_count = data_map.host_send.dimension_0();

    // Pack and send:

    for ( size_t i = 0 , j = 0 ; i < send_msg_count ; ++i ) {
      const int proc  = data_map.host_send(i,0);
      const int count = data_map.host_send(i,1);

      for ( int k = 0 , km = 0 ; k < count ; ++k , ++j ) {
        const int km_end = km + chunk_size ;
        int ki = chunk_size * data_map.host_send_item(j);

        for ( ; km < km_end ; ++km , ++ki ) {
          send_msg_buffer[km] = host_send_buffer[ki];
        }
      }

      // MPI_Ssend blocks until
      // (1) a receive is matched for the message and
      // (2) the send buffer can be re-used.
      //
      // It is suggested that MPI_Ssend will have the best performance:
      // http://www.mcs.anl.gov/research/projects/mpi/sendmode.html .

      MPI_Ssend( send_msg_buffer.ptr_on_device(),
                 count * chunk_size * sizeof(ValueType) , MPI_BYTE ,
                 proc , mpi_tag , data_map.machine.mpi_comm );
    }

    // Wait for receives and verify:

    for ( size_t i = 0 ; i < recv_msg_count ; ++i ) {
      MPI_Status recv_status ;
      int recv_which = 0 ;
      int recv_size  = 0 ;

      MPI_Waitany( recv_msg_count , & recv_request[0] ,
                   & recv_which , & recv_status );

      const int recv_proc = recv_status.MPI_SOURCE ;

      MPI_Get_count( & recv_status , MPI_BYTE , & recv_size );

      // Verify message properly received:

      const int  expected_proc = data_map.host_recv(recv_which,0);
      const int  expected_size = data_map.host_recv(recv_which,1) *
                                 chunk_size * sizeof(ValueType);

      if ( ( expected_proc != recv_proc ) ||
           ( expected_size != recv_size ) ) {
        std::ostringstream msg ;
        msg << "AsyncExchange error:"
            << " P" << comm::rank( data_map.machine )
            << " received from P" << recv_proc
            << " size "     << recv_size
            << " expected " << expected_size
            << " from P"    << expected_proc ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Copy received data to device memory.

    Kokkos::deep_copy( dev_recv_buffer , host_recv_buffer );
  }
};

} // namespace Kokkos

#else /* ! #ifdef KOKKOS_ENABLE_MPI */

namespace Kokkos {

template< class ValueType , class Device >
class AsyncExchange< ValueType, Device , Kokkos::ParallelDataMap > {
public:

  typedef Device                                    execution_space ;
  typedef Kokkos::ParallelDataMap                   data_map_type ;
  typedef Kokkos::View< ValueType* , execution_space >  buffer_dev_type ;
  typedef typename buffer_dev_type::HostMirror      buffer_host_type ;

  buffer_dev_type      dev_buffer ;

public:

  const buffer_dev_type & buffer() const { return dev_buffer ; }

  AsyncExchange( const data_map_type & , const size_t )
  : dev_buffer()
  { }

  //------------------------------------------------------------------------

  void setup() { }

  void send_receive() { }
};

} // namespace Kokkos

#endif /* ! #ifdef KOKKOS_ENABLE_MPI */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_PARALLELDATAMAP_HPP */


