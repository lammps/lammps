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

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshPartition.hpp>

//----------------------------------------------------------------------------

void test_box_partition( bool print )
{
  const size_t np_max = 10000 ;

  const BoxBoundsLinear use_box ;

  BoxType root_box ;

  root_box[0][0] = 0 ; root_box[0][1] = 100 ;
  root_box[1][0] = 0 ; root_box[1][1] = 200 ;
  root_box[2][0] = 0 ; root_box[2][1] = 300 ;

  const size_t cell_total =
    ( root_box[0][1] - root_box[0][0] ) *
    ( root_box[1][1] - root_box[1][0] ) *
    ( root_box[2][1] - root_box[2][0] );

  for ( size_t np = 2 ; np < np_max ; np = 2 * ( np + 1 ) ) {

    std::vector<BoxType> part_boxes( np );

    box_partition_rcb( root_box , part_boxes );

    size_t cell_goal = ( cell_total + np - 1 ) / np ;
    size_t cell_max = 0 ;

    for ( size_t i = 0 ; i < np ; ++i ) {
      cell_max = std::max( cell_max , count( part_boxes[i] ) );
    }

    if ( print ) {
      std::cout << std::endl
                << "box_part( " << np 
                << " ) max( " << cell_max
                << " ) goal( " << cell_goal
                << " ) ratio( " << double(cell_max) / double(cell_goal)
                << " )" << std::endl ;
    }

    const size_t nsample = std::min(np,(size_t)4);
    const size_t stride = ( np + nsample - 1 ) / nsample ;

    for ( size_t my_part = 0 ; my_part < np ; my_part += stride ) {
      BoxType             my_use_box ;
      std::vector<size_t> my_use_id_map ;
      size_t              my_count_interior ;
      size_t              my_count_owned ;
      size_t              my_count_uses ;
      std::vector<size_t> my_recv_counts ;
      std::vector<std::vector<size_t> > my_send_map ;

      size_t count_verify = 0 ;

      box_partition_maps( root_box , part_boxes ,
                          use_box , my_part ,
                          my_use_box , my_use_id_map ,
                          my_count_interior ,
                          my_count_owned ,
                          my_count_uses ,
                          my_recv_counts ,
                          my_send_map );

      count_verify = my_count_owned ;

      if ( print ) {
        std::cout << "  my_part(" << my_part << ") layout { "
                  << "P" << my_part
                  << "(" << my_count_interior
                  << "," << ( my_count_owned - my_count_interior )
                  << ")" ;
      }

      for ( size_t i = 1 ; i < np ; ++i ) {
        if ( my_recv_counts[i] ) {
          count_verify += my_recv_counts[i] ;
          const size_t ip = ( my_part + i ) % np ;

          if ( print ) {
            std::cout << " P" << ip << "(" << my_recv_counts[i] << ")" ;
          }

          // Compare recv & send lists

          BoxType             ip_use_box ;
          std::vector<size_t> ip_use_id_map ;
          size_t              ip_count_interior ;
          size_t              ip_count_owned ;
          size_t              ip_count_uses ;
          std::vector<size_t> ip_recv_counts ;
          std::vector<std::vector<size_t> > ip_send_map ;

          box_partition_maps( root_box , part_boxes ,
                              use_box , ip ,
                              ip_use_box , ip_use_id_map ,
                              ip_count_interior ,
                              ip_count_owned ,
                              ip_count_uses ,
                              ip_recv_counts ,
                              ip_send_map );

          // Sent by ip, received by my_part:

          const BoxType recv_send = intersect( part_boxes[ip] , my_use_box );
          const size_t recv_send_count = count( recv_send );

          const size_t j = ( my_part + np - ip ) % np ;

          if ( recv_send_count != my_recv_counts[i] ||
               recv_send_count != ip_send_map[j].size() ) {
            throw std::runtime_error( std::string("bad recv/send map") );
          }
        }
      }
      if ( print ) { std::cout << " }" << std::endl ; }

      if ( count_verify != my_count_uses ) {
        throw std::runtime_error( std::string("bad partition map") );
      }
    }
  }
}


