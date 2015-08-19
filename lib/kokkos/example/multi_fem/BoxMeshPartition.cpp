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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <BoxMeshPartition.hpp>

//----------------------------------------------------------------------------

namespace {

void box_partition( size_t ip , size_t up ,
                    const BoxType & box ,
                    BoxType * const p_box )
{
  const size_t np = up - ip ;

  if ( 1 == np ) {
    p_box[ip] = box ;
  }
  else {
    // Choose axis with largest count:

    const size_t n0 = box[0][1] - box[0][0] ;
    const size_t n1 = box[1][1] - box[1][0] ;
    const size_t n2 = box[2][1] - box[2][0] ;

    const size_t axis = n2 > n1 ? ( n2 > n0 ? 2 : ( n1 > n0 ? 1 : 0 ) ) :
                                  ( n1 > n0 ? 1 : 0 );

    const size_t n = box[ axis ][1] - box[ axis ][0] ;

    if ( 0 == np % 3 ) {
      const size_t np_part = np / 3 ; // exact

      const size_t nbox_low = (size_t)(( (double) n ) * ( 1.0 / 3.0 ));
      const size_t nbox_mid = (size_t)(( (double) n ) * ( 2.0 / 3.0 ));

      BoxType dbox_low = box ; // P = [ip,ip+np/3) 
      BoxType dbox_mid = box ; // P = [ip+np/3,ip+2*np/3) 
      BoxType dbox_upp = box ; // P = [ip+2*np/3,ip+np) 

      dbox_low[ axis ][1] = box[ axis ][0] + nbox_low ;
      dbox_mid[ axis ][1] = box[ axis ][0] + nbox_mid ;

      dbox_mid[ axis ][0] = dbox_low[ axis ][1];
      dbox_upp[ axis ][0] = dbox_mid[ axis ][1];

      box_partition( ip,           ip +   np_part, dbox_low , p_box );
      box_partition( ip+  np_part, ip + 2*np_part, dbox_mid , p_box );
      box_partition( ip+2*np_part, up,             dbox_upp , p_box );
    }
    else {
      const size_t np_low = np / 2 ; /* Rounded down */
      const size_t nbox_low = (size_t)
        (((double)n) * ( ((double) np_low ) / ((double) np ) ));

      BoxType dbox_low = box ;
      BoxType dbox_upp = box ;

      dbox_low[ axis ][1] = dbox_low[ axis ][0] + nbox_low ; 
      dbox_upp[ axis ][0] = dbox_low[ axis ][1];

      box_partition( ip, ip + np_low, dbox_low , p_box );
      box_partition( ip + np_low, up, dbox_upp , p_box );
    }
  }
}

size_t box_map_offset( const BoxType & local_use ,
                       const size_t global_i ,
                       const size_t global_j ,
                       const size_t global_k )

{
  const size_t max = std::numeric_limits<size_t>::max();

  const size_t n[3] =
    { local_use[0][1] - local_use[0][0] ,
      local_use[1][1] - local_use[1][0] ,
      local_use[2][1] - local_use[2][0] };

  const size_t use[3] = {
    ( global_i >= local_use[0][0] ? global_i - local_use[0][0] : max ) ,
    ( global_j >= local_use[1][0] ? global_j - local_use[1][0] : max ) ,
    ( global_k >= local_use[2][0] ? global_k - local_use[2][0] : max ) };

  const size_t offset =
    ( use[0] < n[0] && use[1] < n[1] && use[2] < n[2] ) ?
    ( use[0] + n[0] * ( use[1] + n[1] * use[2] ) ) : max ;

  if ( offset == max ) {
    std::ostringstream msg ;
    msg << "box_map_offset ERROR: "
        << " use " << local_use
        << " ( " << global_i
        << " , " << global_j
        << " , " << global_k
        << " )" ;
    throw std::runtime_error( msg.str() );
  }

  return offset ;
}

} // namespace

//----------------------------------------------------------------------------

void BoxBoundsLinear::apply(  const BoxType & box_global ,
                              const BoxType & box_part ,
                                    BoxType & box_interior ,
                                    BoxType & box_use ) const
{
  const unsigned ghost = 1 ;

  if ( 0 == count( box_part ) ) {
    box_interior = box_part ;
    box_use      = box_part ;
  }
  else {
    for ( size_t i = 0 ; i < 3 ; ++i ) {

      box_interior[i][0] =
        ( box_part[i][0] == box_global[i][0] )      ? box_part[i][0] : (
        ( box_part[i][0] + ghost < box_part[i][1] ) ? box_part[i][0] + ghost : 
                                                      box_part[i][1] );

      box_interior[i][1] =
        ( box_part[i][1] == box_global[i][1] )      ? box_part[i][1] : (
        ( box_part[i][0] + ghost < box_part[i][1] ) ? box_part[i][1] - ghost :
                                                      box_part[i][0] );

      box_use[i][0] = 
        ( box_part[i][0] > ghost + box_global[i][0] ) ? box_part[i][0] - ghost :
                                                        box_global[i][0] ;
      box_use[i][1] = 
        ( box_part[i][1] + ghost < box_global[i][1] ) ? box_part[i][1] + ghost :
                                                        box_global[i][1] ;
    }
  }
}

void BoxBoundsQuadratic::apply( const BoxType & box_global ,
                                const BoxType & box_part ,
                                      BoxType & box_interior ,
                                      BoxType & box_use ) const
{
  if ( 0 == count( box_part ) ) {
    box_interior = box_part ;
    box_use      = box_part ;
  }
  else {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      const bool odd = ( box_part[i][0] - box_global[i][0] ) & 01 ;

      const unsigned ghost = odd ? 1 : 2 ;

      box_interior[i][0] =
        ( box_part[i][0] == box_global[i][0] )      ? box_part[i][0] : (
        ( box_part[i][0] + ghost < box_part[i][1] ) ? box_part[i][0] + ghost : 
                                                      box_part[i][1] );

      box_interior[i][1] =
        ( box_part[i][1] == box_global[i][1] )      ? box_part[i][1] : (
        ( box_part[i][0] + ghost < box_part[i][1] ) ? box_part[i][1] - ghost :
                                                      box_part[i][0] );

      box_use[i][0] = 
        ( box_part[i][0] > ghost + box_global[i][0] ) ? box_part[i][0] - ghost :
                                                        box_global[i][0] ;
      box_use[i][1] = 
        ( box_part[i][1] + ghost < box_global[i][1] ) ? box_part[i][1] + ghost :
                                                        box_global[i][1] ;
    }
  }
}

//----------------------------------------------------------------------------

void box_partition_rcb( const BoxType        & root_box ,
                        std::vector<BoxType> & part_boxes )
{
  const BoxBoundsLinear use_boxes ;

  const size_t part_count = part_boxes.size();

  box_partition( 0 , part_count , root_box , & part_boxes[0] );

  // Verify partitioning

  size_t total_cell = 0 ;

  for ( size_t i = 0 ; i < part_count ; ++i ) {

    total_cell += count( part_boxes[i] );

    BoxType box_interior , box_use ;

    use_boxes.apply( root_box , part_boxes[i] , box_interior , box_use );

    if ( count( box_use ) < count( part_boxes[i] ) ||
         count( part_boxes[i] ) < count( box_interior ) ||
         part_boxes[i] != intersect( part_boxes[i] , box_use ) ||
         box_interior  != intersect( part_boxes[i] , box_interior )) {

      std::ostringstream msg ;

      msg << "box_partition_rcb ERROR : "
          << "part_boxes[" << i << "] = "
          << part_boxes[i]
          << " use " << box_use
          << " interior " << box_interior
          << std::endl 
          << "  part ^ use " << intersect( part_boxes[i] , box_use )
          << "  part ^ interior " << intersect( part_boxes[i] , box_interior );

      throw std::runtime_error( msg.str() );
    }

    for ( size_t j = i + 1 ; j < part_count ; ++j ) {
      const BoxType tmp = intersect( part_boxes[i] , part_boxes[j] );

      if ( count( tmp ) ) {
        throw std::runtime_error( std::string("box partition intersection") );
      }
    }
  }

  if ( total_cell != count( root_box ) ) {
    throw std::runtime_error( std::string("box partition count") );
  }
}

//----------------------------------------------------------------------------
         
size_t box_map_id( const BoxType & local_use ,
                   const std::vector<size_t> & local_use_id_map ,
                   const size_t global_i ,
                   const size_t global_j ,
                   const size_t global_k )

{
  const size_t offset =
    box_map_offset( local_use , global_i , global_j , global_k );
  return local_use_id_map[ offset ];
}
         
//----------------------------------------------------------------------------

void box_partition_maps( const BoxType              & root_box ,
                         const std::vector<BoxType> & part_boxes ,
                         const BoxBounds            & use_boxes ,
                         const size_t          my_part ,
                         BoxType             & my_use_box ,
                         std::vector<size_t> & my_use_id_map ,
                         size_t              & my_count_interior ,
                         size_t              & my_count_owned ,
                         size_t              & my_count_uses ,
                         std::vector<size_t> & my_part_counts ,
                         std::vector<std::vector<size_t> > & my_send_map )
{
  const size_t np = part_boxes.size();

  if ( np <= my_part ) {
    std::ostringstream msg ;
    msg << "box_partition_maps ERROR : "
        << " np(" << np << ") <= my_part(" << my_part << ")" ;
    throw std::runtime_error( msg.str() );
  }

  const BoxType my_owned_box = part_boxes[my_part];
  BoxType my_interior_box ;


  use_boxes.apply( root_box, my_owned_box, my_interior_box, my_use_box );

  my_count_interior = count( my_interior_box );
  my_count_owned    = count( my_owned_box );
  my_count_uses     = count( my_use_box );

  my_use_id_map.assign( my_count_uses , std::numeric_limits<size_t>::max() );

  // Order ids as { owned-interior , owned-parallel , received_{(p+i)%np} }

  size_t offset_interior = 0 ;
  size_t offset_parallel = my_count_interior ;

  for ( size_t iz = my_owned_box[2][0] ; iz < my_owned_box[2][1] ; ++iz ) {
  for ( size_t iy = my_owned_box[1][0] ; iy < my_owned_box[1][1] ; ++iy ) {
  for ( size_t ix = my_owned_box[0][0] ; ix < my_owned_box[0][1] ; ++ix ) {
    const size_t offset = box_map_offset( my_use_box , ix , iy , iz );
    if ( contain( my_interior_box , ix , iy , iz ) ) {
      my_use_id_map[ offset ] = offset_interior++ ;
    }
    else {
      my_use_id_map[ offset ] = offset_parallel++ ;
    }
  }}}


  my_part_counts.assign( np , (size_t) 0 );
  my_send_map.assign( np , std::vector<size_t>() );

  my_part_counts[0] = my_count_owned ;

  for ( size_t i = 1 ; i < np ; ++i ) {

    const size_t ip = ( my_part + i ) % np ;

    const BoxType p_owned_box = part_boxes[ip];
    BoxType p_use_box , p_interior_box ;
    use_boxes.apply( root_box, p_owned_box, p_interior_box, p_use_box );

    const BoxType recv_box = intersect( my_use_box , p_owned_box );
    const BoxType send_box = intersect( my_owned_box , p_use_box );

    if ( 0 != ( my_part_counts[i] = count( recv_box ) ) ) {
      for ( size_t iz = recv_box[2][0] ; iz < recv_box[2][1] ; ++iz ) {
      for ( size_t iy = recv_box[1][0] ; iy < recv_box[1][1] ; ++iy ) {
      for ( size_t ix = recv_box[0][0] ; ix < recv_box[0][1] ; ++ix ) {
        const size_t offset = box_map_offset( my_use_box , ix , iy , iz );
        my_use_id_map[ offset ] = offset_parallel++ ;
      }}}
    }

    if ( 0 != count( send_box ) ) {
      for ( size_t iz = send_box[2][0] ; iz < send_box[2][1] ; ++iz ) {
      for ( size_t iy = send_box[1][0] ; iy < send_box[1][1] ; ++iy ) {
      for ( size_t ix = send_box[0][0] ; ix < send_box[0][1] ; ++ix ) {
        const size_t offset = box_map_offset( my_use_box , ix , iy , iz );

        my_send_map[ i ].push_back( my_use_id_map[ offset ] );
      }}}
    }
  }
}


