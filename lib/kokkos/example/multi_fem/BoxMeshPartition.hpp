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

#ifndef BOXMESHPARTITION_HPP
#define BOXMESHPARTITION_HPP

#include <cstddef>
#include <utility>
#include <vector>
#include <iostream>

//----------------------------------------------------------------------------

struct BoxType {
  size_t data[3][2] ;

  typedef size_t range_type[2] ;

  inline
  const range_type & operator[]( size_t i ) const { return data[i]; }

  inline
  range_type & operator[]( size_t i ) { return data[i]; }

  inline
  bool operator == ( const BoxType & rhs ) const
  {
    return data[0][0] == rhs.data[0][0] && data[0][1] == rhs.data[0][1] &&
           data[1][0] == rhs.data[1][0] && data[1][1] == rhs.data[2][1] &&
           data[2][0] == rhs.data[2][0] && data[2][1] == rhs.data[2][1] ;
  }

  inline
  bool operator != ( const BoxType & rhs ) const
  {
    return data[0][0] != rhs.data[0][0] || data[0][1] != rhs.data[0][1] ||
           data[1][0] != rhs.data[1][0] || data[1][1] != rhs.data[1][1] ||
           data[2][0] != rhs.data[2][0] || data[2][1] != rhs.data[2][1] ;
  }
};

inline
size_t count( const BoxType & b )
{
  size_t n = 1 ;
  for ( size_t i = 0 ; i < 3 ; ++i ) {
    n *= b[i][1] > b[i][0] ? b[i][1] - b[i][0] : 0 ;
  }
  return n ;
}

inline
bool contain( const BoxType & b , size_t i , size_t j , size_t k )
{
  return b[0][0] <= i && i < b[0][1] &&
         b[1][0] <= j && j < b[1][1] &&
         b[2][0] <= k && k < b[2][1] ;
}

inline
BoxType intersect( const BoxType & x , const BoxType & y )
{
  BoxType z ;
  for ( size_t i = 0 ; i < 3 ; ++i ) {
    z[i][0] = std::max( x[i][0] , y[i][0] );    
    z[i][1] = std::min( x[i][1] , y[i][1] );    
  }

  return z ;
}

inline
std::ostream & operator << ( std::ostream & s , const BoxType & box )
{
  s << "{ "
    << box[0][0] << " " << box[0][1] << " , "
    << box[1][0] << " " << box[1][1] << " , "
    << box[2][0] << " " << box[2][1] << " }" ;
  return s ;
}

//----------------------------------------------------------------------------

class BoxBounds {
public:
  /** \brief  Default bounds to one layer of ghosting */
  virtual
  void apply( const BoxType & box_global ,
              const BoxType & box_part ,
                    BoxType & box_interior ,
                    BoxType & box_use ) const = 0 ;

  virtual ~BoxBounds() {}
  BoxBounds() {}
};

class BoxBoundsLinear : public BoxBounds
{
public:
  /** \brief  Default bounds to one layer of ghosting */
  virtual
  void apply( const BoxType & box_global ,
              const BoxType & box_part ,
                    BoxType & box_interior ,
                    BoxType & box_use ) const ;

  virtual ~BoxBoundsLinear() {}
  BoxBoundsLinear() {}
};

class BoxBoundsQuadratic : public BoxBounds {
public:
  /** \brief  Quadratic mesh: even ordinates have two layers,
   *          odd ordinates have one layer.
   */
  virtual
  void apply( const BoxType & box_global ,
              const BoxType & box_part ,
                    BoxType & box_interior ,
                    BoxType & box_use ) const ;

  virtual ~BoxBoundsQuadratic() {}
  BoxBoundsQuadratic() {}
};

//----------------------------------------------------------------------------
/* Partition box into part_boxes.size() sub-boxes */

void box_partition_rcb( const BoxType        & root_box ,
                        std::vector<BoxType> & part_boxes );

//----------------------------------------------------------------------------
/* Determine local id layout and communication maps for partitioned boxes.
 *
 *  Local ids are laid out as follows:
 *    { [ owned-interior ids not sent ] ,
 *      [ owned-boundary ids to be sent to other processes ] ,
 *      [ received ids from processor ( my_part + 1 ) % part_count ]
 *      [ received ids from processor ( my_part + 2 ) % part_count ]
 *      [ received ids from processor ( my_part + 3 ) % part_count ]
 *      ... };
 *
 *  This layout allows
 *  (1) received data to be copied into a contiguous block of memory
 *  (2) send data to be extracted from a contiguous block of memory.
 */
void box_partition_maps(
  const BoxType              & root_box ,   // [in] Global box
  const std::vector<BoxType> & part_boxes , // [in] Partitioned boxes
  const BoxBounds            & use_boxes ,  // [in] Ghost boundaries
  const size_t          my_part ,           // [in] My local part
  BoxType             & my_use_box ,        // [out] My used box with ghost
  std::vector<size_t> & my_use_id_map ,     // [out] Local ordering map
  size_t              & my_count_interior , // [out] How many interior
  size_t              & my_count_owned ,    // [out] How many owned
  size_t              & my_count_uses ,     // [out] How may used
  std::vector<size_t> & my_part_counts ,    // [out] Partitioning of my_use_id_map
  std::vector<std::vector<size_t> > & my_send_map ); // [out] Send id map

/*  Mapping of cartesian coordinate to local id */
size_t box_map_id( const BoxType             & my_use_box ,
                   const std::vector<size_t> & my_use_id_map ,
                   const size_t global_i ,
                   const size_t global_j ,
                   const size_t global_k );

//----------------------------------------------------------------------------

#endif /* #ifndef BOXMESHPARTITION_HPP */

