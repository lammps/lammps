/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

#include <stddef.h>
#include <limits>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemoryTracking.hpp>

namespace Kokkos {
namespace Impl {
namespace {

//----------------------------------------------------------------------------
// Fast search for result[-1] <= val < result[0].
// Requires result[max] == upper_bound.
// Start with a binary search until the search range is
// less than LINEAR_LIMIT, then switch to linear search.

int upper_bound( const ptrdiff_t * const begin , unsigned length ,
                 const ptrdiff_t val )
{
  enum { LINEAR_LIMIT = 32 };

  // precondition: begin[length-1] == std::numeric_limits<ptrdiff_t>::max()

  const ptrdiff_t * first = begin ;

  while ( LINEAR_LIMIT < length ) {
    unsigned          half   = length >> 1 ;
    const ptrdiff_t * middle = first + half ;

    if ( val < *middle ) {
      length = half ;
    }
    else {
      first   = ++middle ;
      length -= ++half ;
    }
  }

  for ( ; ! ( val < *first ) ; ++first ) {}

  return first - begin ;
}

} // namespace

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

MemoryTracking::MemoryTracking( const std::string & space )
  : m_space( space ), m_tracking(), m_tracking_end()
{
  ptrdiff_t max = std::numeric_limits<ptrdiff_t>::max();
  void * const ptr = reinterpret_cast<void*>( max );

  m_tracking.reserve(64);
  m_tracking_end.reserve(64);

  // Sentinal value of end

  m_tracking.push_back( new MemoryTrackingEntry( "sentinal" , typeid(void) , ptr , 0 ) );
  m_tracking_end.push_back( max );
}

MemoryTracking::~MemoryTracking()
{
  const ptrdiff_t max = std::numeric_limits<ptrdiff_t>::max();

  try {
    if ( 1 < m_tracking.size() ) {
      std::cerr << m_space << " destroyed with memory leaks:" << std::endl ;
      print( std::cerr , std::string("  ") );
    }
    else if ( 1 != m_tracking_end.size() || m_tracking_end.back() != max ) {
      std::cerr << m_space << " corrupted data structure" << std::endl ;
    }

    // Deallocate memory within the try-catch block:
    m_space        = std::string();
    m_tracking     = std::vector<MemoryTrackingEntry*>();
    m_tracking_end = std::vector<ptrdiff_t>();

  } catch( ... ) {}
}

void MemoryTracking::insert( MemoryTrackingEntry * entry )
{
  const ptrdiff_t max = std::numeric_limits<ptrdiff_t>::max();

  const bool ok_exists = ! m_tracking_end.empty();

  const bool ok_range = entry &&
                        0 < entry->begin &&
                            entry->begin < entry->end &&
                                           entry->end < max ;

  int i = -1 ;

  if ( ok_exists && ok_range ) {

    i = upper_bound( & m_tracking_end[0] , m_tracking_end.size() , entry->begin );

    // Guaranteed:
    //   a) entry->begin < m_tracking_end[i]
    //   b) i == 0 || m_tracking_end[i-1] <= entry->begin

    if ( entry->end <= m_tracking[i]->begin ) {

      // Non-overlapping range:
      // m_tracking[i-1].end <= entry->begin < entry->end <= m_tracking[i].begin

      entry->m_count = 1 ;

      m_tracking.insert(     m_tracking.begin() + i , entry );
      m_tracking_end.insert( m_tracking_end.begin() + i , entry->end );
    }
  }

  if ( ! ok_exists || ! ok_range || -1 == i ) {
    std::ostringstream msg ;
    msg << "MemoryTracking(" << m_space << ")::insert( " ;
    entry->print( msg );
    msg << " ) ERROR: " ;

    if ( ! ok_range ) {
      msg << "Invalid memory range" ;
    }
    else {
      msg << "Overlapping memory range with " ;
      m_tracking[i]->print( msg );
    }
    msg << " )" ;
    throw_runtime_exception( msg.str() );
  }
}

void MemoryTracking::increment( const void * ptr )
{
  if ( ptr ) {
    const ptrdiff_t p = reinterpret_cast<ptrdiff_t>( ptr );

    bool error = m_tracking_end.empty();

    if ( ! error ) {

      const int i = upper_bound( & m_tracking_end[0] , m_tracking_end.size() , p );

      error = p < m_tracking[i]->begin ;

      if ( ! error ) {
        ++( m_tracking[i]->m_count );
      }
    }

    if ( error ) {
      std::ostringstream msg ;
      msg << "MemoryTracking(" << m_space
          << ")::increment( " << p << " ) ERROR: Not being tracked" ;
      throw_runtime_exception( msg.str() );
    }
  }
}

void MemoryTracking::decrement( const void * ptr )
{
  if ( ptr ) {
    const ptrdiff_t p = reinterpret_cast<ptrdiff_t>( ptr );

    bool error = m_tracking_end.empty();

    if ( ! error ) {

      const int i = upper_bound( & m_tracking_end[0] , m_tracking_end.size() , p );

      error = p < m_tracking[i]->begin ;

      if ( ! error && ( 0 == --( m_tracking[i]->m_count ) ) ) {
        delete m_tracking[i] ;

        m_tracking.erase(     m_tracking.begin() + i );
        m_tracking_end.erase( m_tracking_end.begin() + i );
      }
    }

    if ( error ) {
      std::ostringstream msg ;
      msg << "MemoryTracking(" << m_space
          << ")::decrement( " << p << " ) ERROR: Not being tracked" 
          << std::endl ;
      std::cerr << msg.str();
    }
  }
}

MemoryTrackingEntry *
MemoryTracking::query( const void * ptr ) const
{
  MemoryTrackingEntry * result = 0 ;

  if ( ptr && ! m_tracking_end.empty() ) {
    const ptrdiff_t p = reinterpret_cast<ptrdiff_t>( ptr );

    const int i = upper_bound( & m_tracking_end[0] , m_tracking_end.size() , p );

    if ( m_tracking[i]->begin <= p ) result = m_tracking[i] ;
  }

  return result ;
}

void MemoryTracking::print( std::ostream & s , const std::string & lead ) const
{
  // Don't print the sentinal value:
  const size_t n = m_tracking.empty() ? 0 : m_tracking.size() - 1 ;

  for ( size_t i = 0 ; i < n ; ++i ) {
    s << lead ;
    m_tracking[i]->print( s );
    s << std::endl ;
  }
}

MemoryTrackingEntry::~MemoryTrackingEntry()
{}

void MemoryTrackingEntry::print( std::ostream & s ) const
{
  s << "{ "
    << "label("  << label << ") "
    << "typeid(" << type.name() << ") "
    << "range[ " << ((void*)begin) << " : " << ((void*)end) << " ) "
    << "count("  << m_count << ") }" ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */


