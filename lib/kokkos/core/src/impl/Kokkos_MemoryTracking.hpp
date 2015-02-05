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

#ifndef KOKKOS_MEMORY_TRACKING_HPP
#define KOKKOS_MEMORY_TRACKING_HPP

#include <cstddef>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Impl {
namespace {

// Fast search for result[-1] <= val < result[0].
// Requires result[max] == upper_bound.
// Start with a binary search until the search range is
// less than LINEAR_LIMIT, then switch to linear search.

int memory_tracking_upper_bound( const ptrdiff_t * const begin
                               , unsigned length
                               , const ptrdiff_t value )
{
  enum { LINEAR_LIMIT = 32 };

  // precondition: begin[length-1] == std::numeric_limits<ptrdiff_t>::max()

  const ptrdiff_t * first = begin ;

  while ( LINEAR_LIMIT < length ) {
    unsigned          half   = length >> 1 ;
    const ptrdiff_t * middle = first + half ;

    if ( value < *middle ) {
      length = half ;
    }
    else {
      first   = ++middle ;
      length -= ++half ;
    }
  }

  for ( ; ! ( value < *first ) ; ++first ) {}

  return first - begin ;
}

template< class AttributeType = size_t >
class MemoryTracking {
public:

  class Entry {
  private:

    friend class MemoryTracking ;

    enum { LABEL_LENGTH = 128 };

    Entry( const Entry & );
    Entry & operator = ( const Entry & );

    ~Entry() {}

    Entry()
     : m_count(0)
     , m_alloc_ptr( reinterpret_cast<void*>( std::numeric_limits<ptrdiff_t>::max() ) )
     , m_alloc_size(0)
     , m_attribute()
     { strcpy( m_label , "sentinel" ); }

    Entry( const std::string & arg_label 
         , void * const        arg_alloc_ptr 
         , size_t const        arg_alloc_size )
      : m_count( 0 )
      , m_alloc_ptr( arg_alloc_ptr )
      , m_alloc_size( arg_alloc_size )
      , m_attribute()
      {
        strncpy( m_label , arg_label.c_str() , LABEL_LENGTH );
        m_label[ LABEL_LENGTH - 1 ] = 0 ;
      }

    char    m_label[ LABEL_LENGTH ] ;
    size_t  m_count ;

  public:

    void * const   m_alloc_ptr ;
    size_t const   m_alloc_size ;
    AttributeType  m_attribute ;

    size_t       count() const { return m_count ; }
    const char * label() const { return m_label ; }

    void print( std::ostream & oss ) const
      {
        oss << "{ \"" << m_label
            << "\" count(" << m_count
            << ") memory[ " << m_alloc_ptr
            << " + " << m_alloc_size
            << " ]" ;
      }
  };

  //------------------------------------------------------------
  /** \brief  Track a memory range defined by the entry.
   *          Return the input entry pointer for success.
   *          Throw exception for failure.
   */
  Entry * insert( const std::string & arg_label
                , void * const  arg_alloc_ptr
                , size_t const  arg_alloc_size
                )
    {
      Entry * result = 0 ;

      const ptrdiff_t alloc_begin = reinterpret_cast<ptrdiff_t>(arg_alloc_ptr);
      const ptrdiff_t alloc_end   = alloc_begin + arg_alloc_size ;

      const bool ok_exist = ! m_tracking_end.empty(); 

      const bool ok_input =
        ok_exist &&
        ( 0 < alloc_begin ) &&
            ( alloc_begin < alloc_end ) &&
                          ( alloc_end < std::numeric_limits<ptrdiff_t>::max() ); 

      const int i = ok_input
                  ? memory_tracking_upper_bound( & m_tracking_end[0] , m_tracking_end.size() , alloc_end )
                  : -1 ;

      const bool ok_range = ( 0 <= i ) && ( alloc_end <= reinterpret_cast<ptrdiff_t>( m_tracking[i]->m_alloc_ptr ) );

      // allocate the new entry only if the vector inserts succeed.
      const bool ok_insert =
        ok_range &&
        ( alloc_end == *m_tracking_end.insert(m_tracking_end.begin()+i,alloc_end) ) &&
        ( 0 == *m_tracking.insert(m_tracking.begin()+i,0) ) &&
        ( 0 != ( result = new Entry(arg_label,arg_alloc_ptr,arg_alloc_size) ) );

      if ( ok_insert ) {
        result->m_count = 1 ;
        m_tracking[i] = result ;
      }
      else {
        std::ostringstream msg ;
        msg << m_space
            << "::insert( " << arg_label
            << " , " << arg_alloc_ptr
            << " , " << arg_alloc_size
            << " ) ERROR : " ;
        if ( ! ok_exist ) {
          msg << " called after return from main()" ;
        }
        else if ( ! ok_input ) {
          msg << " bad allocation range" ;
        }
        else if ( ! ok_range ) {
          msg << " overlapping memory range with"
              << " { " << m_tracking[i]->m_label
              << " , " << m_tracking[i]->m_alloc_ptr
              << " , " << m_tracking[i]->m_alloc_size
              << " }" ;
        }
        else {
          msg << " internal allocation error" ;
        }
        Kokkos::Impl::throw_runtime_exception( msg.str() );
      }

      return result ;
    }

  /** \brief  Decrement the tracked memory range.
   *          If the count is zero then return the originally inserted pointer.
   *          If the count is non zero then return zero.
   */
  void * decrement( void const * const ptr )
    {
      void * result = 0 ;

      if ( ptr ) {
        const bool ok_exist = ! m_tracking_end.empty();

        const int i = ok_exist
                    ? memory_tracking_upper_bound( & m_tracking_end[0] , m_tracking_end.size() , reinterpret_cast<ptrdiff_t>(ptr) )
                    : -1 ;

        const bool ok_found = ( 0 <= i ) && ( reinterpret_cast<ptrdiff_t>( m_tracking[i]->m_alloc_ptr ) <=
                                              reinterpret_cast<ptrdiff_t>(ptr) );

        if ( ok_found ) {
          if ( 0 == --( m_tracking[i]->m_count ) ) {
            result = m_tracking[i]->m_alloc_ptr ;          
            delete m_tracking[i] ;
            m_tracking.erase(     m_tracking.begin() + i );
            m_tracking_end.erase( m_tracking_end.begin() + i );
          }
        }
        else {
          // Don't throw as this is likely called from within a destructor.
          std::cerr << m_space
                    << "::decrement( " << ptr << " ) ERROR : " 
                    << ( ! ok_exist ? " called after return from main()" 
                                    : " memory not being tracked" )
                    << std::endl ;
          std::cerr.flush();
        }
      }
      return result ;
    }

  /** \brief  Increment the tracking count.  */
  void increment( void const * const ptr )
    {
      if ( ptr ) {
        const bool ok_exist = ! m_tracking_end.empty();

        const int i = ok_exist
                    ? memory_tracking_upper_bound( & m_tracking_end[0] , m_tracking_end.size() , reinterpret_cast<ptrdiff_t>(ptr) )
                    : -1 ;

        const bool ok_found = ( 0 <= i ) && ( reinterpret_cast<ptrdiff_t>( m_tracking[i]->m_alloc_ptr ) <=
                                              reinterpret_cast<ptrdiff_t>(ptr) );

        if ( ok_found ) {
          ++( m_tracking[i]->m_count );
        }
        else {
          std::ostringstream msg ;
          msg << m_space
              << "::increment( " << ptr << " ) ERROR : "
              << ( ! ok_exist ? " called after return from main()" 
                              : " memory not being tracked" )
              << std::endl ;
          Kokkos::Impl::throw_runtime_exception( msg.str() );
        }
      }
    }

  /** \brief  Query a tracked memory range.
   *          Return zero for not found.
   */
  Entry * query( void const * const ptr ) const
    {
      const bool ok_exist = ! m_tracking_end.empty();

      const int i = ( ok_exist && ptr )
                  ? memory_tracking_upper_bound( & m_tracking_end[0] , m_tracking_end.size() , reinterpret_cast<ptrdiff_t>(ptr) )
                  : -1 ;

      const bool ok_found = ( 0 <= i ) && ( reinterpret_cast<ptrdiff_t>( m_tracking[i]->m_alloc_ptr ) <=
                                            reinterpret_cast<ptrdiff_t>(ptr) );

      return ok_found ? m_tracking[i] : (Entry *) 0 ;
    }

  /** \brief  Call the 'print' method on all entries. */
  void print( std::ostream & oss , const std::string & lead ) const
    {
      const size_t n = m_tracking.empty() ? 0 : m_tracking.size() - 1 ;
      for ( size_t i = 0 ; i < n ; ++i ) {
        oss << lead ;
        m_tracking[i]->print( oss );
        oss << std::endl ;
      }
    }

  size_t size() const { return m_tracking.size(); }

  template< typename iType >
  MemoryTracking & operator[]( const iType & i ) const
    { return *m_tracking[i]; }

  /** \brief Construct with a name for error messages */
  explicit MemoryTracking( const std::string & space_name )
    : m_space( space_name )
    , m_tracking()
    , m_tracking_end()
    , m_sentinel()
    {
      m_tracking.reserve( 512 );
      m_tracking_end.reserve( 512 );
      m_tracking.push_back( & m_sentinel );
      m_tracking_end.push_back( reinterpret_cast<ptrdiff_t>( m_sentinel.m_alloc_ptr ) );
    }

  /** \brief  Print memory leak warning for all entries. */
  ~MemoryTracking()
    {
      try {
        const ptrdiff_t max = std::numeric_limits<ptrdiff_t>::max();

        if ( 1 < m_tracking.size() ) {
          std::cerr << m_space << " destroyed with memory leaks:" ;
          print( std::cerr , std::string("  ") );
        }
        else if ( m_tracking.empty() || max != m_tracking_end.back() ) {
          std::cerr << m_space << " corrupted data structure" << std::endl ;
        }

        m_space = std::string();
        m_tracking = std::vector<Entry*>();
        m_tracking_end = std::vector<ptrdiff_t>();
      }
      catch( ... ) {}
    }

  const std::string & label() const { return m_space ; }

private:
  MemoryTracking();
  MemoryTracking( const MemoryTracking & );
  MemoryTracking & operator = ( const MemoryTracking & );

  std::string             m_space ;
  std::vector<Entry*>     m_tracking ;
  std::vector<ptrdiff_t>  m_tracking_end ;
  Entry                   m_sentinel ;
};

} /* namespace */
} /* namespace Impl */
} /* namespace Kokkos */

#endif

