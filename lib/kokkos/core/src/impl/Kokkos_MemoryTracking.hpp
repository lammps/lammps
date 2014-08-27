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
#include <utility>
#include <vector>
#include <string>
#include <typeinfo>
#include <iosfwd>

namespace Kokkos {
namespace Impl {

class MemoryTracking ;

class MemoryTrackingEntry {
public:
  const std::string      label ;
  const std::type_info & type ;
  const ptrdiff_t        begin ;
  const ptrdiff_t        end ;
private:
  unsigned m_count ;
protected:

  MemoryTrackingEntry( const std::string    & arg_label ,
                       const std::type_info & arg_type ,
                       const void * const     arg_begin ,
                       const ptrdiff_t        arg_bytes )
    : label( arg_label )
    , type(  arg_type )
    , begin( reinterpret_cast<ptrdiff_t>( arg_begin ) )
    , end(   reinterpret_cast<ptrdiff_t>(
               reinterpret_cast<const unsigned char *>( arg_begin ) + arg_bytes ) )
    , m_count( 0 )
    {}

public:

  unsigned count() const { return m_count ; }

  virtual void print( std::ostream & ) const ;

  virtual ~MemoryTrackingEntry();

private:

  MemoryTrackingEntry();
  MemoryTrackingEntry( const MemoryTrackingEntry & rhs );
  MemoryTrackingEntry & operator = ( const MemoryTrackingEntry & rhs );

  friend class MemoryTracking ;
};


class MemoryTracking {
public:

  /** \brief  Track a memory range defined by the entry.
   *          This entry must be allocated via 'new'.
   */
  void insert( MemoryTrackingEntry * entry );

  /** \brief  Decrement the tracked memory range.
   *          If the count is zero then the entry is deleted
   *          via the 'delete' operator.
   */
  void decrement( const void * ptr );

  /** \brief  Increment the tracking count.  */
  void increment( const void * ptr );

  /** \brief  Query a tracked memory range. */
  MemoryTrackingEntry * query( const void * ptr ) const ;

  /** \brief  Call the 'print' method on all entries. */
  void print( std::ostream & , const std::string & lead ) const ;

  size_t size() const { return m_tracking.size(); }

  template< typename iType >
  MemoryTracking & operator[]( const iType & i ) const
    { return *m_tracking[i]; }

  /** \brief Construct with a name for error messages */
  explicit MemoryTracking( const std::string & space );

  /** \brief  Print memory leak warning for all entries. */
  ~MemoryTracking();

  /** \brief Query if constructed */
  bool exists() const { return ! m_tracking_end.empty(); }

private:
  MemoryTracking();
  MemoryTracking( const MemoryTracking & );
  MemoryTracking & operator = ( const MemoryTracking & );

  std::string                        m_space ;
  std::vector<MemoryTrackingEntry*>  m_tracking ;
  std::vector<ptrdiff_t>             m_tracking_end ;
};

} /* namespace Impl */
} /* namespace Kokkos */

#endif

