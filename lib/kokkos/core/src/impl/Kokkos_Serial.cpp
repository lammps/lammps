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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <sstream>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
namespace {

struct Sentinel {

  void *   m_scratch ;
  unsigned m_reduce_end ;
  unsigned m_shared_end ;

  Sentinel() : m_scratch(0), m_reduce_end(0), m_shared_end(0) {}

  ~Sentinel()
    {
      if ( m_scratch ) { free( m_scratch ); }
      m_scratch = 0 ;
      m_reduce_end = 0 ;
      m_shared_end = 0 ;
    }

  static Sentinel & singleton();
};

Sentinel & Sentinel::singleton()
{
  static Sentinel s ; return s ;
}

inline
unsigned align( unsigned n )
{
  enum { ALIGN = 0x0100 /* 256 */ , MASK = ALIGN - 1 };
  return ( n + MASK ) & ~MASK ;
}

} // namespace

SerialTeamMember::SerialTeamMember( int arg_league_rank
                                  , int arg_league_size
                                  , int arg_shared_size
                                  )
  : m_space( ((char *) Sentinel::singleton().m_scratch) + Sentinel::singleton().m_reduce_end
           , arg_shared_size )
  , m_league_rank( arg_league_rank )
  , m_league_size( arg_league_size )
{}

} // namespace Impl

void * Serial::scratch_memory_resize( unsigned reduce_size , unsigned shared_size )
{
  static Impl::Sentinel & s = Impl::Sentinel::singleton();

  reduce_size = Impl::align( reduce_size );
  shared_size = Impl::align( shared_size );

  if ( ( s.m_reduce_end < reduce_size ) ||
       ( s.m_shared_end < s.m_reduce_end + shared_size ) ) {

    if ( s.m_scratch ) { free( s.m_scratch ); }
  
    if ( s.m_reduce_end < reduce_size ) s.m_reduce_end = reduce_size ;
    if ( s.m_shared_end < s.m_reduce_end + shared_size ) s.m_shared_end = s.m_reduce_end + shared_size ;

    s.m_scratch = malloc( s.m_shared_end );
  }

  return s.m_scratch ;
}

} // namespace Kokkos

