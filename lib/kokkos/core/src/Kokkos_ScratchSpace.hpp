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

#ifndef KOKKOS_SCRATCHSPACE_HPP
#define KOKKOS_SCRATCHSPACE_HPP

#include <stdio.h>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Scratch memory space associated with an execution space.
 *
 */
template< class ExecSpace >
class ScratchMemorySpace {
public:

  // Alignment of memory chunks returned by 'get'
  // must be a power of two
  enum { ALIGN = 8 };

private:

  mutable char * m_iter ;
  char *         m_end ;

  ScratchMemorySpace();
  ScratchMemorySpace & operator = ( const ScratchMemorySpace & );

  enum { MASK = ALIGN - 1 }; // Alignment used by View::shmem_size

public:

  //! Tag this class as a memory space
  typedef ScratchMemorySpace                memory_space ;
  typedef ExecSpace                         execution_space ;
  typedef typename ExecSpace::array_layout  array_layout ;
  typedef typename ExecSpace::size_type     size_type ;

  template< typename IntType >
  KOKKOS_INLINE_FUNCTION static
  IntType align( const IntType & size )
    { return ( size + MASK ) & ~MASK ; }

  template< typename IntType >
  KOKKOS_INLINE_FUNCTION
  void* get_shmem (const IntType& size) const {
    void* tmp = m_iter ;
    if (m_end < (m_iter += align (size))) {
      m_iter -= align (size); // put it back like it was
      printf ("ScratchMemorySpace<...>::get_shmem: Failed to allocate %ld byte(s); remaining capacity is %ld byte(s)\n", long(size), long(m_end-m_iter));
      tmp = 0;
    }
    return tmp;
  }

  template< typename IntType >
  KOKKOS_INLINE_FUNCTION
  ScratchMemorySpace( void * ptr , const IntType & size )
    : m_iter( (char *) ptr )
    , m_end(  m_iter + size )
    {}
};

} // namespace Kokkos

#endif /* #ifndef KOKKOS_SCRATCHSPACE_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

