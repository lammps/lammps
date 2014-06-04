/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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

#ifndef KOKKOS_MEMORYTRAITS_HPP
#define KOKKOS_MEMORYTRAITS_HPP

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Memory access traits for views, an extension point.
 *
 *  These traits should be orthogonal.  If there are dependencies then
 *  the MemoryTraits template must detect and enforce dependencies.
 *
 *  A zero value is the default for a View, indicating that none of
 *  these traits are present.
 */
enum MemoryTraitsFlags
  { Unmanaged  = 0x01
  , RandomAccess = 0x02
  };

template < unsigned T >
struct MemoryTraits {
  enum { Unmanaged  = T & unsigned(Kokkos::Unmanaged) };
  enum { RandomAccess = T & unsigned(Kokkos::RandomAccess) };

  typedef MemoryTraits memory_traits ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

typedef Kokkos::MemoryTraits<0> MemoryManaged ;
typedef Kokkos::MemoryTraits< Kokkos::Unmanaged > MemoryUnmanaged ;
typedef Kokkos::MemoryTraits< Kokkos::Unmanaged | Kokkos::RandomAccess > MemoryRandomAccess ;

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief Memory alignment settings
 *
 *  Sets global value for memory alignment.
 *  Enable compatibility of views from different devices with static stride.
 *  Use compiler flag to enable overwrites.
 */
enum { MEMORY_ALIGNMENT =
#if defined( KOKKOS_MEMORY_ALIGNMENT )
  KOKKOS_MEMORY_ALIGNMENT
#else
  128
#endif
  };

enum { MEMORY_ALIGNMENT_THRESHOLD = 4 };

} //namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_MEMORYTRAITS_HPP */

