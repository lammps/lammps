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

#ifndef KOKKOS_ANONYMOUSSPACE_HPP
#define KOKKOS_ANONYMOUSSPACE_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>
#include <cstddef>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class AnonymousSpace {
public:
  //! Tag this class as a kokkos memory space
  typedef AnonymousSpace  memory_space;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef size_t     size_type;

  //! This memory space preferred device_type
  typedef Kokkos::Device< execution_space, memory_space > device_type;

  /**\brief  Default memory space instance */
  AnonymousSpace() = default;
  AnonymousSpace( AnonymousSpace && rhs ) = default;
  AnonymousSpace( const AnonymousSpace & rhs ) = default;
  AnonymousSpace & operator = ( AnonymousSpace && ) = default;
  AnonymousSpace & operator = ( const AnonymousSpace & ) = default;
  ~AnonymousSpace() = default;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return "Anonymous"; }
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template<typename OtherSpace>
struct MemorySpaceAccess< Kokkos::AnonymousSpace , OtherSpace > {
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy   = true };
};

template<typename OtherSpace>
struct MemorySpaceAccess<  OtherSpace, Kokkos::AnonymousSpace > {
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy   = true };
};

template<>
struct MemorySpaceAccess<  Kokkos::AnonymousSpace, Kokkos::AnonymousSpace > {
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy   = true };
};

template<typename OtherSpace>
struct VerifyExecutionCanAccessMemorySpace<OtherSpace, Kokkos::AnonymousSpace>
{
  enum {value = 1};
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void *) {}
};

template<typename OtherSpace>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::AnonymousSpace, OtherSpace>
{
  enum {value = 1};
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void *) {}
};

} // namespace Impl

} // namespace Kokkos

#endif // #define KOKKOS_ANONYMOUSSPACE_HPP

