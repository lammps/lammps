/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
#endif
#ifndef KOKKOS_MEMORYTRAITS_HPP
#define KOKKOS_MEMORYTRAITS_HPP

#include <impl/Kokkos_Traits.hpp>

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
enum MemoryTraitsFlags {
  Unmanaged    = 0x01,
  RandomAccess = 0x02,
  Atomic       = 0x04,
  Restrict     = 0x08,
  Aligned      = 0x10
};

template <unsigned T>
struct MemoryTraits {
  //! Tag this class as a kokkos memory traits:
  using memory_traits = MemoryTraits<T>;
  enum : bool {
    is_unmanaged = (unsigned(0) != (T & unsigned(Kokkos::Unmanaged)))
  };
  enum : bool {
    is_random_access = (unsigned(0) != (T & unsigned(Kokkos::RandomAccess)))
  };
  enum : bool { is_atomic = (unsigned(0) != (T & unsigned(Kokkos::Atomic))) };
  enum : bool {
    is_restrict = (unsigned(0) != (T & unsigned(Kokkos::Restrict)))
  };
  enum : bool { is_aligned = (unsigned(0) != (T & unsigned(Kokkos::Aligned))) };
};

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

using MemoryManaged   = Kokkos::MemoryTraits<0>;
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
using MemoryRandomAccess =
    Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>;

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

static_assert((0 < int(KOKKOS_MEMORY_ALIGNMENT)) &&
                  (0 == (int(KOKKOS_MEMORY_ALIGNMENT) &
                         (int(KOKKOS_MEMORY_ALIGNMENT) - 1))),
              "KOKKOS_MEMORY_ALIGNMENT must be a power of two");

/** \brief Memory alignment settings
 *
 *  Sets global value for memory alignment.  Must be a power of two!
 *  Enable compatibility of views from different devices with static stride.
 *  Use compiler flag to enable overwrites.
 */
enum : unsigned {
  MEMORY_ALIGNMENT           = KOKKOS_MEMORY_ALIGNMENT,
  MEMORY_ALIGNMENT_THRESHOLD = KOKKOS_MEMORY_ALIGNMENT_THRESHOLD
};

// ------------------------------------------------------------------ //
//  this identifies the default memory trait
//
template <typename Tp>
struct is_default_memory_trait : std::false_type {};

template <>
struct is_default_memory_trait<Kokkos::MemoryTraits<0>> : std::true_type {};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_MEMORYTRAITS_HPP */
