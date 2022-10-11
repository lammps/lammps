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

#ifndef KOKKOS_SERIAL_UNIQUE_TOKEN_HPP
#define KOKKOS_SERIAL_UNIQUE_TOKEN_HPP

#include <Kokkos_UniqueToken.hpp>

namespace Kokkos {
namespace Experimental {

template <>
class UniqueToken<Serial, UniqueTokenScope::Instance> {
 public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief create object size for requested size on given instance
  ///
  /// It is the users responsibility to only acquire size tokens concurrently
  UniqueToken(size_type, execution_space const& = execution_space()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

template <>
class UniqueToken<Serial, UniqueTokenScope::Global> {
 public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
