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
#ifndef KOKKOS_UNIQUE_TOKEN_HPP
#define KOKKOS_UNIQUE_TOKEN_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core_fwd.hpp>

namespace Kokkos {
namespace Experimental {

enum class UniqueTokenScope : int { Instance, Global };

/// \brief class to generate unique ids base on the required amount of
/// concurrency
///
/// This object should behave like a ref-counted object, so that when the last
/// instance is destroy resources are free if needed
template <typename ExecutionSpace = Kokkos::DefaultExecutionSpace,
          UniqueTokenScope        = UniqueTokenScope::Instance>
class UniqueToken {
 public:
  using execution_space = ExecutionSpace;
  using size_type       = typename execution_space::size_type;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space());

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type size() const;

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type acquire() const;

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(size_type) const;
};

/// \brief Instance scope UniqueToken allows for a max size other than
/// execution_space::concurrency()
///
/// This object should behave like a ref-counted object, so that when the last
/// instance is destroyed, resources are free if needed
template <typename ExecutionSpace>
class UniqueToken<ExecutionSpace, UniqueTokenScope::Instance>
    : public UniqueToken<ExecutionSpace, UniqueTokenScope::Global> {
 public:
  using execution_space = ExecutionSpace;
  using size_type       = typename execution_space::size_type;

  /// \brief Create object with specified size
  ///
  /// It is required that max_size is >= the maximum number of concurrent
  /// threads that will attempt to acquire the UniqueToken. This constructor is
  /// most commonly useful when you:
  ///   1) Have a loop bound that may be smaller than
  ///   execution_space::concurrency().
  ///   2) Want a per-team unique token in the range [0,
  ///   execution_space::concurrency() / team_size)
  UniqueToken(size_type max_size, execution_space const& = execution_space());
};

// NOTE There was an agreement amongst developers that "AcquireUniqueToken" is a
// bad name but at this time no one has suggested a better alternative.

/// \brief RAII helper for per-thread unique token values.
///
/// The token value will be acquired at construction and automatically
/// released at destruction.
template <typename ExecutionSpace,
          UniqueTokenScope TokenScope = UniqueTokenScope::Instance>
class AcquireUniqueToken {
 public:
  using exec_space = ExecutionSpace;
  using size_type  = typename exec_space::size_type;
  using token_type = UniqueToken<exec_space, TokenScope>;

 private:
  token_type my_token;
  size_type my_acquired_val;

 public:
  KOKKOS_FUNCTION AcquireUniqueToken(token_type t)
      : my_token(t), my_acquired_val(my_token.acquire()) {}

  KOKKOS_FUNCTION ~AcquireUniqueToken() { my_token.release(my_acquired_val); }

  KOKKOS_FUNCTION size_type value() const { return my_acquired_val; }
};

/// \brief RAII helper for per-team unique token values.
///
/// The token value will be acquired at construction and automatically
/// released at destruction. All threads in a team will share the same
/// token value.
template <typename TeamPolicy>
class AcquireTeamUniqueToken {
 public:
  using exec_space       = typename TeamPolicy::execution_space;
  using token_type       = UniqueToken<exec_space>;
  using size_type        = typename token_type::size_type;
  using team_member_type = typename TeamPolicy::member_type;
  using scratch_view =
      Kokkos::View<size_type, typename exec_space::scratch_memory_space,
                   Kokkos::MemoryUnmanaged>;

 private:
  token_type my_token;
  size_type my_acquired_val;
  scratch_view my_team_acquired_val;
  team_member_type my_team;

 public:
  // NOTE The implementations of the constructor and destructor use
  // `Kokkos::single()` which is an inline function defined in each backend.
  // This creates circular dependency issues.  Moving them to a separate header
  // is less than ideal and should be revisited later.  Having a `UniqueToken`
  // forward declaration was considered but the non-type template parameter
  // makes things complicated because it would require moving the definition of
  // `UniqueTokenScope` enumeration type and its enumerators away which would
  // hurt readability.
  KOKKOS_FUNCTION AcquireTeamUniqueToken(token_type t, team_member_type team);
  KOKKOS_FUNCTION ~AcquireTeamUniqueToken();
  KOKKOS_FUNCTION size_type value() const { return my_acquired_val; }
  static std::size_t shmem_size() { return scratch_view::shmem_size(); }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_UNIQUE_TOKEN_HPP
