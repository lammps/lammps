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

#ifndef KOKKOS_EBO_HPP
#define KOKKOS_EBO_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>
//----------------------------------------------------------------------------


#include <utility>
#include <type_traits>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <int I>
struct NotOnDeviceCtorDisambiguator { };

template <class... Args>
struct NoCtorsNotOnDevice : std::false_type { };

template <class... Args>
struct DefaultCtorNotOnDevice : std::false_type { };

template <>
struct DefaultCtorNotOnDevice<> : std::true_type { };

template <class T, bool Empty, template <class...> class CtorNotOnDevice = NoCtorsNotOnDevice>
struct EBOBaseImpl;

template <class T, template <class...> class CtorNotOnDevice>
struct EBOBaseImpl<T, true, CtorNotOnDevice> {

  /*
   * Workaround for constexpr in C++11: we need to still call T(args...), but we
   * can't do so in the body of a constexpr function (in C++11), and there's no
   * data member to construct into. But we can construct into an argument
   * of a delegating constructor...
   */
  // TODO @minor DSH the destructor gets called too early with this workaround
  struct _constexpr_14_workaround_tag { };
  struct _constexpr_14_workaround_no_device_tag { };
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr EBOBaseImpl(_constexpr_14_workaround_tag, T&&) noexcept { }
  inline constexpr EBOBaseImpl(_constexpr_14_workaround_no_device_tag, T&&) noexcept { }

  template <
    class... Args,
    class _ignored = void,
    typename std::enable_if<
      std::is_void<_ignored>::value
      && std::is_constructible<T, Args...>::value
      && !CtorNotOnDevice<Args...>::value,
      int
    >::type = 0
  >
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr explicit
  EBOBaseImpl(
    Args&&... args
  ) noexcept(noexcept(T(std::forward<Args>(args)...)))
    // still call the constructor
    : EBOBaseImpl(_constexpr_14_workaround_tag{}, T(std::forward<Args>(args)...))
  { }

  template <
    class... Args,
    class _ignored=void,
    typename std::enable_if<
      std::is_void<_ignored>::value
      && std::is_constructible<T, Args...>::value
      && CtorNotOnDevice<Args...>::value,
      long
    >::type = 0
  >
  inline constexpr explicit
  EBOBaseImpl(
    Args&&... args
  ) noexcept(noexcept(T(std::forward<Args>(args)...)))
    // still call the constructor
    : EBOBaseImpl(_constexpr_14_workaround_no_device_tag{}, T(std::forward<Args>(args)...))
  { }

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl const&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl&&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  EBOBaseImpl& operator=(EBOBaseImpl const&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  EBOBaseImpl& operator=(EBOBaseImpl&&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  ~EBOBaseImpl() = default;

  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  T& _ebo_data_member() & {
    return *reinterpret_cast<T*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  constexpr
  T const& _ebo_data_member() const & {
    return *reinterpret_cast<T const*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  T volatile& _ebo_data_member() volatile & {
    return *reinterpret_cast<T volatile*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  T const volatile& _ebo_data_member() const volatile & {
    return *reinterpret_cast<T const volatile*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  T&& _ebo_data_member() && {
    return std::move(*reinterpret_cast<T*>(this));
  }

};

template <class T, template <class...> class CTorsNotOnDevice>
struct EBOBaseImpl<T, false, CTorsNotOnDevice> {

  T m_ebo_object;

  template <
    class... Args,
    class _ignored=void,
    typename std::enable_if<
      std::is_void<_ignored>::value
        && !CTorsNotOnDevice<Args...>::value
        && std::is_constructible<T, Args...>::value,
      int
    >::type = 0
  >
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr explicit
  EBOBaseImpl(
    Args&&... args
  ) noexcept(noexcept(T(std::forward<Args>(args)...)))
    : m_ebo_object(std::forward<Args>(args)...)
  { }

  template <
    class... Args,
    class _ignored=void,
    typename std::enable_if<
      std::is_void<_ignored>::value
        && CTorsNotOnDevice<Args...>::value
        && std::is_constructible<T, Args...>::value,
      long
    >::type = 0
  >
  inline
  constexpr explicit
  EBOBaseImpl(
    Args&&... args
  ) noexcept(noexcept(T(std::forward<Args>(args)...)))
    : m_ebo_object(std::forward<Args>(args)...)
  { }


  // TODO @tasking @minor DSH noexcept in the right places?

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr
  EBOBaseImpl(EBOBaseImpl const&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr
  EBOBaseImpl(EBOBaseImpl&&) noexcept = default;

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  EBOBaseImpl& operator=(EBOBaseImpl const&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  EBOBaseImpl& operator=(EBOBaseImpl&&) = default;

  KOKKOS_FORCEINLINE_FUNCTION
  ~EBOBaseImpl() = default;

  KOKKOS_INLINE_FUNCTION
  T& _ebo_data_member() & {
    return m_ebo_object;
  }

  KOKKOS_INLINE_FUNCTION
  T const& _ebo_data_member() const & {
    return m_ebo_object;
  }

  KOKKOS_INLINE_FUNCTION
  T volatile& _ebo_data_member() volatile & {
    return m_ebo_object;
  }

  KOKKOS_INLINE_FUNCTION
  T const volatile& _ebo_data_member() const volatile & {
    return m_ebo_object;
  }

  KOKKOS_INLINE_FUNCTION
  T&& _ebo_data_member() && {
    return m_ebo_object;
  }

};

/**
 *
 * @tparam T
 */
template <class T, template <class...> class CtorsNotOnDevice=NoCtorsNotOnDevice>
struct StandardLayoutNoUniqueAddressMemberEmulation
  : EBOBaseImpl<T, std::is_empty<T>::value, CtorsNotOnDevice>
{
private:

  using ebo_base_t = EBOBaseImpl<T, std::is_empty<T>::value, CtorsNotOnDevice>;

public:

  using ebo_base_t::ebo_base_t;

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  T& no_unique_address_data_member() & {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr
  T const& no_unique_address_data_member() const & {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  T volatile& no_unique_address_data_member() volatile & {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  T const volatile& no_unique_address_data_member() const volatile & {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  T&& no_unique_address_data_member() && {
    return this->ebo_base_t::_ebo_data_member();
  }
};

/**
 *
 * @tparam T
 */
template <class T, template <class...> class CtorsNotOnDevice=NoCtorsNotOnDevice>
class NoUniqueAddressMemberEmulation
  : private StandardLayoutNoUniqueAddressMemberEmulation<T, CtorsNotOnDevice>
{
private:

  using base_t = StandardLayoutNoUniqueAddressMemberEmulation<T, CtorsNotOnDevice>;

public:

  using base_t::base_t;
  using base_t::no_unique_address_data_member;

};


} // end namespace Impl
} // end namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#endif /* #ifndef KOKKOS_EBO_HPP */

