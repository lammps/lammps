// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
struct NotOnDeviceCtorDisambiguator {};

template <class... Args>
struct NoCtorsNotOnDevice : std::false_type {};

template <class... Args>
struct DefaultCtorNotOnDevice : std::false_type {};

template <>
struct DefaultCtorNotOnDevice<> : std::true_type {};

template <class T, bool Empty,
          template <class...> class CtorNotOnDevice = NoCtorsNotOnDevice>
struct EBOBaseImpl;

template <class T, template <class...> class CtorNotOnDevice>
struct EBOBaseImpl<T, true, CtorNotOnDevice> {
  template <class... Args, class _ignored = void,
            std::enable_if_t<std::is_void_v<_ignored> &&
                                 std::is_constructible_v<T, Args...> &&
                                 !CtorNotOnDevice<Args...>::value,
                             int> = 0>
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit EBOBaseImpl(
      Args&&...) noexcept {}

  template <class... Args, class _ignored = void,
            std::enable_if_t<std::is_void_v<_ignored> &&
                                 std::is_constructible_v<T, Args...> &&
                                 CtorNotOnDevice<Args...>::value,
                             long> = 0>
  inline constexpr explicit EBOBaseImpl(Args&&...) noexcept {}

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl& operator=(EBOBaseImpl const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl& operator=(EBOBaseImpl&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~EBOBaseImpl() = default;

  KOKKOS_INLINE_FUNCTION
  constexpr T& _ebo_data_member() & { return *reinterpret_cast<T*>(this); }

  KOKKOS_INLINE_FUNCTION
  constexpr T const& _ebo_data_member() const& {
    return *reinterpret_cast<T const*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  T volatile& _ebo_data_member() volatile& {
    return *reinterpret_cast<T volatile*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  T const volatile& _ebo_data_member() const volatile& {
    return *reinterpret_cast<T const volatile*>(this);
  }

  KOKKOS_INLINE_FUNCTION
  constexpr T&& _ebo_data_member() && {
    return std::move(*reinterpret_cast<T*>(this));
  }
};

template <class T, template <class...> class CTorsNotOnDevice>
struct EBOBaseImpl<T, false, CTorsNotOnDevice> {
  T m_ebo_object;

  // NOLINTBEGIN(modernize-type-traits)
  template <class... Args, class _ignored = void,
            std::enable_if_t<std::is_void_v<_ignored> &&
                                 !CTorsNotOnDevice<Args...>::value &&
                                 std::is_constructible_v<T, Args...>,
                             int> = 0>
  // NOLINTEND(modernize-type-traits)
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit EBOBaseImpl(
      Args&&... args) noexcept(noexcept(T(std::forward<Args>(args)...)))
      : m_ebo_object(std::forward<Args>(args)...) {}

  template <class... Args, class _ignored = void,
            std::enable_if_t<std::is_void_v<_ignored> &&
                                 CTorsNotOnDevice<Args...>::value &&
                                 std::is_constructible_v<T, Args...>,
                             long> = 0>
  inline constexpr explicit EBOBaseImpl(Args&&... args) noexcept(
      noexcept(T(std::forward<Args>(args)...)))
      : m_ebo_object(std::forward<Args>(args)...) {}

  // TODO @tasking @minor DSH noexcept in the right places?

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl(EBOBaseImpl&&) noexcept = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl& operator=(EBOBaseImpl const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr EBOBaseImpl& operator=(EBOBaseImpl&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~EBOBaseImpl() = default;

  KOKKOS_INLINE_FUNCTION
  T& _ebo_data_member() & { return m_ebo_object; }

  KOKKOS_INLINE_FUNCTION
  T const& _ebo_data_member() const& { return m_ebo_object; }

  KOKKOS_INLINE_FUNCTION
  T volatile& _ebo_data_member() volatile& { return m_ebo_object; }

  KOKKOS_INLINE_FUNCTION
  T const volatile& _ebo_data_member() const volatile& { return m_ebo_object; }

  KOKKOS_INLINE_FUNCTION
  T&& _ebo_data_member() && { return m_ebo_object; }
};

/**
 *
 * @tparam T
 */
template <class T,
          template <class...> class CtorsNotOnDevice = NoCtorsNotOnDevice>
struct StandardLayoutNoUniqueAddressMemberEmulation
    : EBOBaseImpl<T, std::is_empty_v<T>, CtorsNotOnDevice> {
 private:
  using ebo_base_t = EBOBaseImpl<T, std::is_empty_v<T>, CtorsNotOnDevice>;

 public:
  using ebo_base_t::ebo_base_t;

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr T& no_unique_address_data_member() & {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr T const& no_unique_address_data_member() const& {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  T volatile& no_unique_address_data_member() volatile& {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  T const volatile& no_unique_address_data_member() const volatile& {
    return this->ebo_base_t::_ebo_data_member();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr T&& no_unique_address_data_member() && {
    return this->ebo_base_t::_ebo_data_member();
  }
};

/**
 *
 * @tparam T
 */
template <class T,
          template <class...> class CtorsNotOnDevice = NoCtorsNotOnDevice>
class NoUniqueAddressMemberEmulation
    : private StandardLayoutNoUniqueAddressMemberEmulation<T,
                                                           CtorsNotOnDevice> {
 private:
  using base_t =
      StandardLayoutNoUniqueAddressMemberEmulation<T, CtorsNotOnDevice>;

 public:
  using base_t::base_t;
  using base_t::no_unique_address_data_member;
};

template <class InstanceType>
class InstanceStorage
    : private NoUniqueAddressMemberEmulation<InstanceType,
                                             DefaultCtorNotOnDevice> {
 private:
  using base_t =
      NoUniqueAddressMemberEmulation<InstanceType, DefaultCtorNotOnDevice>;

 protected:
  constexpr explicit InstanceStorage() : base_t() {}

  KOKKOS_INLINE_FUNCTION
  constexpr explicit InstanceStorage(InstanceType const& instance)
      : base_t(instance) {}

  KOKKOS_INLINE_FUNCTION
  constexpr explicit InstanceStorage(InstanceType&& instance)
      : base_t(std::move(instance)) {}

  KOKKOS_INLINE_FUNCTION
  InstanceType& instance() & { return this->no_unique_address_data_member(); }

  KOKKOS_INLINE_FUNCTION
  InstanceType const& instance() const& {
    return this->no_unique_address_data_member();
  }

  KOKKOS_INLINE_FUNCTION
  InstanceType&& instance() && {
    return std::move(*this).no_unique_address_data_member();
  }
};

}  // end namespace Impl
}  // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EBO_HPP */
