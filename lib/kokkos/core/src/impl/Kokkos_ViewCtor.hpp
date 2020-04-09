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

#ifndef KOKKOS_EXPERIMENTAL_IMPL_VIEW_CTOR_PROP_HPP
#define KOKKOS_EXPERIMENTAL_IMPL_VIEW_CTOR_PROP_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/* For backward compatibility */

struct ViewAllocateWithoutInitializing {
  const std::string label;

  ViewAllocateWithoutInitializing() : label() {}

  explicit ViewAllocateWithoutInitializing(const std::string &arg_label)
      : label(arg_label) {}

  explicit ViewAllocateWithoutInitializing(const char *const arg_label)
      : label(arg_label) {}
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct WithoutInitializing_t {};
struct AllowPadding_t {};
struct NullSpace_t {};

//----------------------------------------------------------------------------
/**\brief Whether a type can be used for a view label */

template <typename>
struct is_view_label : public std::false_type {};

template <>
struct is_view_label<std::string> : public std::true_type {};

template <unsigned N>
struct is_view_label<char[N]> : public std::true_type {};

template <unsigned N>
struct is_view_label<const char[N]> : public std::true_type {};

//----------------------------------------------------------------------------

template <typename... P>
struct ViewCtorProp;

// Forward declare
template <typename Specialize, typename T>
struct CommonViewAllocProp;

/* Common value_type stored as ViewCtorProp
 */
template <typename Specialize, typename T>
struct ViewCtorProp<void, CommonViewAllocProp<Specialize, T> > {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  using type = CommonViewAllocProp<Specialize, T>;

  KOKKOS_INLINE_FUNCTION
  ViewCtorProp(const type &arg) : value(arg) {}
  KOKKOS_INLINE_FUNCTION
  ViewCtorProp(type &&arg) : value(arg) {}

  type value;
};

/*  std::integral_constant<unsigned,I> are dummy arguments
 *  that avoid duplicate base class errors
 */
template <unsigned I>
struct ViewCtorProp<void, std::integral_constant<unsigned, I> > {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  template <typename P>
  KOKKOS_INLINE_FUNCTION ViewCtorProp(const P &) {}
};

/* Property flags have constexpr value */
template <typename P>
struct ViewCtorProp<typename std::enable_if<
                        std::is_same<P, AllowPadding_t>::value ||
                        std::is_same<P, WithoutInitializing_t>::value>::type,
                    P> {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  typedef P type;

  ViewCtorProp(const type &) {}

  static constexpr type value = type();
};

/* Map input label type to std::string */
template <typename Label>
struct ViewCtorProp<typename std::enable_if<is_view_label<Label>::value>::type,
                    Label> {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  typedef std::string type;

  ViewCtorProp(const type &arg) : value(arg) {}
  ViewCtorProp(type &&arg) : value(arg) {}

  type value;
};

template <typename Space>
struct ViewCtorProp<typename std::enable_if<
                        Kokkos::Impl::is_memory_space<Space>::value ||
                        Kokkos::Impl::is_execution_space<Space>::value>::type,
                    Space> {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  typedef Space type;

  ViewCtorProp(const type &arg) : value(arg) {}

  type value;
};

template <typename T>
struct ViewCtorProp<void, T *> {
  ViewCtorProp()                     = default;
  ViewCtorProp(const ViewCtorProp &) = default;
  ViewCtorProp &operator=(const ViewCtorProp &) = default;

  typedef T *type;

  KOKKOS_INLINE_FUNCTION
  ViewCtorProp(const type arg) : value(arg) {}

  type value;
};

template <typename... P>
struct ViewCtorProp : public ViewCtorProp<void, P>... {
 private:
  typedef Kokkos::Impl::has_condition<void, Kokkos::Impl::is_memory_space, P...>
      var_memory_space;

  typedef Kokkos::Impl::has_condition<void, Kokkos::Impl::is_execution_space,
                                      P...>
      var_execution_space;

  struct VOIDDUMMY {};

  typedef Kokkos::Impl::has_condition<VOIDDUMMY, std::is_pointer, P...>
      var_pointer;

 public:
  /* Flags for the common properties */
  enum { has_memory_space = var_memory_space::value };
  enum { has_execution_space = var_execution_space::value };
  enum { has_pointer = var_pointer::value };
  enum { has_label = Kokkos::Impl::has_type<std::string, P...>::value };
  enum { allow_padding = Kokkos::Impl::has_type<AllowPadding_t, P...>::value };
  enum {
    initialize = !Kokkos::Impl::has_type<WithoutInitializing_t, P...>::value
  };

  typedef typename var_memory_space::type memory_space;
  typedef typename var_execution_space::type execution_space;
  typedef typename var_pointer::type pointer_type;

  /*  Copy from a matching argument list.
   *  Requires  std::is_same< P , ViewCtorProp< void , Args >::value ...
   */
  template <typename... Args>
  inline ViewCtorProp(Args const &... args) : ViewCtorProp<void, P>(args)... {}

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION ViewCtorProp(pointer_type arg0, Args const &... args)
      : ViewCtorProp<void, pointer_type>(arg0),
        ViewCtorProp<void, typename ViewCtorProp<void, Args>::type>(args)... {}

  /* Copy from a matching property subset */
  template <typename... Args>
  ViewCtorProp(ViewCtorProp<Args...> const &arg)
      : ViewCtorProp<void, Args>(((ViewCtorProp<void, Args> const &)arg))... {}
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
