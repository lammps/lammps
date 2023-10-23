//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXPERIMENTAL_VIEWHOOKS_HPP
#define KOKKOS_EXPERIMENTAL_VIEWHOOKS_HPP

namespace Kokkos {
namespace Experimental {

namespace Impl {
template <typename View>
using copy_subscription_function_type = void (*)(View &, const View &);

template <template <typename> class Invoker, typename... Subscribers>
struct invoke_subscriber_impl;

template <template <typename> class Invoker>
struct invoke_subscriber_impl<Invoker> {
  template <typename ViewType>
  static void invoke(ViewType &, const ViewType &) {}
};

template <template <typename> class Invoker, typename Subscriber,
          typename... RemSubscribers>
struct invoke_subscriber_impl<Invoker, Subscriber, RemSubscribers...> {
  template <typename ViewType>
  static void invoke(ViewType &self, const ViewType &other) {
    Invoker<Subscriber>::call(self, other);
    invoke_subscriber_impl<Invoker, RemSubscribers...>::invoke(self, other);
  }
};

template <typename Subscriber>
struct copy_constructor_invoker {
  template <typename View>
  static void call(View &self, const View &other) {
    Subscriber::copy_constructed(self, other);
  }
};

template <typename Subscriber>
struct move_constructor_invoker {
  template <typename View>
  static void call(View &self, const View &other) {
    Subscriber::move_constructed(self, other);
  }
};

template <typename Subscriber>
struct copy_assignment_operator_invoker {
  template <typename View>
  static void call(View &self, const View &other) {
    Subscriber::copy_assigned(self, other);
  }
};

template <typename Subscriber>
struct move_assignment_operator_invoker {
  template <typename View>
  static void call(View &self, const View &other) {
    Subscriber::move_assigned(self, other);
  }
};
}  // namespace Impl

struct EmptyViewHooks {
  using hooks_policy = EmptyViewHooks;

  template <typename View>
  static void copy_construct(View &, const View &) {}
  template <typename View>
  static void copy_assign(View &, const View &) {}
  template <typename View>
  static void move_construct(View &, const View &) {}
  template <typename View>
  static void move_assign(View &, const View &) {}
};

template <class... Subscribers>
struct SubscribableViewHooks {
  using hooks_policy = SubscribableViewHooks<Subscribers...>;

  template <typename View>
  static void copy_construct(View &self, const View &other) {
    Impl::invoke_subscriber_impl<Impl::copy_constructor_invoker,
                                 Subscribers...>::invoke(self, other);
  }
  template <typename View>
  static void copy_assign(View &self, const View &other) {
    Impl::invoke_subscriber_impl<Impl::copy_assignment_operator_invoker,
                                 Subscribers...>::invoke(self, other);
  }
  template <typename View>
  static void move_construct(View &self, const View &other) {
    Impl::invoke_subscriber_impl<Impl::move_constructor_invoker,
                                 Subscribers...>::invoke(self, other);
  }
  template <typename View>
  static void move_assign(View &self, const View &other) {
    Impl::invoke_subscriber_impl<Impl::move_assignment_operator_invoker,
                                 Subscribers...>::invoke(self, other);
  }
};

using DefaultViewHooks = EmptyViewHooks;

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_EXPERIMENTAL_VIEWHOOKS_HPP
