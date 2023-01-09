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
