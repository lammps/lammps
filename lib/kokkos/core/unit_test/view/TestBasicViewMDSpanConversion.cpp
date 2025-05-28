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

#include <Kokkos_Core.hpp>
#include <type_traits>

using Kokkos::Impl::BV::BasicView;
#if 0  // TODO: after View is using BasicView this should be true
static_assert(
    std::is_convertible_v<
        Kokkos::View<long long ****, Kokkos::LayoutRight, Kokkos::Serial>,
        BasicView<long long, Kokkos::dextents<size_t, 4>,
                  Kokkos::Experimental::layout_right_padded<>,
                  Kokkos::Impl::CheckedReferenceCountedAccessor<
                    long long, Kokkos::HostSpace>>>);
#endif

static_assert(std::is_convertible_v<
              BasicView<long long, Kokkos::dextents<size_t, 4>,
                        Kokkos::Experimental::layout_right_padded<>,
                        Kokkos::Impl::CheckedReferenceCountedAccessor<
                            long long, Kokkos::HostSpace>>,
              BasicView<const long long, Kokkos::dextents<size_t, 4>,
                        Kokkos::Experimental::layout_right_padded<>,
                        Kokkos::Impl::CheckedReferenceCountedAccessor<
                            const long long, Kokkos::HostSpace>>>);
#if 0  // TODO: after View is using BasicView this should be true
static_assert(
    std::is_convertible_v<
        Kokkos::View<long long ****, Kokkos::LayoutRight, Kokkos::Serial>,
        BasicView<const long long, Kokkos::dextents<size_t, 4>,
                  Kokkos::Experimental::layout_right_padded<>,
                  Kokkos::Impl::CheckedReferenceCountedAccessor<
                    const long long, Kokkos::HostSpace>>>);

using test_atomic_view = Kokkos::View<double *, Kokkos::Serial,
                                      Kokkos::MemoryTraits<Kokkos::Atomic>>;
static_assert(std::is_same_v<
              decltype(std::declval<test_atomic_view>()(std::declval<int>())),
              desul::AtomicRef<double, desul::MemoryOrderRelaxed,
                               desul::MemoryScopeDevice>>);
#endif

static_assert(std::is_convertible_v<Kokkos::default_accessor<double>,
                                    Kokkos::Impl::ReferenceCountedAccessor<
                                        double, Kokkos::HostSpace,
                                        Kokkos::default_accessor<double>>>);

static_assert(std::is_constructible_v<Kokkos::default_accessor<const double>,
                                      Kokkos::default_accessor<double>>);

static_assert(std::is_convertible_v<Kokkos::default_accessor<double>,
                                    Kokkos::default_accessor<const double>>);

static_assert(
    std::is_constructible_v<
        Kokkos::Impl::ReferenceCountedAccessor<
            const double, Kokkos::HostSpace,
            Kokkos::default_accessor<const double>>,
        Kokkos::Impl::ReferenceCountedAccessor<
            double, Kokkos::HostSpace, Kokkos::default_accessor<double>>>);

static_assert(std::is_convertible_v<
              Kokkos::Impl::ReferenceCountedAccessor<
                  double, Kokkos::HostSpace, Kokkos::default_accessor<double>>,
              Kokkos::Impl::ReferenceCountedAccessor<
                  const double, Kokkos::HostSpace,
                  Kokkos::default_accessor<const double>>>);

static_assert(std::is_constructible_v<Kokkos::default_accessor<const double>,
                                      Kokkos::Impl::ReferenceCountedAccessor<
                                          double, Kokkos::HostSpace,
                                          Kokkos::default_accessor<double>>>);

static_assert(
    std::is_convertible_v<
        Kokkos::Impl::SpaceAwareAccessor<
            Kokkos::HostSpace,
            Kokkos::Impl::ReferenceCountedAccessor<
                double, Kokkos::HostSpace, Kokkos::default_accessor<double>>>,
        Kokkos::Impl::SpaceAwareAccessor<
            Kokkos::HostSpace, Kokkos::default_accessor<const double>>>);
