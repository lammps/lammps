// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstddef>
#include <type_traits>

using Kokkos::Impl::BV::BasicView;
#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
static_assert(
    std::is_convertible_v<
        Kokkos::View<long long ****, Kokkos::LayoutRight, Kokkos::HostSpace>,
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
#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
static_assert(
    std::is_convertible_v<
        Kokkos::View<long long ****, Kokkos::LayoutRight, Kokkos::HostSpace>,
        BasicView<const long long, Kokkos::dextents<size_t, 4>,
                  Kokkos::Experimental::layout_right_padded<>,
                  Kokkos::Impl::CheckedReferenceCountedAccessor<
                      const long long, Kokkos::HostSpace>>>);
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
