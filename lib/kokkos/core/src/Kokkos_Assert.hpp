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

#ifndef KOKKOS_ASSERT_HPP
#define KOKKOS_ASSERT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Abort.hpp>

#if !defined(NDEBUG) || defined(KOKKOS_ENFORCE_CONTRACTS) || \
    defined(KOKKOS_ENABLE_DEBUG)
#define KOKKOS_EXPECTS(...)                                                    \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Expected precondition `" #__VA_ARGS__                             \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
#define KOKKOS_ENSURES(...)                                                    \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Ensured postcondition `" #__VA_ARGS__                             \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
// some projects already define this for themselves, so don't mess
// them up
#ifndef KOKKOS_ASSERT
#define KOKKOS_ASSERT(...)                                                     \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Asserted condition `" #__VA_ARGS__                                \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
#endif  // ifndef KOKKOS_ASSERT
#else   // not debug mode
#define KOKKOS_EXPECTS(...)
#define KOKKOS_ENSURES(...)
#ifndef KOKKOS_ASSERT
#define KOKKOS_ASSERT(...)
#endif  // ifndef KOKKOS_ASSERT
#endif  // end debug mode ifdefs

#endif /* #ifndef KOKKOS_ASSERT_HPP */
