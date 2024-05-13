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

#ifndef KOKKOS_OPENACC_INSTANCE_HPP
#define KOKKOS_OPENACC_INSTANCE_HPP

#include <impl/Kokkos_InitializationSettings.hpp>

#include <openacc.h>

#include <cstdint>
#include <iosfwd>
#include <string>

namespace Kokkos::Experimental::Impl {

class OpenACCInternal {
  bool m_is_initialized = false;

  OpenACCInternal(const OpenACCInternal&) = default;
  OpenACCInternal& operator=(const OpenACCInternal&) = default;

 public:
  static int m_acc_device_num;
  int m_async_arg = acc_async_noval;

  OpenACCInternal() = default;

  static OpenACCInternal& singleton();

  bool verify_is_initialized(const char* const label) const;

  void initialize(int async_arg = acc_async_noval);
  void finalize();
  bool is_initialized() const;

  void print_configuration(std::ostream& os, bool verbose = false) const;

  void fence(std::string const& name) const;

  uint32_t instance_id() const noexcept;
};

}  // namespace Kokkos::Experimental::Impl

#endif
