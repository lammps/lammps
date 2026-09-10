// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>  // show_warnings
#include <impl/Kokkos_Error.hpp>

void Kokkos::Impl::throw_runtime_exception(const std::string &msg) {
  throw std::runtime_error(msg);
}

void Kokkos::Impl::throw_bad_alloc(std::string_view memory_space_name,
                                   std::size_t size, std::string_view label) {
  throw Kokkos::Experimental::BadAlloc(memory_space_name, size, label);
}

void Kokkos::Impl::log_warning(const std::string &msg) {
  if (show_warnings()) {
    std::cerr << msg << std::flush;
  }
}

std::string Kokkos::Impl::human_memory_size(size_t arg_bytes) {
  double bytes   = arg_bytes;
  const double K = 1024;
  const double M = K * 1024;
  const double G = M * 1024;
  const double T = G * 1024;

  std::ostringstream out;
  if (bytes < K) {
    out << std::setprecision(4) << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << std::setprecision(4) << bytes << " KiB";
  } else if (bytes < G) {
    bytes /= M;
    out << std::setprecision(4) << bytes << " MiB";
  } else if (bytes < T) {
    bytes /= G;
    out << std::setprecision(4) << bytes << " GiB";
  } else {
    bytes /= T;
    out << std::setprecision(4) << bytes << " TiB";
  }
  return out.str();
}
