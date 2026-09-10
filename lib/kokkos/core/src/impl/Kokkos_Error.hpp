// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_ERROR_HPP
#define KOKKOS_IMPL_ERROR_HPP

#include <string>
#include <string_view>
#include <stdexcept>
#include <sstream>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Abort.hpp>
#include <Kokkos_Assert.hpp>

namespace Kokkos::Impl {

[[noreturn]] void throw_runtime_exception(const std::string& msg);
[[noreturn]] void throw_bad_alloc(std::string_view memory_space_name,
                                  std::size_t size, std::string_view label);
void log_warning(const std::string& msg);

std::string human_memory_size(size_t bytes);

}  // namespace Kokkos::Impl

namespace Kokkos::Experimental {

/// Exception thrown when a Kokkos memory allocation fails.
///
/// Inherits from std::runtime_error for backward compatibility — existing
/// catch blocks that catch std::runtime_error still work.  Provides structured
/// access to the memory space, requested allocation size, and label of the
/// allocation that failed.
class BadAlloc : public std::runtime_error {
  std::string m_memory_space_name;
  std::size_t m_allocation_size;
  std::string m_label;

  static std::string make_message(std::string_view memory_space_name,
                                  std::size_t size, std::string_view label) {
    std::stringstream ss;
    ss << "Kokkos ERROR: " << memory_space_name
       << " memory space failed to allocate "
       << Kokkos::Impl::human_memory_size(size) << " (label=\"" << label
       << "\").";
    return ss.str();
  }

 public:
  BadAlloc(std::string_view memory_space_name, std::size_t size,
           std::string_view label)
      : std::runtime_error(make_message(memory_space_name, size, label)),
        m_memory_space_name(memory_space_name),
        m_allocation_size(size),
        m_label(label) {}

  const std::string& memory_space_name() const noexcept {
    return m_memory_space_name;
  }
  std::size_t allocation_size() const noexcept { return m_allocation_size; }
  const std::string& label() const noexcept { return m_label; }
};

}  // namespace Kokkos::Experimental

#endif
