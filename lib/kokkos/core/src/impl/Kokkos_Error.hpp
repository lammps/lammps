// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_ERROR_HPP
#define KOKKOS_IMPL_ERROR_HPP

#include <string>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Abort.hpp>
#include <Kokkos_Assert.hpp>

namespace Kokkos::Impl {

[[noreturn]] void throw_runtime_exception(const std::string &msg);
[[noreturn]] void throw_bad_alloc(std::string_view memory_space_name,
                                  std::size_t size, std::string_view label);
void log_warning(const std::string &msg);

std::string human_memory_size(size_t bytes);

}  // namespace Kokkos::Impl

#endif
