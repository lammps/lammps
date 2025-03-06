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

#ifndef KOKKOS_COMMAND_LINE_PARSING_HPP
#define KOKKOS_COMMAND_LINE_PARSING_HPP

#include <string>
#include <regex>

namespace Kokkos {
namespace Impl {
bool is_unsigned_int(const char* str);
bool check_arg(char const* arg, char const* expected);
bool check_arg_bool(char const* arg, char const* name, bool& val);
bool check_arg_int(char const* arg, char const* name, int& val);
bool check_arg_str(char const* arg, char const* name, std::string& val);
bool check_env_bool(char const* name, bool& val);
bool check_env_int(char const* name, int& val);
void warn_deprecated_environment_variable(std::string deprecated);
void warn_deprecated_environment_variable(std::string deprecated,
                                          std::string use_instead);
void warn_deprecated_command_line_argument(std::string deprecated);
void warn_deprecated_command_line_argument(std::string deprecated,
                                           std::string use_instead);
void warn_not_recognized_command_line_argument(std::string not_recognized);
void do_not_warn_not_recognized_command_line_argument(std::regex ignore);
}  // namespace Impl
}  // namespace Kokkos

#endif
