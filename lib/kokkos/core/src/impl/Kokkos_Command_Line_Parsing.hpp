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
