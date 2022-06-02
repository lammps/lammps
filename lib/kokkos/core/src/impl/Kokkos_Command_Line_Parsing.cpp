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
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <impl/Kokkos_Command_Line_Parsing.hpp>
/** Duplicates of Kokkos_Error.cpp/hpp, reproduced here
 * for use in non-Kokkos applications
 */
namespace {
void traceback_callstack(std::ostream& msg) {
  msg << std::endl << "Traceback functionality not available" << std::endl;
}
void throw_runtime_exception(const std::string& msg) {
  std::ostringstream o;
  o << msg;
  traceback_callstack(o);
  throw std::runtime_error(o.str());
}
}  // namespace

bool Kokkos::Impl::is_unsigned_int(const char* str) {
  const size_t len = strlen(str);
  for (size_t i = 0; i < len; ++i) {
    if (!isdigit(str[i])) {
      return false;
    }
  }
  return true;
}

bool Kokkos::Impl::check_arg(char const* arg, char const* expected) {
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  if (arg_len < exp_len) return false;
  if (std::strncmp(arg, expected, exp_len) != 0) return false;
  if (arg_len == exp_len) return true;

  if (std::isalnum(arg[exp_len]) || arg[exp_len] == '-' ||
      arg[exp_len] == '_') {
    return false;
  }
  return true;
}

bool Kokkos::Impl::check_int_arg(char const* arg, char const* expected,
                                 int* value) {
  if (!check_arg(arg, expected)) return false;
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  bool okay           = true;
  if (arg_len == exp_len || arg[exp_len] != '=') okay = false;
  char const* number = arg + exp_len + 1;
  if (!Kokkos::Impl::is_unsigned_int(number) || strlen(number) == 0)
    okay = false;
  *value = std::stoi(number);
  if (!okay) {
    std::ostringstream ss;
    ss << "Error: expecting an '=INT' after command line argument '" << expected
       << "'";
    ss << ". Raised by Kokkos::initialize(int narg, char* argc[]).";
    throw_runtime_exception(ss.str());
  }
  return true;
}
bool Kokkos::Impl::check_str_arg(char const* arg, char const* expected,
                                 std::string& value) {
  if (!check_arg(arg, expected)) return false;
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  bool okay           = true;
  if (arg_len == exp_len || arg[exp_len] != '=') okay = false;
  char const* remain = arg + exp_len + 1;
  value              = remain;
  if (!okay) {
    std::ostringstream ss;
    ss << "Error: expecting an '=STRING' after command line argument '"
       << expected << "'";
    ss << ". Raised by Kokkos::initialize(int narg, char* argc[]).";
    throw_runtime_exception(ss.str());
  }
  return true;
}
void Kokkos::Impl::warn_deprecated_command_line_argument(std::string deprecated,
                                                         std::string valid) {
  std::cerr
      << "Warning: command line argument '" << deprecated
      << "' is deprecated. Use '" << valid
      << "' instead. Raised by Kokkos::initialize(int narg, char* argc[])."
      << std::endl;
}
