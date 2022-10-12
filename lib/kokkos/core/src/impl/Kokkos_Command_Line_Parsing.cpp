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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <impl/Kokkos_Command_Line_Parsing.hpp>
#include <impl/Kokkos_Error.hpp>

#include <cstring>
#include <iostream>
#include <regex>
#include <string>
#include <sstream>

namespace {

auto const regex_true = std::regex(
    "(yes|true|1)", std::regex_constants::icase | std::regex_constants::egrep);

auto const regex_false = std::regex(
    "(no|false|0)", std::regex_constants::icase | std::regex_constants::egrep);

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

bool Kokkos::Impl::check_env_bool(char const* name, bool& val) {
  char const* var = std::getenv(name);

  if (!var) {
    return false;
  }

  if (std::regex_match(var, regex_true)) {
    val = true;
    return true;
  }

  if (!std::regex_match(var, regex_false)) {
    std::stringstream ss;
    ss << "Error: cannot convert environment variable '" << name << "=" << var
       << "' to a boolean."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  val = false;
  return true;
}

bool Kokkos::Impl::check_env_int(char const* name, int& val) {
  char const* var = std::getenv(name);

  if (!var) {
    return false;
  }

  errno = 0;
  char* var_end;
  val = std::strtol(var, &var_end, 10);

  if (var == var_end) {
    std::stringstream ss;
    ss << "Error: cannot convert environment variable '" << name << '=' << var
       << "' to an integer."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for environment variable '" << name << '='
       << var << "' falls out of range."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  return true;
}

bool Kokkos::Impl::check_arg_bool(char const* arg, char const* name,
                                  bool& val) {
  auto const len = std::strlen(name);
  if (std::strncmp(arg, name, len) != 0) {
    return false;
  }
  auto const arg_len = strlen(arg);
  if (arg_len == len) {
    val = true;  // --kokkos-foo without =BOOL interpreted as fool=true
    return true;
  }
  if (arg_len <= len + 1 || arg[len] != '=') {
    std::stringstream ss;
    ss << "Error: command line argument '" << arg
       << "' is not recognized as a valid boolean."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  std::advance(arg, len + 1);
  if (std::regex_match(arg, regex_true)) {
    val = true;
    return true;
  }
  if (!std::regex_match(arg, regex_false)) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << name << "=" << arg
       << "' to a boolean."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }
  val = false;
  return true;
}

bool Kokkos::Impl::check_arg_int(char const* arg, char const* name, int& val) {
  auto const len = std::strlen(name);
  if (std::strncmp(arg, name, len) != 0) {
    return false;
  }
  auto const arg_len = strlen(arg);
  if (arg_len <= len + 1 || arg[len] != '=') {
    std::stringstream ss;
    ss << "Error: command line argument '" << arg
       << "' is not recognized as a valid integer."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  std::advance(arg, len + 1);

  errno = 0;
  char* arg_end;
  val = std::strtol(arg, &arg_end, 10);

  if (arg == arg_end) {
    std::stringstream ss;
    ss << "Error: cannot convert command line argument '" << name << '=' << arg
       << "' to an integer."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  if (errno == ERANGE) {
    std::stringstream ss;
    ss << "Error: converted value for command line argument '" << name << '='
       << arg << "' falls out of range."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  return true;
}

bool Kokkos::Impl::check_arg_str(char const* arg, char const* name,
                                 std::string& val) {
  auto const len = std::strlen(name);
  if (std::strncmp(arg, name, len) != 0) {
    return false;
  }
  auto const arg_len = strlen(arg);
  if (arg_len <= len + 1 || arg[len] != '=') {
    std::stringstream ss;
    ss << "Error: command line argument '" << arg
       << "' is not recognized as a valid string."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  std::advance(arg, len + 1);

  val = arg;
  return true;
}

void Kokkos::Impl::warn_deprecated_environment_variable(
    std::string deprecated) {
  std::cerr << "Warning: environment variable '" << deprecated
            << "' is deprecated."
            << " Raised by Kokkos::initialize()." << std::endl;
}

void Kokkos::Impl::warn_deprecated_environment_variable(
    std::string deprecated, std::string use_instead) {
  std::cerr << "Warning: environment variable '" << deprecated
            << "' is deprecated."
            << " Use '" << use_instead << "' instead."
            << " Raised by Kokkos::initialize()." << std::endl;
}

void Kokkos::Impl::warn_deprecated_command_line_argument(
    std::string deprecated) {
  std::cerr << "Warning: command line argument '" << deprecated
            << "' is deprecated."
            << " Raised by Kokkos::initialize()." << std::endl;
}

void Kokkos::Impl::warn_deprecated_command_line_argument(
    std::string deprecated, std::string use_instead) {
  std::cerr << "Warning: command line argument '" << deprecated
            << "' is deprecated."
            << " Use '" << use_instead << "' instead."
            << " Raised by Kokkos::initialize()." << std::endl;
}

namespace {
std::vector<std::regex> do_not_warn_regular_expressions{
    std::regex{"--kokkos-tool.*", std::regex::egrep},
};
}

void Kokkos::Impl::do_not_warn_not_recognized_command_line_argument(
    std::regex ignore) {
  do_not_warn_regular_expressions.push_back(std::move(ignore));
}

void Kokkos::Impl::warn_not_recognized_command_line_argument(
    std::string not_recognized) {
  for (auto const& ignore : do_not_warn_regular_expressions) {
    if (std::regex_match(not_recognized, ignore)) {
      return;
    }
  }
  std::cerr << "Warning: command line argument '" << not_recognized
            << "' is not recognized."
            << " Raised by Kokkos::initialize()." << std::endl;
}
