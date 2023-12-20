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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include "Kokkos_Macros.hpp"
#include "Kokkos_Stacktrace.hpp"

#ifdef KOKKOS_IMPL_ENABLE_STACKTRACE
// backtrace() function for retrieving the stacktrace
#include <execinfo.h>
#endif
#ifdef KOKKOS_IMPL_ENABLE_CXXABI
#include <cxxabi.h>
#endif  // KOKKOS_ENABLE_CXXABI

#include <exception>
#include <iostream>
#include <tuple>
#include <vector>

namespace Kokkos {
namespace Impl {
#ifndef KOKKOS_IMPL_ENABLE_STACKTRACE
int backtrace(void**, int) { return 0; }
char** backtrace_symbols(void* const*, int) { return nullptr; }
#endif

std::string demangle(const std::string& name) {
#ifndef KOKKOS_IMPL_ENABLE_CXXABI
  return name;
#else
  size_t found_end = name.find_first_of("+)", 0, 2);
  if (found_end == std::string::npos) {
    found_end = name.size();
  }
  size_t found_parenthesis = name.find_first_of("(");
  size_t start             = found_parenthesis + 1;
  if (found_parenthesis == std::string::npos) start = 0;

  std::string s = name.substr(start, found_end - start);

  if (s.length() != 0) {
    int status          = 0;
    char* output_buffer = nullptr;
    size_t length       = s.length();
    char* d = abi::__cxa_demangle(s.c_str(), output_buffer, &length, &status);
    if (d != nullptr) {
      s = d;
      free(d);
    }
  }

  // Special cases for "main" and "start" on Mac
  if (s.length() == 0) {
    if (name == "main" || name == "start") {
      s = name;
    }
  }
  return s;
#endif  // KOKKOS_ENABLE_CXXABI
}

class Stacktrace {
 public:
  Stacktrace()                  = delete;
  Stacktrace(const Stacktrace&) = delete;
  Stacktrace& operator=(const Stacktrace&) = delete;
  Stacktrace(Stacktrace&&)                 = delete;
  Stacktrace& operator=(Stacktrace&&) = delete;
  ~Stacktrace()                       = delete;

  // These are public only to avoid wasting an extra stacktrace line.
  // See save_stacktrace below.
  static constexpr int capacity = 100;
  static void* buffer[capacity];
  static int length;

  static std::vector<std::string> lines() {
    char** symbols = backtrace_symbols(buffer, length);
    if (symbols == nullptr) {
      return {};
    } else {
      std::vector<std::string> trace(length);
      for (int i = 0; i < length; ++i) {
        if (symbols[i] != nullptr) {
          trace[i] = std::string(symbols[i]);
        }
      }
      free(symbols);
      return trace;
    }
  }
};

int Stacktrace::length = 0;
void* Stacktrace::buffer[Stacktrace::capacity];

void save_stacktrace() {
  Stacktrace::length = backtrace(Stacktrace::buffer, Stacktrace::capacity);
}

size_t find_first_non_whitespace(const std::string& s, const size_t start_pos) {
  constexpr size_t num_ws_chars = 3;
  const char ws_chars[]         = "\n\t ";
  return s.find_first_not_of(ws_chars, start_pos, num_ws_chars);
}

size_t find_first_whitespace(const std::string& s, const size_t start_pos) {
  constexpr size_t num_ws_chars = 3;
  const char ws_chars[]         = "\n\t ";
  return s.find_first_of(ws_chars, start_pos, num_ws_chars);
}

template <class Callback>
void for_each_token(const std::string& s, Callback c) {
  size_t cur = find_first_non_whitespace(s, 0);
  while (cur != std::string::npos) {
    const size_t end   = find_first_whitespace(s, cur);
    const bool last    = (end == std::string::npos);
    const size_t count = last ? end : size_t(end - cur);
    c(s.substr(cur, count));
    cur = find_first_non_whitespace(s, end);
  }
}

// Search the whole backtrace, column by column, for "main".
// This tells us what column has the function names.
// While we're doing that, figure out the longest column,
// so we can compute spacing correctly.

struct main_column_info {
  bool found_main;
  size_t main_col;
};

main_column_info find_main_column(const std::vector<std::string>& traceback) {
  bool found_main = false;
  size_t main_col = 0;
  for (auto&& entry : traceback) {
    size_t col_count = 0;
    for_each_token(entry, [&](const std::string& s) {
      const size_t pos = s.find("main");
      if (pos != std::string::npos) {
        found_main = true;
        main_col   = col_count;
      }
      ++col_count;
    });
    if (found_main) {
      break;
    }
  }

  return main_column_info{found_main, main_col};
}

void demangle_and_print_traceback_entry(std::ostream& out,
                                        const std::string& traceback_entry,
                                        const bool found_main,
                                        const size_t main_col) {
  std::vector<std::string> tokens;
  size_t cur_col = 0;

  // Print the address column first
  for_each_token(traceback_entry, [&](const std::string& s) {
    if (!(found_main && cur_col == main_col)) {
      out << s;
    }
    ++cur_col;
  });

  out << " ";

  // Then the function name
  cur_col = 0;
  for_each_token(traceback_entry, [&](const std::string& s) {
    if (found_main && cur_col == main_col) {
      out << demangle(s);
    }
    ++cur_col;
  });
}

void demangle_and_print_traceback(std::ostream& out,
                                  const std::vector<std::string>& traceback) {
  const auto result = find_main_column(traceback);
  for (auto&& entry : traceback) {
    demangle_and_print_traceback_entry(out, entry, result.found_main,
                                       result.main_col);
    out << std::endl;
  }
}

void print_saved_stacktrace(std::ostream& out) {
  auto lines = Stacktrace::lines();
  for (auto&& entry : lines) {
    out << entry << std::endl;
  }
}

void print_demangled_saved_stacktrace(std::ostream& out) {
  demangle_and_print_traceback(out, Stacktrace::lines());
}

std::function<void()> user_terminate_handler_post_ = nullptr;

void kokkos_terminate_handler() {
  using std::cerr;
  using std::endl;

  cerr << "Kokkos observes that std::terminate has been called.  "
          "Here is the last saved stack trace.  Note that this does not "
          "necessarily show what called std::terminate."
       << endl
       << endl;
  print_demangled_saved_stacktrace(std::cerr);

  if (user_terminate_handler_post_ != nullptr) {
    user_terminate_handler_post_();
  } else {
    std::abort();
  }
}

void set_kokkos_terminate_handler(std::function<void()> user_post) {
  user_terminate_handler_post_ = user_post;
  std::set_terminate(kokkos_terminate_handler);
}

}  // namespace Impl
}  // namespace Kokkos
