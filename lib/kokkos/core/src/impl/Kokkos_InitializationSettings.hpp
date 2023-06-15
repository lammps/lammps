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

#ifndef KOKKOS_INITIALIZATION_SETTINGS_HPP
#define KOKKOS_INITIALIZATION_SETTINGS_HPP

#include <Kokkos_Macros.hpp>

#include <climits>
#include <string>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
struct InitArguments {
  int num_threads;
  int num_numa;
  int device_id;
  int ndevices;
  int skip_device;
  bool disable_warnings;
  bool tune_internals;
  bool tool_help        = false;
  std::string tool_lib  = {};
  std::string tool_args = {};

  KOKKOS_DEPRECATED_WITH_COMMENT("Use InitializationSettings instead!")
  InitArguments(int nt = -1, int nn = -1, int dv = -1, bool dw = false,
                bool ti = false)
      : num_threads{nt},
        num_numa{nn},
        device_id{dv},
        ndevices{-1},
        skip_device{9999},
        disable_warnings{dw},
        tune_internals{ti} {}
};
#endif

namespace Impl {
// FIXME_CXX17 replace with std::optional
template <class>
struct InitializationSettingsHelper;
template <>
struct InitializationSettingsHelper<int> {
  using value_type   = int;
  using storage_type = int;

  static constexpr storage_type unspecified = INT_MIN;
};
template <>
struct InitializationSettingsHelper<bool> {
  using value_type   = bool;
  using storage_type = char;

  static constexpr storage_type unspecified = CHAR_MAX;
  static_assert(static_cast<storage_type>(true) != unspecified &&
                    static_cast<storage_type>(false) != unspecified,
                "");
};
template <>
struct InitializationSettingsHelper<std::string> {
  using value_type   = std::string;
  using storage_type = std::string;

  // prefer c-string to avoid static initialization order nightmare
  static constexpr char unspecified[] =
      "some string we don't expect user would ever provide";
};
}  // namespace Impl

class InitializationSettings {
#define KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER(NAME) \
  impl_do_not_use_i_really_mean_it_##NAME##_

#define KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE(NAME) impl_##NAME##_type

#define KOKKOS_IMPL_DECLARE(TYPE, NAME)                                      \
 private:                                                                    \
  using KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE(NAME) = TYPE;                 \
  Impl::InitializationSettingsHelper<TYPE>::storage_type                     \
      KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER(NAME) =                              \
          Impl::InitializationSettingsHelper<TYPE>::unspecified;             \
                                                                             \
 public:                                                                     \
  InitializationSettings& set_##NAME(                                        \
      Impl::InitializationSettingsHelper<TYPE>::value_type NAME) {           \
    KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER(NAME) = NAME;                          \
    return *this;                                                            \
  }                                                                          \
  bool has_##NAME() const noexcept {                                         \
    return KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER(NAME) !=                        \
           Impl::InitializationSettingsHelper<                               \
               KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE(NAME)>::unspecified;   \
  }                                                                          \
  KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE(NAME) get_##NAME() const noexcept { \
    return KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER(NAME);                          \
  }                                                                          \
  static_assert(true, "no-op to require trailing semicolon")

 public:
  KOKKOS_IMPL_DECLARE(int, num_threads);
  KOKKOS_IMPL_DECLARE(int, device_id);
  KOKKOS_IMPL_DECLARE(std::string, map_device_id_by);
  KOKKOS_IMPL_DECLARE(int, num_devices);  // deprecated
  KOKKOS_IMPL_DECLARE(int, skip_device);  // deprecated
  KOKKOS_IMPL_DECLARE(bool, disable_warnings);
  KOKKOS_IMPL_DECLARE(bool, print_configuration);
  KOKKOS_IMPL_DECLARE(bool, tune_internals);
  KOKKOS_IMPL_DECLARE(bool, tools_help);
  KOKKOS_IMPL_DECLARE(std::string, tools_libs);
  KOKKOS_IMPL_DECLARE(std::string, tools_args);

#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE
#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER
#undef KOKKOS_IMPL_DECLARE

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
 public:
  InitializationSettings() = default;

  InitializationSettings(InitArguments const& old) {
    if (old.num_threads != -1) {
      set_num_threads(old.num_threads);
    }
    if (old.device_id != -1) {
      set_device_id(old.device_id);
    }
    if (old.ndevices != -1) {
      set_num_devices(old.ndevices);
    }
    if (old.skip_device != 9999) {
      set_skip_device(old.skip_device);
    }
    if (old.disable_warnings) {
      set_disable_warnings(true);
    }
    if (old.tune_internals) {
      set_tune_internals(true);
    }
    if (old.tool_help) {
      set_tools_help(true);
    }
    if (!old.tool_lib.empty()) {
      set_tools_libs(old.tool_lib);
    }
    if (!old.tool_args.empty()) {
      set_tools_args(old.tool_args);
    }
  }
#endif
};

}  // namespace Kokkos

#endif
