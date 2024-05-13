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

#ifndef KOKKOS_PARSE_COMMAND_LINE_ARGUMENTS_AND_ENVIRONMENT_VARIABLES_HPP
#define KOKKOS_PARSE_COMMAND_LINE_ARGUMENTS_AND_ENVIRONMENT_VARIABLES_HPP

// These declaration are only provided for testing purposes
namespace Kokkos {
class InitializationSettings;
namespace Impl {
void parse_command_line_arguments(int& argc, char* argv[],
                                  InitializationSettings& settings);
void parse_environment_variables(InitializationSettings& settings);
}  // namespace Impl
}  // namespace Kokkos

#endif
