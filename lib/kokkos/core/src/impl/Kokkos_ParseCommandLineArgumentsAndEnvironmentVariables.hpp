// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
