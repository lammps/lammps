// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// as we use min and max as names in our headers and these are defined as macros
// in windows.h, we added a pragma push to disable them at the beginning of
// Kokkos_Core and pop them back into existence at the end.
#include <windows.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
