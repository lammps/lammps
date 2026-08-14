// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CORE_HPP
#define KOKKOS_CORE_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

//----------------------------------------------------------------------------
// In the case windows.h is included before Kokkos_Core.hpp there might be
// errors due to the potentially defined macros with name "min" and "max" in
// windows.h. These collide with the use of "min" and "max" in names inside
// Kokkos. The macros will be redefined at the end of Kokkos_Core.hpp
#if defined(min)
#pragma push_macro("min")
#undef min
#define KOKKOS_IMPL_PUSH_MACRO_MIN
#endif
#if defined(max)
#pragma push_macro("max")
#undef max
#define KOKKOS_IMPL_PUSH_MACRO_MAX
#endif

//----------------------------------------------------------------------------
// Include the execution space header files for the enabled execution spaces.

#include <Kokkos_Core_fwd.hpp>

#include <KokkosCore_Config_DeclareBackend.hpp>

#include <Kokkos_Half.hpp>
#include <Kokkos_AnonymousSpace.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Clamp.hpp>
#include <Kokkos_MinMax.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_MathematicalSpecialFunctions.hpp>
#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_BitManipulation.hpp>
#include <Kokkos_Swap.hpp>
#include <Kokkos_MemoryPool.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_Tuners.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_CopyViews.hpp>
#include <impl/Kokkos_TeamMDPolicy.hpp>
#include <impl/Kokkos_InitializeFinalize.hpp>
#include <impl/Kokkos_ScopeGuard.hpp>
#include <impl/Kokkos_PartitionSpace.hpp>
#include <impl/Kokkos_CStyleMemoryManagement.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <Kokkos_Crs.hpp>
#include <Kokkos_WorkGraphPolicy.hpp>
// Including this in Kokkos_Parallel_Reduce.hpp led to a circular dependency
// because Kokkos::Sum is used in Kokkos_Combined_Reducer.hpp and the default.
// The real answer is to finally break up Kokkos_Parallel_Reduce.hpp into
// smaller parts...
#include <impl/Kokkos_Combined_Reducer.hpp>
// Yet another workaround to deal with circular dependency issues because the
// implementation of the RAII wrapper is using Kokkos::single.
#include <Kokkos_AcquireUniqueTokenImpl.hpp>

//----------------------------------------------------------------------------
// Redefinition of the macros min and max if we pushed them at entry of
// Kokkos_Core.hpp
#if defined(KOKKOS_IMPL_PUSH_MACRO_MIN)
#pragma pop_macro("min")
#undef KOKKOS_IMPL_PUSH_MACRO_MIN
#endif
#if defined(KOKKOS_IMPL_PUSH_MACRO_MAX)
#pragma pop_macro("max")
#undef KOKKOS_IMPL_PUSH_MACRO_MAX
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#endif
