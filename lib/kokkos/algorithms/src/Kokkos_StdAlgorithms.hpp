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

#ifndef KOKKOS_STD_ALGORITHMS_HPP
#define KOKKOS_STD_ALGORITHMS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_STD_ALGORITHMS
#endif

/// \file Kokkos_StdAlgorithms.hpp
/// \brief Kokkos counterparts for Standard C++ Library algorithms

#include "std_algorithms/impl/Kokkos_Constraints.hpp"
#include "std_algorithms/impl/Kokkos_RandomAccessIterator.hpp"
#include "std_algorithms/Kokkos_BeginEnd.hpp"

// distance
#include "std_algorithms/Kokkos_Distance.hpp"

// note that we categorize below the headers
// following the std classification.

// modifying ops
#include "std_algorithms/Kokkos_Swap.hpp"
#include "std_algorithms/Kokkos_IterSwap.hpp"

// non-modifying sequence
#include "std_algorithms/Kokkos_AdjacentFind.hpp"
#include "std_algorithms/Kokkos_Count.hpp"
#include "std_algorithms/Kokkos_CountIf.hpp"
#include "std_algorithms/Kokkos_AllOf.hpp"
#include "std_algorithms/Kokkos_AnyOf.hpp"
#include "std_algorithms/Kokkos_NoneOf.hpp"
#include "std_algorithms/Kokkos_Equal.hpp"
#include "std_algorithms/Kokkos_Find.hpp"
#include "std_algorithms/Kokkos_FindIf.hpp"
#include "std_algorithms/Kokkos_FindIfNot.hpp"
#include "std_algorithms/Kokkos_FindEnd.hpp"
#include "std_algorithms/Kokkos_FindFirstOf.hpp"
#include "std_algorithms/Kokkos_ForEach.hpp"
#include "std_algorithms/Kokkos_ForEachN.hpp"
#include "std_algorithms/Kokkos_LexicographicalCompare.hpp"
#include "std_algorithms/Kokkos_Mismatch.hpp"
#include "std_algorithms/Kokkos_Search.hpp"
#include "std_algorithms/Kokkos_SearchN.hpp"

// modifying sequence
#include "std_algorithms/Kokkos_Fill.hpp"
#include "std_algorithms/Kokkos_FillN.hpp"
#include "std_algorithms/Kokkos_Replace.hpp"
#include "std_algorithms/Kokkos_ReplaceIf.hpp"
#include "std_algorithms/Kokkos_ReplaceCopyIf.hpp"
#include "std_algorithms/Kokkos_ReplaceCopy.hpp"
#include "std_algorithms/Kokkos_Copy.hpp"
#include "std_algorithms/Kokkos_CopyN.hpp"
#include "std_algorithms/Kokkos_CopyBackward.hpp"
#include "std_algorithms/Kokkos_CopyIf.hpp"
#include "std_algorithms/Kokkos_Transform.hpp"
#include "std_algorithms/Kokkos_Generate.hpp"
#include "std_algorithms/Kokkos_GenerateN.hpp"
#include "std_algorithms/Kokkos_Reverse.hpp"
#include "std_algorithms/Kokkos_ReverseCopy.hpp"
#include "std_algorithms/Kokkos_Move.hpp"
#include "std_algorithms/Kokkos_MoveBackward.hpp"
#include "std_algorithms/Kokkos_SwapRanges.hpp"
#include "std_algorithms/Kokkos_Unique.hpp"
#include "std_algorithms/Kokkos_UniqueCopy.hpp"
#include "std_algorithms/Kokkos_Rotate.hpp"
#include "std_algorithms/Kokkos_RotateCopy.hpp"
#include "std_algorithms/Kokkos_Remove.hpp"
#include "std_algorithms/Kokkos_RemoveIf.hpp"
#include "std_algorithms/Kokkos_RemoveCopy.hpp"
#include "std_algorithms/Kokkos_RemoveCopyIf.hpp"
#include "std_algorithms/Kokkos_ShiftLeft.hpp"
#include "std_algorithms/Kokkos_ShiftRight.hpp"

// sorting
#include "std_algorithms/Kokkos_IsSortedUntil.hpp"
#include "std_algorithms/Kokkos_IsSorted.hpp"

// min/max element
#include "std_algorithms/Kokkos_MinElement.hpp"
#include "std_algorithms/Kokkos_MaxElement.hpp"
#include "std_algorithms/Kokkos_MinMaxElement.hpp"

// partitioning
#include "std_algorithms/Kokkos_IsPartitioned.hpp"
#include "std_algorithms/Kokkos_PartitionCopy.hpp"
#include "std_algorithms/Kokkos_PartitionPoint.hpp"

// numeric
#include "std_algorithms/Kokkos_AdjacentDifference.hpp"
#include "std_algorithms/Kokkos_Reduce.hpp"
#include "std_algorithms/Kokkos_TransformReduce.hpp"
#include "std_algorithms/Kokkos_ExclusiveScan.hpp"
#include "std_algorithms/Kokkos_TransformExclusiveScan.hpp"
#include "std_algorithms/Kokkos_InclusiveScan.hpp"
#include "std_algorithms/Kokkos_TransformInclusiveScan.hpp"

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_STD_ALGORITHMS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_STD_ALGORITHMS
#endif
#endif
