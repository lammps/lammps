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
