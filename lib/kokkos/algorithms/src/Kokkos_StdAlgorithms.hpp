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

/// \file Kokkos_StdAlgorithms.hpp
/// \brief Kokkos counterparts for Standard C++ Library algorithms

#include <std_algorithms/Kokkos_Constraints.hpp>
#include <std_algorithms/Kokkos_RandomAccessIterator.hpp>
#include <std_algorithms/Kokkos_BeginEnd.hpp>

// distance
#include <std_algorithms/Kokkos_Distance.hpp>

// move, swap, iter_swap
#include "std_algorithms/Kokkos_ModifyingOperations.hpp"

// find, find_if, find_if_not
// for_each, for_each_n
// mismatch
// equal
// count_if, count
// all_of, any_of, none_of
// adjacent_find
// lexicographical_compare
// search, search_n
// find_first_of, find_end
#include <std_algorithms/Kokkos_NonModifyingSequenceOperations.hpp>

// replace, replace_copy_if, replace_copy, replace_if
// copy, copy_n, copy_backward, copy_if
// fill, fill_n
// transform
// generate, generate_n
// reverse, reverse_copy
// move, move_backward
// swap_ranges
// unique, unique_copy
// rotate, rotate_copy
// remove, remove_if, remove_copy, remove_copy_if
// shift_left, shift_right
#include <std_algorithms/Kokkos_ModifyingSequenceOperations.hpp>

// is_sorted_until, is_sorted
#include <std_algorithms/Kokkos_SortingOperations.hpp>

// min_element, max_element, minmax_element
#include <std_algorithms/Kokkos_MinMaxElementOperations.hpp>

// is_partitioned, partition_copy, partition_point
#include <std_algorithms/Kokkos_PartitioningOperations.hpp>

// adjacent_difference
// reduce, transform_reduce
// exclusive_scan, transform_exclusive_scan
// inclusive_scan, transform_inclusive_scan
#include <std_algorithms/Kokkos_Numeric.hpp>

#endif
