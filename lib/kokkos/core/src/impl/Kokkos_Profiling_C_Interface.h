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

#ifndef KOKKOS_PROFILING_C_INTERFACE_HPP
#define KOKKOS_PROFILING_C_INTERFACE_HPP

#ifdef __cplusplus
#include <cstddef>
#include <cstdint>
#else
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#endif

#define KOKKOSP_INTERFACE_VERSION 20200625

// Profiling

struct Kokkos_Profiling_KokkosPDeviceInfo {
  size_t deviceID;
};

struct Kokkos_Profiling_SpaceHandle {
  char name[64];
};

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_initFunction)(
    const int, const uint64_t, const uint32_t,
    struct Kokkos_Profiling_KokkosPDeviceInfo*);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_finalizeFunction)();
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_beginFunction)(const char*, const uint32_t,
                                               uint64_t*);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_endFunction)(uint64_t);

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_pushFunction)(const char*);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_popFunction)();

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_allocateDataFunction)(
    const struct Kokkos_Profiling_SpaceHandle, const char*, const void*,
    const uint64_t);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_deallocateDataFunction)(
    const struct Kokkos_Profiling_SpaceHandle, const char*, const void*,
    const uint64_t);

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_createProfileSectionFunction)(const char*,
                                                              uint32_t*);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_startProfileSectionFunction)(const uint32_t);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_stopProfileSectionFunction)(const uint32_t);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_destroyProfileSectionFunction)(const uint32_t);

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_profileEventFunction)(const char*);

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_beginDeepCopyFunction)(
    struct Kokkos_Profiling_SpaceHandle, const char*, const void*,
    struct Kokkos_Profiling_SpaceHandle, const char*, const void*, uint64_t);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_endDeepCopyFunction)();
typedef void (*Kokkos_Profiling_beginFenceFunction)(const char*, const uint32_t,
                                                    uint64_t*);
typedef void (*Kokkos_Profiling_endFenceFunction)(uint64_t);

// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_dualViewSyncFunction)(const char*,
                                                      const void* const, bool);
// NOLINTNEXTLINE(modernize-use-using): C compatibility
typedef void (*Kokkos_Profiling_dualViewModifyFunction)(const char*,
                                                        const void* const,
                                                        bool);

// Tuning

#define KOKKOS_TOOLS_TUNING_STRING_LENGTH 64
typedef char Kokkos_Tools_Tuning_String[KOKKOS_TOOLS_TUNING_STRING_LENGTH];
union Kokkos_Tools_VariableValue_ValueUnion {
  int64_t int_value;
  double double_value;
  Kokkos_Tools_Tuning_String string_value;
};

union Kokkos_Tools_VariableValue_ValueUnionSet {
  int64_t* int_value;
  double* double_value;
  Kokkos_Tools_Tuning_String* string_value;
};

struct Kokkos_Tools_ValueSet {
  size_t size;
  union Kokkos_Tools_VariableValue_ValueUnionSet values;
};

enum Kokkos_Tools_OptimizationType {
  Kokkos_Tools_Minimize,
  Kokkos_Tools_Maximize
};

struct Kokkos_Tools_OptimzationGoal {
  size_t type_id;
  enum Kokkos_Tools_OptimizationType goal;
};

struct Kokkos_Tools_ValueRange {
  union Kokkos_Tools_VariableValue_ValueUnion lower;
  union Kokkos_Tools_VariableValue_ValueUnion upper;
  union Kokkos_Tools_VariableValue_ValueUnion step;
  bool openLower;
  bool openUpper;
};

enum Kokkos_Tools_VariableInfo_ValueType {
  kokkos_value_double,
  kokkos_value_int64,
  kokkos_value_string,
};

enum Kokkos_Tools_VariableInfo_StatisticalCategory {
  kokkos_value_categorical,  // unordered distinct objects
  kokkos_value_ordinal,      // ordered distinct objects
  kokkos_value_interval,  // ordered distinct objects for which distance matters
  kokkos_value_ratio  // ordered distinct objects for which distance matters,
                      // division matters, and the concept of zero exists
};

enum Kokkos_Tools_VariableInfo_CandidateValueType {
  kokkos_value_set,       // I am one of [2,3,4,5]
  kokkos_value_range,     // I am somewhere in [2,12)
  kokkos_value_unbounded  // I am [text/int/float], but we don't know at
                          // declaration time what values are appropriate. Only
                          // valid for Context Variables
};

union Kokkos_Tools_VariableInfo_SetOrRange {
  struct Kokkos_Tools_ValueSet set;
  struct Kokkos_Tools_ValueRange range;
};

struct Kokkos_Tools_VariableInfo {
  enum Kokkos_Tools_VariableInfo_ValueType type;
  enum Kokkos_Tools_VariableInfo_StatisticalCategory category;
  enum Kokkos_Tools_VariableInfo_CandidateValueType valueQuantity;
  union Kokkos_Tools_VariableInfo_SetOrRange candidates;
  void* toolProvidedInfo;
};

struct Kokkos_Tools_VariableValue {
  size_t type_id;
  union Kokkos_Tools_VariableValue_ValueUnion value;
  struct Kokkos_Tools_VariableInfo* metadata;
};

typedef void (*Kokkos_Tools_outputTypeDeclarationFunction)(
    const char*, const size_t, struct Kokkos_Tools_VariableInfo* info);
typedef void (*Kokkos_Tools_inputTypeDeclarationFunction)(
    const char*, const size_t, struct Kokkos_Tools_VariableInfo* info);

typedef void (*Kokkos_Tools_requestValueFunction)(
    const size_t, const size_t, const struct Kokkos_Tools_VariableValue*,
    const size_t count, struct Kokkos_Tools_VariableValue*);
typedef void (*Kokkos_Tools_contextBeginFunction)(const size_t);
typedef void (*Kokkos_Tools_contextEndFunction)(
    const size_t, struct Kokkos_Tools_VariableValue);
typedef void (*Kokkos_Tools_optimizationGoalDeclarationFunction)(
    const size_t, const struct Kokkos_Tools_OptimzationGoal goal);

typedef void (*function_pointer)();

struct Kokkos_Profiling_EventSet {
  Kokkos_Profiling_initFunction init;
  Kokkos_Profiling_finalizeFunction finalize;
  Kokkos_Profiling_beginFunction begin_parallel_for;
  Kokkos_Profiling_endFunction end_parallel_for;
  Kokkos_Profiling_beginFunction begin_parallel_reduce;
  Kokkos_Profiling_endFunction end_parallel_reduce;
  Kokkos_Profiling_beginFunction begin_parallel_scan;
  Kokkos_Profiling_endFunction end_parallel_scan;
  Kokkos_Profiling_pushFunction push_region;
  Kokkos_Profiling_popFunction pop_region;
  Kokkos_Profiling_allocateDataFunction allocate_data;
  Kokkos_Profiling_deallocateDataFunction deallocate_data;
  Kokkos_Profiling_createProfileSectionFunction create_profile_section;
  Kokkos_Profiling_startProfileSectionFunction start_profile_section;
  Kokkos_Profiling_stopProfileSectionFunction stop_profile_section;
  Kokkos_Profiling_destroyProfileSectionFunction destroy_profile_section;
  Kokkos_Profiling_profileEventFunction profile_event;
  Kokkos_Profiling_beginDeepCopyFunction begin_deep_copy;
  Kokkos_Profiling_endDeepCopyFunction end_deep_copy;
  Kokkos_Profiling_beginFenceFunction begin_fence;
  Kokkos_Profiling_endFenceFunction end_fence;
  Kokkos_Profiling_dualViewSyncFunction sync_dual_view;
  Kokkos_Profiling_dualViewModifyFunction modify_dual_view;
  char profiling_padding[12 * sizeof(function_pointer)];
  Kokkos_Tools_outputTypeDeclarationFunction declare_output_type;
  Kokkos_Tools_inputTypeDeclarationFunction declare_input_type;
  Kokkos_Tools_requestValueFunction request_output_values;
  Kokkos_Tools_contextBeginFunction begin_tuning_context;
  Kokkos_Tools_contextEndFunction end_tuning_context;
  Kokkos_Tools_optimizationGoalDeclarationFunction declare_optimization_goal;
  char padding[234 *
               sizeof(function_pointer)];  // allows us to add another 256
                                           // events to the Tools interface
                                           // without changing struct layout
};

#endif  // KOKKOS_PROFILING_C_INTERFACE_HPP
