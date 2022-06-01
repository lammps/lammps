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

#ifndef KOKKOSP_INTERFACE_HPP
#define KOKKOSP_INTERFACE_HPP

#include <cinttypes>
#include <cstddef>
#include <climits>

#include <cstdlib>

// NOTE: in this Kokkos::Profiling block, do not define anything that shouldn't
// exist should Profiling be disabled

namespace Kokkos {
namespace Tools {
namespace Experimental {

constexpr const uint32_t NumReservedDeviceIDs = 1;

enum SpecialSynchronizationCases : int {
  GlobalDeviceSynchronization     = 1,
  DeepCopyResourceSynchronization = 2,
};

enum struct DeviceType {
  Serial,
  OpenMP,
  Cuda,
  HIP,
  OpenMPTarget,
  HPX,
  Threads,
  SYCL,
  Unknown
};

struct ExecutionSpaceIdentifier {
  DeviceType type;
  uint32_t device_id;
  uint32_t instance_id;
};

constexpr const uint32_t num_type_bits     = 8;
constexpr const uint32_t num_device_bits   = 7;
constexpr const uint32_t num_instance_bits = 17;
constexpr const uint32_t num_avail_bits    = sizeof(uint32_t) * CHAR_BIT;

inline DeviceType devicetype_from_uint32t(const uint32_t in) {
  switch (in) {
    case 0: return DeviceType::Serial;
    case 1: return DeviceType::OpenMP;
    case 2: return DeviceType::Cuda;
    case 3: return DeviceType::HIP;
    case 4: return DeviceType::OpenMPTarget;
    case 5: return DeviceType::HPX;
    case 6: return DeviceType::Threads;
    case 7: return DeviceType::SYCL;
    default: return DeviceType::Unknown;  // TODO: error out?
  }
}

inline ExecutionSpaceIdentifier identifier_from_devid(const uint32_t in) {
  constexpr const uint32_t shift = num_avail_bits - num_type_bits;

  return {devicetype_from_uint32t(in >> shift), /*First 8 bits*/
          (~((uint32_t(-1)) << num_device_bits)) &
              (in >> num_instance_bits),                  /*Next 7 bits */
          (~((uint32_t(-1)) << num_instance_bits)) & in}; /*Last 17 bits*/
}

template <typename ExecutionSpace>
struct DeviceTypeTraits;

template <typename ExecutionSpace>
constexpr uint32_t device_id_root() {
  constexpr auto device_id =
      static_cast<uint32_t>(DeviceTypeTraits<ExecutionSpace>::id);
  return (device_id << num_instance_bits);
}
template <typename ExecutionSpace>
inline uint32_t device_id(ExecutionSpace const& space) noexcept {
  return device_id_root<ExecutionSpace>() + space.impl_instance_id();
}
}  // namespace Experimental
}  // namespace Tools
}  // end namespace Kokkos

#if defined(KOKKOS_ENABLE_LIBDL)
// We check at configure time that libdl is available.
#include <dlfcn.h>
#endif

#include <impl/Kokkos_Profiling_DeviceInfo.hpp>
#include <impl/Kokkos_Profiling_C_Interface.h>

namespace Kokkos {
namespace Tools {

using SpaceHandle = Kokkos_Profiling_SpaceHandle;

}  // namespace Tools

namespace Tools {

namespace Experimental {
using EventSet = Kokkos_Profiling_EventSet;
static_assert(sizeof(EventSet) / sizeof(Kokkos_Tools_functionPointer) == 275,
              "sizeof EventSet has changed, this is an error on the part of a "
              "Kokkos developer");
static_assert(sizeof(Kokkos_Tools_ToolSettings) / sizeof(bool) == 256,
              "sizeof EventSet has changed, this is an error on the part of a "
              "Kokkos developer");
static_assert(sizeof(Kokkos_Tools_ToolProgrammingInterface) /
                      sizeof(Kokkos_Tools_functionPointer) ==
                  32,
              "sizeof EventSet has changed, this is an error on the part of a "
              "Kokkos developer");

using toolInvokedFenceFunction = Kokkos_Tools_toolInvokedFenceFunction;
using provideToolProgrammingInterfaceFunction =
    Kokkos_Tools_provideToolProgrammingInterfaceFunction;
using requestToolSettingsFunction = Kokkos_Tools_requestToolSettingsFunction;
using ToolSettings                = Kokkos_Tools_ToolSettings;
using ToolProgrammingInterface    = Kokkos_Tools_ToolProgrammingInterface;
}  // namespace Experimental
using initFunction           = Kokkos_Profiling_initFunction;
using finalizeFunction       = Kokkos_Profiling_finalizeFunction;
using parseArgsFunction      = Kokkos_Profiling_parseArgsFunction;
using printHelpFunction      = Kokkos_Profiling_printHelpFunction;
using beginFunction          = Kokkos_Profiling_beginFunction;
using endFunction            = Kokkos_Profiling_endFunction;
using pushFunction           = Kokkos_Profiling_pushFunction;
using popFunction            = Kokkos_Profiling_popFunction;
using allocateDataFunction   = Kokkos_Profiling_allocateDataFunction;
using deallocateDataFunction = Kokkos_Profiling_deallocateDataFunction;
using createProfileSectionFunction =
    Kokkos_Profiling_createProfileSectionFunction;
using startProfileSectionFunction =
    Kokkos_Profiling_startProfileSectionFunction;
using stopProfileSectionFunction = Kokkos_Profiling_stopProfileSectionFunction;
using destroyProfileSectionFunction =
    Kokkos_Profiling_destroyProfileSectionFunction;
using profileEventFunction    = Kokkos_Profiling_profileEventFunction;
using beginDeepCopyFunction   = Kokkos_Profiling_beginDeepCopyFunction;
using endDeepCopyFunction     = Kokkos_Profiling_endDeepCopyFunction;
using beginFenceFunction      = Kokkos_Profiling_beginFenceFunction;
using endFenceFunction        = Kokkos_Profiling_endFenceFunction;
using dualViewSyncFunction    = Kokkos_Profiling_dualViewSyncFunction;
using dualViewModifyFunction  = Kokkos_Profiling_dualViewModifyFunction;
using declareMetadataFunction = Kokkos_Profiling_declareMetadataFunction;

}  // namespace Tools

}  // namespace Kokkos

// Profiling

namespace Kokkos {

namespace Profiling {

/** The Profiling namespace is being renamed to Tools.
 * This is reexposing the contents of what used to be the Profiling
 * Interface with their original names, to avoid breaking old code
 */

namespace Experimental {

using Kokkos::Tools::Experimental::device_id;
using Kokkos::Tools::Experimental::DeviceType;
using Kokkos::Tools::Experimental::DeviceTypeTraits;

}  // namespace Experimental

using Kokkos::Tools::allocateDataFunction;
using Kokkos::Tools::beginDeepCopyFunction;
using Kokkos::Tools::beginFunction;
using Kokkos::Tools::createProfileSectionFunction;
using Kokkos::Tools::deallocateDataFunction;
using Kokkos::Tools::destroyProfileSectionFunction;
using Kokkos::Tools::endDeepCopyFunction;
using Kokkos::Tools::endFunction;
using Kokkos::Tools::finalizeFunction;
using Kokkos::Tools::initFunction;
using Kokkos::Tools::parseArgsFunction;
using Kokkos::Tools::popFunction;
using Kokkos::Tools::printHelpFunction;
using Kokkos::Tools::profileEventFunction;
using Kokkos::Tools::pushFunction;
using Kokkos::Tools::SpaceHandle;
using Kokkos::Tools::startProfileSectionFunction;
using Kokkos::Tools::stopProfileSectionFunction;

}  // namespace Profiling
}  // namespace Kokkos

// Tuning

namespace Kokkos {
namespace Tools {
namespace Experimental {
using ValueSet            = Kokkos_Tools_ValueSet;
using ValueRange          = Kokkos_Tools_ValueRange;
using StatisticalCategory = Kokkos_Tools_VariableInfo_StatisticalCategory;
using ValueType           = Kokkos_Tools_VariableInfo_ValueType;
using CandidateValueType  = Kokkos_Tools_VariableInfo_CandidateValueType;
using SetOrRange          = Kokkos_Tools_VariableInfo_SetOrRange;
using VariableInfo        = Kokkos_Tools_VariableInfo;
using OptimizationGoal    = Kokkos_Tools_OptimzationGoal;
using TuningString        = Kokkos_Tools_Tuning_String;
using VariableValue       = Kokkos_Tools_VariableValue;

using outputTypeDeclarationFunction =
    Kokkos_Tools_outputTypeDeclarationFunction;
using inputTypeDeclarationFunction = Kokkos_Tools_inputTypeDeclarationFunction;
using requestValueFunction         = Kokkos_Tools_requestValueFunction;
using contextBeginFunction         = Kokkos_Tools_contextBeginFunction;
using contextEndFunction           = Kokkos_Tools_contextEndFunction;
using optimizationGoalDeclarationFunction =
    Kokkos_Tools_optimizationGoalDeclarationFunction;
}  // end namespace Experimental
}  // end namespace Tools

}  // end namespace Kokkos

#endif
