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

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)

// FIXME_OPENMPTARGET - macro for workaround implementation in UniqueToken
// constructor. undef'ed at the end
#define KOKKOS_IMPL_OPENMPTARGET_WORKAROUND

#include <Kokkos_OpenMPTarget.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_UniqueToken.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Instance.hpp>

#include <sstream>

namespace Kokkos {
namespace Experimental {
namespace Impl {
void OpenMPTargetInternal::fence() {}
int OpenMPTargetInternal::concurrency() { return 128000; }
const char* OpenMPTargetInternal::name() { return "OpenMPTarget"; }
void OpenMPTargetInternal::print_configuration(std::ostream& /*stream*/,
                                               const bool) {
  // FIXME_OPENMPTARGET
  printf("Using OpenMPTarget\n");
}

void OpenMPTargetInternal::impl_finalize() {
  m_is_initialized = false;
  Kokkos::Impl::OpenMPTargetExec space;
  if (space.m_lock_array != nullptr) space.clear_lock_array();

  if (space.m_uniquetoken_ptr != nullptr)
    Kokkos::kokkos_free<Kokkos::Experimental::OpenMPTargetSpace>(
        space.m_uniquetoken_ptr);
}
void OpenMPTargetInternal::impl_initialize() { m_is_initialized = true; }
int OpenMPTargetInternal::impl_is_initialized() {
  return m_is_initialized ? 1 : 0;
}

OpenMPTargetInternal* OpenMPTargetInternal::impl_singleton() {
  static OpenMPTargetInternal self;
  return &self;
}

}  // Namespace Impl

OpenMPTarget::OpenMPTarget()
    : m_space_instance(Impl::OpenMPTargetInternal::impl_singleton()) {}

const char* OpenMPTarget::name() {
  return Impl::OpenMPTargetInternal::impl_singleton()->name();
}
void OpenMPTarget::print_configuration(std::ostream& stream,
                                       const bool detail) {
  m_space_instance->print_configuration(stream, detail);
}

int OpenMPTarget::concurrency() {
  return Impl::OpenMPTargetInternal::impl_singleton()->concurrency();
}
void OpenMPTarget::fence() {
  Impl::OpenMPTargetInternal::impl_singleton()->fence();
}

void OpenMPTarget::impl_initialize() { m_space_instance->impl_initialize(); }
void OpenMPTarget::impl_finalize() { m_space_instance->impl_finalize(); }
int OpenMPTarget::impl_is_initialized() {
  return Impl::OpenMPTargetInternal::impl_singleton()->impl_is_initialized();
}
}  // Namespace Experimental

namespace Impl {
int g_openmptarget_space_factory_initialized =
    Kokkos::Impl::initialize_space_factory<OpenMPTargetSpaceInitializer>(
        "160_OpenMPTarget");

void OpenMPTargetSpaceInitializer::initialize(const InitArguments& args) {
  // Prevent "unused variable" warning for 'args' input struct.  If
  // Serial::initialize() ever needs to take arguments from the input
  // struct, you may remove this line of code.
  (void)args;

  if (std::is_same<Kokkos::Experimental::OpenMPTarget,
                   Kokkos::DefaultExecutionSpace>::value) {
    Kokkos::Experimental::OpenMPTarget().impl_initialize();
    // std::cout << "Kokkos::initialize() fyi: OpenMP enabled and initialized"
    // << std::endl ;
  } else {
    // std::cout << "Kokkos::initialize() fyi: OpenMP enabled but not
    // initialized" << std::endl ;
  }
}

void OpenMPTargetSpaceInitializer::finalize(const bool all_spaces) {
  if (std::is_same<Kokkos::Experimental::OpenMPTarget,
                   Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
    if (Kokkos::Experimental::OpenMPTarget().impl_is_initialized())
      Kokkos::Experimental::OpenMPTarget().impl_finalize();
  }
}

void OpenMPTargetSpaceInitializer::fence() {
  Kokkos::Experimental::OpenMPTarget::fence();
}

void OpenMPTargetSpaceInitializer::print_configuration(std::ostream& msg,
                                                       const bool detail) {
  msg << "OpenMPTarget Execution Space:" << std::endl;
  msg << "  KOKKOS_ENABLE_OPENMPTARGET: ";
  msg << "yes" << std::endl;

  msg << "\nOpenMPTarget Runtime Configuration:" << std::endl;
  Kokkos::Experimental::OpenMPTarget().print_configuration(msg, detail);
}

}  // namespace Impl
}  // Namespace Kokkos

namespace Kokkos {
namespace Experimental {

UniqueToken<Kokkos::Experimental::OpenMPTarget,
            Kokkos::Experimental::UniqueTokenScope::Global>::
    UniqueToken(Kokkos::Experimental::OpenMPTarget const&) {
#ifdef KOKKOS_IMPL_OPENMPTARGET_WORKAROUND
  uint32_t* ptr = Kokkos::Impl::OpenMPTargetExec::m_uniquetoken_ptr;
  int count     = Kokkos::Experimental::OpenMPTarget().concurrency();
  if (ptr == nullptr) {
    int size = count * sizeof(uint32_t);
    ptr      = static_cast<uint32_t*>(
        Kokkos::kokkos_malloc<Kokkos::Experimental::OpenMPTargetSpace>(
            "Kokkos::OpenMPTarget::m_uniquetoken_ptr", size));
    std::vector<uint32_t> h_buf(count, 0);
    OMPT_SAFE_CALL(omp_target_memcpy(ptr, h_buf.data(), size, 0, 0,
                                     omp_get_default_device(),
                                     omp_get_initial_device()));

    Kokkos::Impl::OpenMPTargetExec::m_uniquetoken_ptr = ptr;
  }
#else
// FIXME_OPENMPTARGET - 2 versions of non-working implementations to fill `ptr`
// with 0's
// Version 1 - Creating a target region and filling the
// pointer Error - CUDA error: named symbol not found
#pragma omp target teams distribute parallel for is_device_ptr(ptr) \
    map(to                                                          \
        : size)
  for (int i = 0; i < count; ++i) ptr[i] = 0;

  // Version 2 : Allocating a view on the device and filling it with a scalar
  // value of 0.
  Kokkos::View<uint32_t*, Kokkos::Experimental::OpenMPTargetSpace> ptr_view(
      ptr, count);
  Kokkos::deep_copy(ptr_view, 0);
#endif
  m_buffer = ptr;
  m_count  = count;
}
}  // namespace Experimental
}  // namespace Kokkos

#undef KOKKOS_IMPL_OPENMPTARGET_WORKAROUND
#endif  // defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)
