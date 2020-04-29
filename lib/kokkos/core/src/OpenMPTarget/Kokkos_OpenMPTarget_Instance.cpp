#include <OpenMPTarget/Kokkos_OpenMPTarget_Instance.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {
void OpenMPTargetInternal::fence() {}
int OpenMPTargetInternal::concurrency() { return 128000; }
const char* OpenMPTargetInternal::name() { return "OpenMPTarget"; }
void OpenMPTargetInternal::print_configuration(std::ostream& stream,
                                               const bool) {
  printf("Using OpenMPTarget\n");
}

void OpenMPTargetInternal::impl_finalize() { m_is_initialized = false; }
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
}  // Namespace Kokkos
