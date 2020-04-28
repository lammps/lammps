#include <Kokkos_OpenMPTarget.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

class OpenMPTargetInternal {
 private:
  OpenMPTargetInternal()                            = default;
  OpenMPTargetInternal(const OpenMPTargetInternal&) = default;
  OpenMPTargetInternal& operator=(const OpenMPTargetInternal&) = default;

 public:
  void fence();

  /** \brief  Return the maximum amount of concurrency.  */
  int concurrency();

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream&, const bool detail = false);

  static const char* name();

  //! Free any resources being consumed by the device.
  void impl_finalize();

  //! Has been initialized
  int impl_is_initialized();

  //! Initialize, telling the CUDA run-time library which device to use.
  void impl_initialize();

  static OpenMPTargetInternal* impl_singleton();

 private:
  bool m_is_initialized = false;
};
}  // Namespace Impl
}  // Namespace Experimental
}  // Namespace Kokkos
