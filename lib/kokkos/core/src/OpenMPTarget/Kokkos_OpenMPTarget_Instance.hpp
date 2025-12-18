// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMPTARGET_INSTANCE_HPP
#define KOKKOS_OPENMPTARGET_INSTANCE_HPP

namespace Kokkos {
namespace Experimental {
namespace Impl {

enum class openmp_fence_is_static { yes, no };

class OpenMPTargetInternal {
 private:
  OpenMPTargetInternal()                                       = default;
  OpenMPTargetInternal(const OpenMPTargetInternal&)            = delete;
  OpenMPTargetInternal& operator=(const OpenMPTargetInternal&) = delete;

 public:
  void fence(openmp_fence_is_static is_static = openmp_fence_is_static::no);
  void fence(const std::string& name,
             openmp_fence_is_static is_static = openmp_fence_is_static::no);

  /** \brief  Return the maximum amount of concurrency.  */
  int concurrency() const;

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose) const;

  static const char* name();

  //! Free any resources being consumed by the device.
  void impl_finalize();

  uint32_t impl_get_instance_id() const noexcept;
  //! Initialize, telling the CUDA run-time library which device to use.
  void impl_initialize();

  static OpenMPTargetInternal* impl_singleton();

  static void verify_is_process(const char* const);

  void* get_scratch_ptr();
  void clear_scratch();
  void resize_scratch(int64_t team_reduce_bytes, int64_t team_shared_bytes,
                      int64_t thread_local_bytes, int64_t league_size);

  void* m_scratch_ptr = nullptr;
  std::mutex m_mutex_scratch_ptr;
  int64_t m_scratch_size      = 0;
  uint32_t* m_uniquetoken_ptr = nullptr;

 private:
  uint32_t m_instance_id = Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::OpenMPTarget>(reinterpret_cast<uintptr_t>(this));
};
}  // Namespace Impl
}  // Namespace Experimental
}  // Namespace Kokkos

#endif  // KOKKOS_OPENMPTARGET_INSTANCE_HPP
