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

  //! Has been initialized
  int impl_is_initialized();
  uint32_t impl_get_instance_id() const noexcept;
  //! Initialize, telling the CUDA run-time library which device to use.
  void impl_initialize();

  static OpenMPTargetInternal* impl_singleton();

  static void verify_is_process(const char* const);
  static void verify_initialized(const char* const);

  void* get_scratch_ptr();
  void clear_scratch();
  void resize_scratch(int64_t team_reduce_bytes, int64_t team_shared_bytes,
                      int64_t thread_local_bytes, int64_t league_size);

  void* m_scratch_ptr = nullptr;
  std::mutex m_mutex_scratch_ptr;
  int64_t m_scratch_size      = 0;
  uint32_t* m_uniquetoken_ptr = nullptr;

 private:
  bool m_is_initialized  = false;
  uint32_t m_instance_id = Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::OpenMPTarget>(reinterpret_cast<uintptr_t>(this));
};
}  // Namespace Impl
}  // Namespace Experimental
}  // Namespace Kokkos

#endif  // KOKKOS_OPENMPTARGET_INSTANCE_HPP
