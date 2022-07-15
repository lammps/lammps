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

#ifndef KOKKOS_SYCL_INSTANCE_HPP_
#define KOKKOS_SYCL_INSTANCE_HPP_

#include <optional>
#include <CL/sycl.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Profiling.hpp>
namespace Kokkos {
namespace Experimental {
namespace Impl {

class SYCLInternal {
 public:
  using size_type = int;

  SYCLInternal() = default;
  ~SYCLInternal();

  SYCLInternal(const SYCLInternal&) = delete;
  SYCLInternal& operator=(const SYCLInternal&) = delete;
  SYCLInternal& operator=(SYCLInternal&&) = delete;
  SYCLInternal(SYCLInternal&&)            = delete;

  void* scratch_space(const std::size_t size);
  void* scratch_flags(const std::size_t size);
  void* resize_team_scratch_space(std::int64_t bytes,
                                  bool force_shrink = false);

  uint32_t impl_get_instance_id() const;
  int m_syclDev = 0;

  size_t m_maxWorkgroupSize   = 0;
  uint32_t m_maxConcurrency   = 0;
  uint64_t m_maxShmemPerBlock = 0;

  std::size_t m_scratchSpaceCount = 0;
  size_type* m_scratchSpace       = nullptr;
  std::size_t m_scratchFlagsCount = 0;
  size_type* m_scratchFlags       = nullptr;
  // mutex to access shared memory
  mutable std::mutex m_mutexScratchSpace;

  int64_t m_team_scratch_current_size = 0;
  void* m_team_scratch_ptr            = nullptr;
  mutable std::mutex m_team_scratch_mutex;

  uint32_t m_instance_id = Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::SYCL>(reinterpret_cast<uintptr_t>(this));
  std::optional<sycl::queue> m_queue;

  // Using std::vector<std::optional<sycl::queue>> reveals a compiler bug when
  // compiling for the CUDA backend. Storing pointers instead works around this.
  static std::vector<std::optional<sycl::queue>*> all_queues;
  // We need a mutex for thread safety when modifying all_queues.
  static std::mutex mutex;

  // USMObjectMem is a reusable buffer for a single object
  // in USM memory
  template <sycl::usm::alloc Kind>
  class USMObjectMem {
   public:
    void reset();

    void reset(sycl::queue q, uint32_t instance_id) {
      m_instance_id = instance_id;
      reset();
      m_q.emplace(std::move(q));
    }
    USMObjectMem() = default;
    explicit USMObjectMem(sycl::queue q, uint32_t instance_id) noexcept
        : m_q(std::move(q)), m_instance_id(instance_id) {}

    USMObjectMem(USMObjectMem const&) = delete;
    USMObjectMem(USMObjectMem&&)      = delete;
    USMObjectMem& operator=(USMObjectMem&&) = delete;
    USMObjectMem& operator=(USMObjectMem const&) = delete;

    ~USMObjectMem() { reset(); };

    void* data() noexcept { return m_data; }
    const void* data() const noexcept { return m_data; }

    size_t capacity() const noexcept { return m_capacity; }

    // reserve() allocates space for at least n bytes
    // returns the new capacity
    size_t reserve(size_t n);

   private:
    using AllocationSpace = std::conditional_t<
        Kind == sycl::usm::alloc::device,
        Kokkos::Experimental::SYCLDeviceUSMSpace,
        std::conditional_t<Kind == sycl::usm::alloc::shared,
                           Kokkos::Experimental::SYCLSharedUSMSpace,
                           Kokkos::Experimental::SYCLHostUSMSpace>>;

   public:
    // Performs either sycl::memcpy (for USM device memory) or std::memcpy
    // (otherwise) and returns a reference to the copied object.
    template <typename T>
    T& copy_from(const T& t) {
      m_mutex.lock();
      fence();
      reserve(sizeof(T));
      if constexpr (sycl::usm::alloc::device == Kind) {
        std::memcpy(static_cast<void*>(m_staging.get()), std::addressof(t),
                    sizeof(T));
        m_copy_event = m_q->memcpy(m_data, m_staging.get(), sizeof(T));
      } else
        std::memcpy(m_data, std::addressof(t), sizeof(T));
      return *reinterpret_cast<T*>(m_data);
    }

    void fence() {
      SYCLInternal::fence(
          m_last_event,
          "Kokkos::Experimental::SYCLInternal::USMObject fence to wait for "
          "last event to finish",
          m_instance_id);
    }

    void register_event(sycl::event event) {
      assert(m_last_event
                 .get_info<sycl::info::event::command_execution_status>() ==
             sycl::info::event_command_status::complete);
      m_last_event = event;
      m_mutex.unlock();
    }

    sycl::event get_copy_event() const { return m_copy_event; }

   private:
    // USMObjectMem class invariants
    // All four expressions below must evaluate to true:
    //
    //  !m_data == (m_capacity == 0)
    //      m_q || !m_data
    //
    //  The above invariants mean that:
    //  if m_data != nullptr then m_capacity != 0 && m_q != nullopt
    //  if m_data == nullptr then m_capacity == 0

    sycl::event m_copy_event;

    std::optional<sycl::queue> m_q;
    void* m_data = nullptr;
    std::unique_ptr<char[]> m_staging;

    size_t m_capacity = 0;
    sycl::event m_last_event;

    uint32_t m_instance_id;

    // mutex to access the underlying memory
    mutable std::mutex m_mutex;
  };

  // An indirect kernel is one where the functor to be executed is explicitly
  // copied to USM memory before being executed, to get around the
  // trivially copyable limitation of SYCL.
  using IndirectKernelMem = USMObjectMem<sycl::usm::alloc::host>;
  IndirectKernelMem& get_indirect_kernel_mem();

  bool was_finalized = false;

  static SYCLInternal& singleton();

  int verify_is_initialized(const char* const label) const;

  void initialize(const sycl::device& d);

  void initialize(const sycl::queue& q);

  int is_initialized() const { return m_queue.has_value(); }

  void finalize();

 private:
  // fence(...) takes any type with a .wait_and_throw() method
  // (sycl::event and sycl::queue)
  template <typename WAT>
  static void fence_helper(WAT& wat, const std::string& name,
                           uint32_t instance_id);

  const static size_t m_usm_pool_size = 4;
  std::vector<IndirectKernelMem> m_indirectKernelMem{m_usm_pool_size};

  size_t m_pool_next{0};

 public:
  static void fence(sycl::queue& q, const std::string& name,
                    uint32_t instance_id) {
    fence_helper(q, name, instance_id);
  }
  static void fence(sycl::event& e, const std::string& name,
                    uint32_t instance_id) {
    fence_helper(e, name, instance_id);
  }
};

// FIXME_SYCL the limit is 2048 bytes for all arguments handed to a kernel,
// assume for now that the rest doesn't need more than 248 bytes.
#if defined(SYCL_DEVICE_COPYABLE) && defined(KOKKOS_ARCH_INTEL_GPU)
template <typename Functor, typename Storage,
          bool ManualCopy = (sizeof(Functor) >= 1800)>
class SYCLFunctionWrapper;
#else
template <typename Functor, typename Storage,
          bool ManualCopy = (sizeof(Functor) >= 1800 ||
                             !std::is_trivially_copyable_v<Functor>)>
class SYCLFunctionWrapper;
#endif

#if defined(SYCL_DEVICE_COPYABLE) && defined(KOKKOS_ARCH_INTEL_GPU)
template <typename Functor, typename Storage>
class SYCLFunctionWrapper<Functor, Storage, false> {
  // We need a union here so that we can avoid calling a constructor for m_f
  // and can controll all the special member functions.
  union TrivialWrapper {
    TrivialWrapper(){};

    TrivialWrapper(const Functor& f) { std::memcpy(&m_f, &f, sizeof(m_f)); }

    TrivialWrapper(const TrivialWrapper& other) {
      std::memcpy(&m_f, &other.m_f, sizeof(m_f));
    }
    TrivialWrapper(TrivialWrapper&& other) {
      std::memcpy(&m_f, &other.m_f, sizeof(m_f));
    }
    TrivialWrapper& operator=(const TrivialWrapper& other) {
      std::memcpy(&m_f, &other.m_f, sizeof(m_f));
      return *this;
    }
    TrivialWrapper& operator=(TrivialWrapper&& other) {
      std::memcpy(&m_f, &other.m_f, sizeof(m_f));
      return *this;
    }
    ~TrivialWrapper(){};

    Functor m_f;
  } m_functor;

 public:
  SYCLFunctionWrapper(const Functor& functor, Storage&) : m_functor(functor) {}

  const Functor& get_functor() const { return m_functor.m_f; }

  sycl::event get_copy_event() const { return {}; }

  static void register_event(sycl::event) {}
};
#else
template <typename Functor, typename Storage>
class SYCLFunctionWrapper<Functor, Storage, false> {
  const Functor m_functor;

 public:
  SYCLFunctionWrapper(const Functor& functor, Storage&) : m_functor(functor) {}

  const Functor& get_functor() const { return m_functor; }

  sycl::event get_copy_event() const { return {}; }

  static void register_event(sycl::event) {}
};
#endif

template <typename Functor, typename Storage>
class SYCLFunctionWrapper<Functor, Storage, true> {
  std::reference_wrapper<const Functor> m_kernelFunctor;
  std::reference_wrapper<Storage> m_storage;

 public:
  SYCLFunctionWrapper(const Functor& functor, Storage& storage)
      : m_kernelFunctor(storage.copy_from(functor)), m_storage(storage) {}

  std::reference_wrapper<const Functor> get_functor() const {
    return m_kernelFunctor;
  }

  sycl::event get_copy_event() const {
    return m_storage.get().get_copy_event();
  }

  void register_event(sycl::event event) {
    m_storage.get().register_event(event);
  }
};

template <typename Functor, typename Storage>
auto make_sycl_function_wrapper(const Functor& functor, Storage& storage) {
  return SYCLFunctionWrapper<Functor, Storage>(functor, storage);
}
}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#if defined(SYCL_DEVICE_COPYABLE) && defined(KOKKOS_ARCH_INTEL_GPU)
template <typename Functor, typename Storage>
struct sycl::is_device_copyable<
    Kokkos::Experimental::Impl::SYCLFunctionWrapper<Functor, Storage, false>>
    : std::true_type {};

template <typename Functor, typename Storage>
struct sycl::is_device_copyable<
    const Kokkos::Experimental::Impl::SYCLFunctionWrapper<Functor, Storage,
                                                          false>>
    : std::true_type {};
#endif
#endif
