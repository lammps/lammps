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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_HPX_HPP
#define KOKKOS_HPX_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_HPX)

#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_HostSpace.hpp>
#include <cstddef>
#include <iosfwd>

#ifdef KOKKOS_ENABLE_HBWSPACE
#include <Kokkos_HBWSpace.hpp>
#endif

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <impl/Kokkos_ConcurrentBitset.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_TaskQueue.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

#include <hpx/local/barrier.hpp>
#include <hpx/local/condition_variable.hpp>
#include <hpx/local/execution.hpp>
#include <hpx/local/future.hpp>
#include <hpx/local/mutex.hpp>
#include <hpx/local/thread.hpp>

#include <Kokkos_UniqueToken.hpp>

#include <iosfwd>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace Kokkos {
namespace Impl {
class hpx_thread_buffer {
  static constexpr std::size_t m_cache_line_size = 64;

  std::size_t m_num_threads      = 0;
  std::size_t m_size_per_thread  = 0;
  std::size_t m_extra_space      = 0;
  std::size_t m_size_total       = 0;
  std::unique_ptr<char[]> m_data = nullptr;

  static constexpr void pad_to_cache_line(std::size_t &size) {
    size = ((size + m_cache_line_size - 1) / m_cache_line_size) *
           m_cache_line_size;
  }

 public:
  hpx_thread_buffer()                          = default;
  ~hpx_thread_buffer()                         = default;
  hpx_thread_buffer(const hpx_thread_buffer &) = delete;
  hpx_thread_buffer(hpx_thread_buffer &&)      = delete;
  hpx_thread_buffer &operator=(const hpx_thread_buffer &) = delete;
  hpx_thread_buffer &operator=(hpx_thread_buffer) = delete;

  void resize(const std::size_t num_threads, const std::size_t size_per_thread,
              const std::size_t extra_space = 0) noexcept;
  void *get(std::size_t thread_num) const noexcept;
  void *get_extra_space() const noexcept;
};

template <typename T>
struct hpx_range {
  T begin;
  T end;
};

template <typename T>
constexpr T get_num_chunks(const T offset, const T chunk_size, const T max) {
  return (max - offset + chunk_size - 1) / chunk_size;
}

template <typename T>
constexpr hpx_range<T> get_chunk_range(const T i_chunk, const T offset,
                                       const T chunk_size, const T max) {
  const T begin = offset + i_chunk * chunk_size;
  const T end   = (std::min)(begin + chunk_size, max);
  return {begin, end};
}

template <typename Policy>
constexpr bool is_light_weight_policy() {
  constexpr Kokkos::Experimental::WorkItemProperty::HintLightWeight_t
      light_weight = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
  return (typename Policy::work_item_property() & light_weight) == light_weight;
}
}  // namespace Impl

namespace Experimental {
class HPX {
 public:
  static constexpr uint32_t impl_default_instance_id() { return 1; }

 private:
  static bool m_hpx_initialized;
  static std::atomic<uint32_t> m_next_instance_id;

 public:
  enum class instance_mode { default_, independent };

 private:
  static uint32_t m_active_parallel_region_count;
  static hpx::spinlock m_active_parallel_region_count_mutex;
  static hpx::condition_variable_any m_active_parallel_region_count_cond;

  struct instance_data {
    instance_data()  = default;
    ~instance_data() = default;
    instance_data(uint32_t instance_id) : m_instance_id(instance_id) {}
    instance_data(uint32_t instance_id,
                  hpx::execution::experimental::unique_any_sender<> &&sender)
        : m_instance_id(instance_id), m_sender{std::move(sender)} {}

    instance_data(const instance_data &) = delete;
    instance_data(instance_data &&)      = delete;
    instance_data &operator=(const instance_data &) = delete;
    instance_data &operator=(instance_data) = delete;

    uint32_t m_instance_id{HPX::impl_default_instance_id()};
    hpx::execution::experimental::unique_any_sender<> m_sender{
        hpx::execution::experimental::just()};
    Kokkos::Impl::hpx_thread_buffer m_buffer;
    hpx::spinlock m_sender_mutex;
  };

  static void default_instance_deleter(instance_data *) {}
  static instance_data m_default_instance_data;
  Kokkos::Impl::HostSharedPtr<instance_data> m_instance_data;

 public:
  using execution_space      = HPX;
  using memory_space         = HostSpace;
  using device_type          = Kokkos::Device<execution_space, memory_space>;
  using array_layout         = LayoutRight;
  using size_type            = memory_space::size_type;
  using scratch_memory_space = ScratchMemorySpace<HPX>;

  HPX()
      : m_instance_data(Kokkos::Impl::HostSharedPtr<instance_data>(
            &m_default_instance_data, &default_instance_deleter)) {}
  ~HPX() = default;
  HPX(instance_mode mode)
      : m_instance_data(
            mode == instance_mode::independent
                ? (Kokkos::Impl::HostSharedPtr<instance_data>(
                      new instance_data(m_next_instance_id++)))
                : Kokkos::Impl::HostSharedPtr<instance_data>(
                      &m_default_instance_data, &default_instance_deleter)) {}
  HPX(hpx::execution::experimental::unique_any_sender<> &&sender)
      : m_instance_data(Kokkos::Impl::HostSharedPtr<instance_data>(
            new instance_data(m_next_instance_id++, std::move(sender)))) {}

  HPX(HPX &&other)      = default;
  HPX(const HPX &other) = default;

  HPX &operator=(HPX &&) = default;
  HPX &operator=(const HPX &) = default;

  void print_configuration(std::ostream &os, bool /*verbose*/ = false) const;
  instance_data &impl_get_instance_data() const noexcept {
    KOKKOS_EXPECTS(m_instance_data.get());
    return *m_instance_data.get();
  }
  uint32_t impl_instance_id() const noexcept {
    return impl_get_instance_data().m_instance_id;
  }

  static bool &impl_get_in_parallel() noexcept;

  struct impl_in_parallel_scope {
    impl_in_parallel_scope() noexcept;
    ~impl_in_parallel_scope() noexcept;
    impl_in_parallel_scope(impl_in_parallel_scope &&)      = delete;
    impl_in_parallel_scope(impl_in_parallel_scope const &) = delete;
    impl_in_parallel_scope &operator=(impl_in_parallel_scope &&) = delete;
    impl_in_parallel_scope &operator=(impl_in_parallel_scope const &) = delete;
  };

  struct impl_not_in_parallel_scope {
    impl_not_in_parallel_scope() noexcept;
    ~impl_not_in_parallel_scope() noexcept;
    impl_not_in_parallel_scope(impl_not_in_parallel_scope &&)      = delete;
    impl_not_in_parallel_scope(impl_not_in_parallel_scope const &) = delete;
    impl_not_in_parallel_scope &operator=(impl_not_in_parallel_scope &&) =
        delete;
    impl_not_in_parallel_scope &operator=(impl_not_in_parallel_scope const &) =
        delete;
  };

  static bool in_parallel(HPX const & = HPX()) noexcept {
    return impl_get_in_parallel();
  }

  static void impl_decrement_active_parallel_region_count();
  static void impl_increment_active_parallel_region_count();

  void impl_instance_fence_locked(const std::string &name) const;
  void impl_instance_fence(const std::string &name) const;
  static void impl_static_fence(const std::string &name);

  void fence(
      const std::string &name =
          "Kokkos::Experimental::HPX::fence: Unnamed Instance Fence") const {
    impl_instance_fence(name);
  }

  static bool is_asynchronous(HPX const & = HPX()) noexcept {
#if defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH)
    return true;
#else
    return false;
#endif
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  template <typename F>
  KOKKOS_DEPRECATED static void partition_master(
      F const &, int requested_num_partitions = 0, int = 0) {
    if (requested_num_partitions > 1) {
      Kokkos::abort(
          "Kokkos::Experimental::HPX::partition_master: can't partition an "
          "HPX instance\n");
    }
  }
#endif

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency();
#else
  int concurrency() const;
#endif
  static void impl_initialize(InitializationSettings const &);
  static bool impl_is_initialized() noexcept;
  static void impl_finalize();
  static int impl_thread_pool_size() noexcept;
  static int impl_thread_pool_rank() noexcept;
  static int impl_thread_pool_size(int depth);

  static int impl_max_hardware_threads() noexcept {
    return hpx::threads::hardware_concurrency();
  }

  static int impl_hardware_thread_id() noexcept {
    return hpx::get_worker_thread_num();
  }

  Kokkos::Impl::hpx_thread_buffer &impl_get_buffer() const noexcept {
    return impl_get_instance_data().m_buffer;
  }

  hpx::execution::experimental::unique_any_sender<> &impl_get_sender() const
      noexcept {
    return impl_get_instance_data().m_sender;
  }

  hpx::execution::experimental::any_sender<> get_sender() const noexcept {
    std::lock_guard l(impl_get_sender_mutex());
    auto &s      = impl_get_sender();
    auto split_s = hpx::execution::experimental::split(std::move(s));
    s            = split_s;
    return hpx::execution::experimental::any_sender<>{split_s};
  }

  hpx::future<void> impl_get_future() const noexcept {
    return hpx::execution::experimental::make_future(get_sender());
  }

  hpx::spinlock &impl_get_sender_mutex() const noexcept {
    return impl_get_instance_data().m_sender_mutex;
  }

  template <typename I>
  void impl_bulk_plain_erased(
      [[maybe_unused]] bool force_synchronous, bool is_light_weight_policy,
      std::function<void(I)> &&f, I const n,
      hpx::threads::thread_stacksize stacksize =
          hpx::threads::thread_stacksize::default_) const {
    Kokkos::Experimental::HPX::impl_increment_active_parallel_region_count();

    namespace ex = hpx::execution::experimental;

    auto &sen = impl_get_sender();
    auto &mut = impl_get_sender_mutex();

    std::lock_guard<hpx::spinlock> l(mut);
    hpx::util::ignore_lock(&mut);

    {
      if (n == 1 && is_light_weight_policy &&
          (hpx::threads::get_self_ptr() != nullptr)) {
        sen = std::move(sen) | ex::then(hpx::bind_front(std::move(f), 0)) |
              ex::then(Kokkos::Experimental::HPX::
                           impl_decrement_active_parallel_region_count) |
              ex::ensure_started();
      } else {
        sen = std::move(sen) |
              ex::transfer(
                  ex::with_stacksize(ex::thread_pool_scheduler{}, stacksize)) |
              ex::bulk(n, std::move(f)) |
              ex::then(Kokkos::Experimental::HPX::
                           impl_decrement_active_parallel_region_count) |
              ex::ensure_started();
      }
    }

#if defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH)
    if (force_synchronous)
#endif
    {
      impl_instance_fence_locked(
          "Kokkos::Experimental::HPX: fence due to forced synchronizations");
    }
  }

  template <typename Functor, typename Index>
  void impl_bulk_plain(bool force_synchronous, bool is_light_weight_policy,
                       Functor const &functor, Index const n,
                       hpx::threads::thread_stacksize stacksize =
                           hpx::threads::thread_stacksize::default_) const {
    impl_bulk_plain_erased(force_synchronous, is_light_weight_policy,
                           {[functor](Index i) {
                             impl_in_parallel_scope p;
                             functor.execute_range(i);
                           }},
                           n, stacksize);
  }

  template <typename Index>
  void impl_bulk_setup_finalize_erased(
      [[maybe_unused]] bool force_synchronous, bool is_light_weight_policy,
      std::function<void(Index)> &&f, std::function<void()> &&f_setup,
      std::function<void()> &&f_finalize, Index const n,
      hpx::threads::thread_stacksize stacksize =
          hpx::threads::thread_stacksize::default_) const {
    Kokkos::Experimental::HPX::impl_increment_active_parallel_region_count();

    namespace ex = hpx::execution::experimental;
    using hpx::threads::thread_stacksize;

    auto &sen = impl_get_sender();
    auto &mut = impl_get_sender_mutex();

    std::lock_guard<hpx::spinlock> l(mut);
    hpx::util::ignore_lock(&mut);

    {
      if (n == 1 && is_light_weight_policy &&
          (hpx::threads::get_self_ptr() != nullptr)) {
        sen = std::move(sen) | ex::then(std::move(f_setup)) |
              ex::then(hpx::bind_front(std::move(f), 0)) |
              ex::then(std::move(f_finalize)) |
              ex::then(Kokkos::Experimental::HPX::
                           impl_decrement_active_parallel_region_count) |
              ex::ensure_started();
      } else {
        sen = std::move(sen) |
              ex::transfer(
                  ex::with_stacksize(ex::thread_pool_scheduler{}, stacksize)) |
              ex::then(std::move(f_setup)) | ex::bulk(n, std::move(f)) |
              ex::then(std::move(f_finalize)) |
              ex::then(Kokkos::Experimental::HPX::
                           impl_decrement_active_parallel_region_count) |
              ex::ensure_started();
      }
    }

#if defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH)
    if (force_synchronous)
#endif
    {
      impl_instance_fence_locked(
          "Kokkos::Experimental::HPX: fence due to forced syncronizations");
    }
  }

  template <typename Functor, typename Index>
  void impl_bulk_setup_finalize(
      bool force_synchronous, bool is_light_weight_policy,
      Functor const &functor, Index const n,
      hpx::threads::thread_stacksize stacksize =
          hpx::threads::thread_stacksize::default_) const {
    impl_bulk_setup_finalize_erased(force_synchronous, is_light_weight_policy,
                                    {[functor](Index i) {
                                      impl_in_parallel_scope p;
                                      functor.execute_range(i);
                                    }},
                                    {[functor]() {
                                      impl_in_parallel_scope p;
                                      functor.setup();
                                    }},
                                    {[functor]() {
                                      impl_in_parallel_scope p;
                                      functor.finalize();
                                    }},
                                    n, stacksize);
  }

  static constexpr const char *name() noexcept { return "HPX"; }

 private:
  friend bool operator==(HPX const &lhs, HPX const &rhs) {
    return lhs.impl_instance_id() == rhs.impl_instance_id();
  }
  friend bool operator!=(HPX const &lhs, HPX const &rhs) {
    return !(lhs == rhs);
  }
};

extern template void HPX::impl_bulk_plain_erased<int>(
    bool, bool, std::function<void(int)> &&, int const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_plain_erased<unsigned int>(
    bool, bool, std::function<void(unsigned int)> &&, unsigned int const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_plain_erased<long>(
    bool, bool, std::function<void(long)> &&, long const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_plain_erased<std::size_t>(
    bool, bool, std::function<void(std::size_t)> &&, std::size_t const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_setup_finalize_erased<int>(
    bool, bool, std::function<void(int)> &&, std::function<void()> &&,
    std::function<void()> &&, int const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_setup_finalize_erased<unsigned int>(
    bool, bool, std::function<void(unsigned int)> &&, std::function<void()> &&,
    std::function<void()> &&, unsigned int const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_setup_finalize_erased<long>(
    bool, bool, std::function<void(long)> &&, std::function<void()> &&,
    std::function<void()> &&, long const,
    hpx::threads::thread_stacksize stacksize) const;

extern template void HPX::impl_bulk_setup_finalize_erased<std::size_t>(
    bool, bool, std::function<void(std::size_t)> &&, std::function<void()> &&,
    std::function<void()> &&, std::size_t const,
    hpx::threads::thread_stacksize stacksize) const;
}  // namespace Experimental

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::Experimental::HPX> {
  static constexpr DeviceType id = DeviceType::HPX;
  static int device_id(const Kokkos::Experimental::HPX &) { return 0; }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {
template <>
struct MemorySpaceAccess<Kokkos::Experimental::HPX::memory_space,
                         Kokkos::Experimental::HPX::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
template <>
class UniqueToken<HPX, UniqueTokenScope::Instance> {
 private:
  using buffer_type = Kokkos::View<uint32_t *, Kokkos::HostSpace>;
  int m_count;
  buffer_type m_buffer_view;
  uint32_t volatile *m_buffer;

 public:
  using execution_space = HPX;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const & = execution_space()) noexcept
      : m_count(execution_space::impl_max_hardware_threads()),
        m_buffer_view(buffer_type()),
        m_buffer(nullptr) {}

  UniqueToken(size_type max_size, execution_space const & = execution_space())
      : m_count(max_size > execution_space::impl_max_hardware_threads()
                    ? execution_space::impl_max_hardware_threads()
                    : max_size),
        m_buffer_view(
            max_size > execution_space::impl_max_hardware_threads()
                ? buffer_type()
                : buffer_type("UniqueToken::m_buffer_view",
                              ::Kokkos::Impl::concurrent_bitset::buffer_bound(
                                  m_count))),
        m_buffer(m_buffer_view.data()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return m_count; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((
        if (m_buffer == nullptr) {
          return execution_space::impl_hardware_thread_id();
        } else {
          const ::Kokkos::pair<int, int> result =
              ::Kokkos::Impl::concurrent_bitset::acquire_bounded(
                  m_buffer, m_count, ::Kokkos::Impl::clock_tic() % m_count);

          if (result.first < 0) {
            ::Kokkos::abort(
                "UniqueToken<HPX> failure to acquire tokens, no tokens "
                "available");
          }
          return result.first;
        }))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int i) const noexcept {
    KOKKOS_IF_ON_HOST((if (m_buffer != nullptr) {
      ::Kokkos::Impl::concurrent_bitset::release(m_buffer, i);
    }))

    KOKKOS_IF_ON_DEVICE(((void)i;))
  }
};

template <>
class UniqueToken<HPX, UniqueTokenScope::Global> {
 public:
  using execution_space = HPX;
  using size_type       = int;
  UniqueToken(execution_space const & = execution_space()) noexcept {}

  // NOTE: Currently this assumes that there is no oversubscription.
  // hpx::get_num_worker_threads can't be used directly because it may yield
  // it's task (problematic if called after hpx::get_worker_thread_num).
  int size() const noexcept { return HPX::impl_max_hardware_threads(); }
  int acquire() const noexcept { return HPX::impl_hardware_thread_id(); }
  void release(int) const noexcept {}
};
}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

struct HPXTeamMember {
 public:
  using execution_space = Kokkos::Experimental::HPX;
  using scratch_memory_space =
      Kokkos::ScratchMemorySpace<Kokkos::Experimental::HPX>;
  using team_handle = HPXTeamMember;

 private:
  scratch_memory_space m_team_shared;

  int m_league_size;
  int m_league_rank;
  int m_team_size;
  int m_team_rank;

 public:
  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space &team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space &team_scratch(const int) const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space &thread_scratch(const int) const {
    return m_team_shared.set_team_thread_mode(0, team_size(), team_rank());
  }

  KOKKOS_INLINE_FUNCTION int league_rank() const noexcept {
    return m_league_rank;
  }

  KOKKOS_INLINE_FUNCTION int league_size() const noexcept {
    return m_league_size;
  }

  KOKKOS_INLINE_FUNCTION int team_rank() const noexcept { return m_team_rank; }
  KOKKOS_INLINE_FUNCTION int team_size() const noexcept { return m_team_size; }

  template <class... Properties>
  constexpr KOKKOS_INLINE_FUNCTION HPXTeamMember(
      const TeamPolicyInternal<Kokkos::Experimental::HPX, Properties...>
          &policy,
      const int team_rank, const int league_rank, void *scratch,
      size_t scratch_size) noexcept
      : m_team_shared(scratch, scratch_size, scratch, scratch_size),
        m_league_size(policy.league_size()),
        m_league_rank(league_rank),
        m_team_size(policy.team_size()),
        m_team_rank(team_rank) {}

  KOKKOS_INLINE_FUNCTION
  void team_barrier() const {}

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType &, const int &) const {}

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(const Closure &closure,
                                             ValueType &value,
                                             const int &) const {
    closure(value);
  }

  template <class ValueType, class JoinOp>
  KOKKOS_INLINE_FUNCTION ValueType team_reduce(const ValueType &value,
                                               const JoinOp &) const {
    return value;
  }

  template <class ReducerType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<is_reducer<ReducerType>::value>
  team_reduce(const ReducerType &) const {}

  template <typename Type>
  KOKKOS_INLINE_FUNCTION Type
  team_scan(const Type &value, Type *const global_accum = nullptr) const {
    if (global_accum) {
      Kokkos::atomic_fetch_add(global_accum, value);
    }

    return 0;
  }
};

template <class... Properties>
class TeamPolicyInternal<Kokkos::Experimental::HPX, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  using traits = PolicyTraits<Properties...>;

  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using member_type = HPXTeamMember;

  //! Execution space of this execution policy:
  using execution_space = Kokkos::Experimental::HPX;

 private:
  typename traits::execution_space m_space{};
  int m_league_size;
  int m_team_size;
  std::size_t m_team_scratch_size[2];
  std::size_t m_thread_scratch_size[2];
  int m_chunk_size;

 public:
  // NOTE: Max size is 1 for simplicity. In most cases more than 1 is not
  // necessary on CPU. Implement later if there is a need.
  template <class FunctorType>
  inline static int team_size_max(const FunctorType &) {
    return 1;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType &) {
    return 1;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType &, const int &) {
    return 1;
  }

  template <class FunctorType>
  int team_size_max(const FunctorType &, const ParallelForTag &) const {
    return 1;
  }

  template <class FunctorType>
  int team_size_max(const FunctorType &, const ParallelReduceTag &) const {
    return 1;
  }

  template <class FunctorType, class ReducerType>
  int team_size_max(const FunctorType &, const ReducerType &,
                    const ParallelReduceTag &) const {
    return 1;
  }

  template <class FunctorType>
  int team_size_recommended(const FunctorType &, const ParallelForTag &) const {
    return 1;
  }

  template <class FunctorType>
  int team_size_recommended(const FunctorType &,
                            const ParallelReduceTag &) const {
    return 1;
  }

  template <class FunctorType, class ReducerType>
  int team_size_recommended(const FunctorType &, const ReducerType &,
                            const ParallelReduceTag &) const {
    return 1;
  }

  static int vector_length_max() { return 1; }

  inline int impl_vector_length() noexcept { return 1; }
  inline bool impl_auto_team_size() noexcept { return false; }
  inline bool impl_auto_vector_length() noexcept { return false; }
  inline void impl_set_vector_length(int) noexcept {}
  inline void impl_set_team_size(int) noexcept {}

 private:
  inline void init(const int league_size_request, const int team_size_request) {
    m_league_size           = league_size_request;
    const int max_team_size = 1;  // TODO: Can't use team_size_max(...) because
                                  // it requires a functor as argument.
    m_team_size =
        team_size_request > max_team_size ? max_team_size : team_size_request;

    if (m_chunk_size > 0) {
      if (!Impl::is_integral_power_of_two(m_chunk_size))
        Kokkos::abort("TeamPolicy blocking granularity must be power of two");
    } else {
      int new_chunk_size = 1;
      while (new_chunk_size * 4 * m_space.concurrency() < m_league_size) {
        new_chunk_size *= 2;
      }

      if (new_chunk_size < 128) {
        new_chunk_size = 1;
        while ((new_chunk_size * m_space.concurrency() < m_league_size) &&
               (new_chunk_size < 128))
          new_chunk_size *= 2;
      }

      m_chunk_size = new_chunk_size;
    }
  }

 public:
  inline int team_size() const { return m_team_size; }
  inline int league_size() const { return m_league_size; }

  size_t scratch_size(const int &level, int team_size_ = -1) const {
    if (team_size_ < 0) {
      team_size_ = m_team_size;
    }
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  inline static int scratch_size_max(int level) {
    return (level == 0 ? 1024 * 32 :  // Roughly L1 size
                20 * 1024 * 1024);    // Limit to keep compatibility with CUDA
  }

 public:
  template <class ExecSpace, class... OtherProperties>
  friend class TeamPolicyInternal;

  const typename traits::execution_space &space() const { return m_space; }

  template <class... OtherProperties>
  TeamPolicyInternal(const TeamPolicyInternal<Kokkos::Experimental::HPX,
                                              OtherProperties...> &p) {
    m_space                  = p.m_space;
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
  }

  TeamPolicyInternal(const typename traits::execution_space &space,
                     int league_size_request, int team_size_request,
                     int /* vector_length_request */ = 1)
      : m_space{space},
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, team_size_request);
  }

  TeamPolicyInternal(const typename traits::execution_space &space,
                     int league_size_request, const Kokkos::AUTO_t &,
                     int /* vector_length_request */ = 1)
      : m_space{space},
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, 1);
  }

  TeamPolicyInternal(const typename traits::execution_space &space,
                     int league_size_request,
                     const Kokkos::AUTO_t &, /* team_size_request */
                     const Kokkos::AUTO_t & /* vector_length_request */)
      : m_space{space},
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, 1);
  }

  TeamPolicyInternal(const typename traits::execution_space &space,
                     int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t & /* vector_length_request */
                     )
      : m_space{space},
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, team_size_request);
  }

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t &, /* team_size_request */
                     const Kokkos::AUTO_t & /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, 1);
  }

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t & /* vector_length_request */
                     )
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, team_size_request);
  }

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int /* vector_length_request */ = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, team_size_request);
  }

  TeamPolicyInternal(int league_size_request, const Kokkos::AUTO_t &,
                     int /* vector_length_request */ = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0) {
    init(league_size_request, 1);
  }

  inline int chunk_size() const { return m_chunk_size; }

  inline TeamPolicyInternal &set_chunk_size(
      typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  inline TeamPolicyInternal &set_scratch_size(const int &level,
                                              const PerTeamValue &per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  inline TeamPolicyInternal &set_scratch_size(
      const int &level, const PerThreadValue &per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  inline TeamPolicyInternal &set_scratch_size(
      const int &level, const PerTeamValue &per_team,
      const PerThreadValue &per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Experimental::HPX> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  void execute_range(const Member i_chunk) const {
    const auto r = get_chunk_range(i_chunk, m_policy.begin(),
                                   m_policy.chunk_size(), m_policy.end());
    for (Member i = r.begin; i < r.end; ++i) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor(i);
      } else {
        m_functor(WorkTag{}, i);
      }
    }
  }

  void execute() const {
    const Member num_chunks =
        get_num_chunks(m_policy.begin(), m_policy.chunk_size(), m_policy.end());
    m_policy.space().impl_bulk_plain(false, is_light_weight_policy<Policy>(),
                                     *this, num_chunks,
                                     hpx::threads::thread_stacksize::nostack);
  }

  inline ParallelFor(const FunctorType &arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Experimental::HPX> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using WorkTag       = typename MDRangePolicy::work_tag;
  using Member        = typename Policy::member_type;
  using iterate_type =
      typename Kokkos::Impl::HostIterateTile<MDRangePolicy, FunctorType,
                                             WorkTag, void>;

  const iterate_type m_iter;
  const Policy m_policy;

 public:
  void execute_range(const Member i_chunk) const {
    const auto r = get_chunk_range(i_chunk, m_policy.begin(),
                                   m_policy.chunk_size(), m_policy.end());
    for (Member i = r.begin; i < r.end; ++i) {
      m_iter(i);
    }
  }

  void execute() const {
    const Member num_chunks =
        get_num_chunks(m_policy.begin(), m_policy.chunk_size(), m_policy.end());
    m_iter.m_rp.space().impl_bulk_plain(
        false, is_light_weight_policy<MDRangePolicy>(), *this, num_chunks,
        hpx::threads::thread_stacksize::nostack);
  }

  inline ParallelFor(const FunctorType &arg_functor, MDRangePolicy arg_policy)
      : m_iter(arg_policy, arg_functor),
        m_policy(Policy(0, arg_policy.m_num_tiles).set_chunk_size(1)) {}
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy &, const Functor &) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {
template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::Experimental::HPX> {
 private:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  using value_type     = typename ReducerType::value_type;
  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_force_synchronous;

 public:
  void setup() const {
    const ReducerType &reducer   = m_functor_reducer.get_reducer();
    const std::size_t value_size = reducer.value_size();
    const int num_worker_threads = m_policy.space().concurrency();

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    buffer.resize(num_worker_threads, value_size);

    for (int t = 0; t < num_worker_threads; ++t) {
      reducer.init(reinterpret_cast<pointer_type>(buffer.get(t)));
    }
  }

  void execute_range(const Member i_chunk) const {
    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    reference_type update =
        ReducerType::reference(reinterpret_cast<pointer_type>(
            buffer.get(Kokkos::Experimental::HPX::impl_hardware_thread_id())));
    const auto r = get_chunk_range(i_chunk, m_policy.begin(),
                                   m_policy.chunk_size(), m_policy.end());
    for (Member i = r.begin; i < r.end; ++i) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor_reducer.get_functor()(i, update);
      } else {
        m_functor_reducer.get_functor()(WorkTag{}, i, update);
      }
    }
  }

  void finalize() const {
    hpx_thread_buffer &buffer    = m_policy.space().impl_get_buffer();
    const ReducerType &reducer   = m_functor_reducer.get_reducer();
    const int num_worker_threads = m_policy.space().concurrency();
    for (int i = 1; i < num_worker_threads; ++i) {
      reducer.join(reinterpret_cast<pointer_type>(buffer.get(0)),
                   reinterpret_cast<pointer_type>(buffer.get(i)));
    }

    pointer_type final_value_ptr =
        reinterpret_cast<pointer_type>(buffer.get(0));

    reducer.final(final_value_ptr);

    if (m_result_ptr != nullptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = final_value_ptr[j];
      }
    }
  }

  void execute() const {
    if (m_policy.end() <= m_policy.begin()) {
      if (m_result_ptr) {
        const ReducerType &reducer = m_functor_reducer.get_reducer();
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
      return;
    }

    const Member num_chunks =
        get_num_chunks(m_policy.begin(), m_policy.chunk_size(), m_policy.end());
    m_policy.space().impl_bulk_setup_finalize(
        m_force_synchronous, is_light_weight_policy<Policy>(), *this,
        num_chunks, hpx::threads::thread_stacksize::nostack);
  }

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType &arg_functor_reducer,
                        Policy arg_policy, const ViewType &arg_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_view.data()),
        m_force_synchronous(!arg_view.impl_track().has_record()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "HPX reduce result must be a View accessible from HostSpace");
  }
};

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>,
                     Kokkos::Experimental::HPX> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;

  using Policy  = typename MDRangePolicy::impl_range_policy;
  using WorkTag = typename MDRangePolicy::work_tag;
  using Member  = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;
  using iterate_type   = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, CombinedFunctorReducerType, WorkTag, reference_type>;

  const iterate_type m_iter;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_force_synchronous;

 public:
  void setup() const {
    const ReducerType &reducer   = m_iter.m_func.get_reducer();
    const std::size_t value_size = reducer.value_size();
    const int num_worker_threads = m_policy.space().concurrency();

    hpx_thread_buffer &buffer = m_iter.m_rp.space().impl_get_buffer();
    buffer.resize(num_worker_threads, value_size);

    for (int t = 0; t < num_worker_threads; ++t) {
      reducer.init(reinterpret_cast<pointer_type>(buffer.get(t)));
    }
  }

  void execute_range(const Member i_chunk) const {
    hpx_thread_buffer &buffer = m_iter.m_rp.space().impl_get_buffer();
    reference_type update =
        ReducerType::reference(reinterpret_cast<pointer_type>(
            buffer.get(Kokkos::Experimental::HPX::impl_hardware_thread_id())));
    const auto r = get_chunk_range(i_chunk, m_policy.begin(),
                                   m_policy.chunk_size(), m_policy.end());
    for (Member i = r.begin; i < r.end; ++i) {
      m_iter(i, update);
    }
  }

  void finalize() const {
    hpx_thread_buffer &buffer    = m_iter.m_rp.space().impl_get_buffer();
    ReducerType reducer          = m_iter.m_func.get_reducer();
    const int num_worker_threads = m_policy.space().concurrency();
    for (int i = 1; i < num_worker_threads; ++i) {
      reducer.join(reinterpret_cast<pointer_type>(buffer.get(0)),
                   reinterpret_cast<pointer_type>(buffer.get(i)));
    }

    pointer_type final_value_ptr =
        reinterpret_cast<pointer_type>(buffer.get(0));

    reducer.final(final_value_ptr);

    if (m_result_ptr != nullptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = final_value_ptr[j];
      }
    }
  }

  void execute() const {
    const Member num_chunks =
        get_num_chunks(m_policy.begin(), m_policy.chunk_size(), m_policy.end());
    m_iter.m_rp.space().impl_bulk_setup_finalize(
        m_force_synchronous, is_light_weight_policy<MDRangePolicy>(), *this,
        num_chunks, hpx::threads::thread_stacksize::nostack);
  }

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType &arg_functor_reducer,
                        MDRangePolicy arg_policy, const ViewType &arg_view)
      : m_iter(arg_policy, arg_functor_reducer),
        m_policy(Policy(0, arg_policy.m_num_tiles).set_chunk_size(1)),
        m_result_ptr(arg_view.data()),
        m_force_synchronous(!arg_view.impl_track().has_record()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "HPX reduce result must be a View accessible from HostSpace");
  }

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy &, const Functor &) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::HPX> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;
  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::SCAN, Policy, FunctorType, void>;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;
  using value_type     = typename Analysis::value_type;
  using barrier_type   = hpx::barrier<>;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  void setup() const {
    const int num_worker_threads = m_policy.space().concurrency();
    const std::size_t value_size = Analysis::value_size(m_functor);

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    buffer.resize(num_worker_threads, 2 * value_size, sizeof(barrier_type));

    new (buffer.get_extra_space()) barrier_type(num_worker_threads);
  }

  void execute_chunk(const Member i_begin, const Member i_end,
                     reference_type update, const bool final) const {
    for (Member i = i_begin; i < i_end; ++i) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor(i, update, final);
      } else {
        m_functor(WorkTag{}, i, update, final);
      }
    }
  }

  void execute_range(int t) const {
    const int num_worker_threads = m_policy.space().concurrency();
    const int value_count        = Analysis::value_count(m_functor);
    const std::size_t value_size = Analysis::value_size(m_functor);

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    typename Analysis::Reducer final_reducer(m_functor);
    barrier_type &barrier =
        *static_cast<barrier_type *>(buffer.get_extra_space());
    reference_type update_sum =
        final_reducer.init(reinterpret_cast<pointer_type>(buffer.get(t)));

    const WorkRange range(m_policy, t, num_worker_threads);
    execute_chunk(range.begin(), range.end(), update_sum, false);

    {
      // Since arrive_and_wait may yield and resume on another worker thread we
      // set in_parallel = false on the current thread before suspending and set
      // it again to true when we resume.
      Kokkos::Experimental::HPX::impl_not_in_parallel_scope p;
      barrier.arrive_and_wait();
    }

    if (t == 0) {
      final_reducer.init(reinterpret_cast<pointer_type>(
          static_cast<char *>(buffer.get(0)) + value_size));

      for (int i = 1; i < num_worker_threads; ++i) {
        pointer_type ptr_1_prev =
            reinterpret_cast<pointer_type>(buffer.get(i - 1));
        pointer_type ptr_2_prev = reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(i - 1)) + value_size);
        pointer_type ptr_2 = reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(i)) + value_size);

        for (int j = 0; j < value_count; ++j) {
          ptr_2[j] = ptr_2_prev[j];
        }

        final_reducer.join(ptr_2, ptr_1_prev);
      }
    }

    {
      // Since arrive_and_wait may yield and resume on another worker thread we
      // set in_parallel = false on the current thread before suspending and set
      // it again to true when we resume.
      Kokkos::Experimental::HPX::impl_not_in_parallel_scope p;
      barrier.arrive_and_wait();
    }

    reference_type update_base =
        Analysis::Reducer::reference(reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(t)) + value_size));

    execute_chunk(range.begin(), range.end(), update_base, true);
  }

  void finalize() const {
    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    static_cast<barrier_type *>(buffer.get_extra_space())->~barrier_type();
  }

  void execute() const {
    const int num_worker_threads = m_policy.space().concurrency();
    m_policy.space().impl_bulk_setup_finalize(
        false, is_light_weight_policy<Policy>(), *this, num_worker_threads,
        hpx::threads::thread_stacksize::small_);
  }

  inline ParallelScan(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::HPX> {
 private:
  using Policy         = Kokkos::RangePolicy<Traits...>;
  using WorkTag        = typename Policy::work_tag;
  using WorkRange      = typename Policy::WorkRange;
  using Member         = typename Policy::member_type;
  using Analysis       = FunctorAnalysis<FunctorPatternInterface::SCAN, Policy,
                                   FunctorType, ReturnType>;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;
  using value_type     = typename Analysis::value_type;
  using barrier_type   = hpx::barrier<>;

  const FunctorType m_functor;
  const Policy m_policy;
  pointer_type m_result_ptr;

 public:
  void setup() const {
    const int num_worker_threads = m_policy.space().concurrency();
    const std::size_t value_size = Analysis::value_size(m_functor);

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    buffer.resize(num_worker_threads, 2 * value_size, sizeof(barrier_type));

    new (buffer.get_extra_space()) barrier_type(num_worker_threads);
  }

  void execute_chunk(const Member i_begin, const Member i_end,
                     reference_type update, const bool final) const {
    for (Member i = i_begin; i < i_end; ++i) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor(i, update, final);
      } else {
        m_functor(WorkTag{}, i, update, final);
      }
    }
  }

  void execute_range(int t) const {
    const int num_worker_threads = m_policy.space().concurrency();
    const int value_count        = Analysis::value_count(m_functor);
    const std::size_t value_size = Analysis::value_size(m_functor);

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    typename Analysis::Reducer final_reducer(m_functor);
    barrier_type &barrier =
        *static_cast<barrier_type *>(buffer.get_extra_space());
    reference_type update_sum =
        final_reducer.init(reinterpret_cast<pointer_type>(buffer.get(t)));

    const WorkRange range(m_policy, t, num_worker_threads);
    execute_chunk(range.begin(), range.end(), update_sum, false);

    {
      // Since arrive_and_wait may yield and resume on another worker thread we
      // set in_parallel = false on the current thread before suspending and set
      // it again to true when we resume.
      Kokkos::Experimental::HPX::impl_not_in_parallel_scope p;
      barrier.arrive_and_wait();
    }

    if (t == 0) {
      final_reducer.init(reinterpret_cast<pointer_type>(
          static_cast<char *>(buffer.get(0)) + value_size));

      for (int i = 1; i < num_worker_threads; ++i) {
        pointer_type ptr_1_prev =
            reinterpret_cast<pointer_type>(buffer.get(i - 1));
        pointer_type ptr_2_prev = reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(i - 1)) + value_size);
        pointer_type ptr_2 = reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(i)) + value_size);

        for (int j = 0; j < value_count; ++j) {
          ptr_2[j] = ptr_2_prev[j];
        }

        final_reducer.join(ptr_2, ptr_1_prev);
      }
    }

    {
      // Since arrive_and_wait may yield and resume on another worker thread we
      // set in_parallel = false on the current thread before suspending and set
      // it again to true when we resume.
      Kokkos::Experimental::HPX::impl_not_in_parallel_scope p;
      barrier.arrive_and_wait();
    }

    reference_type update_base =
        Analysis::Reducer::reference(reinterpret_cast<pointer_type>(
            static_cast<char *>(buffer.get(t)) + value_size));

    execute_chunk(range.begin(), range.end(), update_base, true);

    if (t == num_worker_threads - 1) {
      *m_result_ptr = update_base;
    }
  }

  void finalize() const {
    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    static_cast<barrier_type *>(buffer.get_extra_space())->~barrier_type();
  }

  void execute() const {
    const int num_worker_threads = m_policy.space().concurrency();
    m_policy.space().impl_bulk_setup_finalize(
        false, is_light_weight_policy<Policy>(), *this, num_worker_threads,
        hpx::threads::thread_stacksize::small_);
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType &arg_functor,
                        const Policy &arg_policy,
                        const ViewType &arg_result_view)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::HPX parallel_scan result must be host-accessible!");
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {
template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Experimental::HPX> {
 private:
  using Policy  = TeamPolicyInternal<Kokkos::Experimental::HPX, Properties...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;
  using memory_space = Kokkos::HostSpace;

  const FunctorType m_functor;
  const Policy m_policy;
  const int m_league;
  const std::size_t m_shared;

 public:
  void setup() const {
    const int num_worker_threads = m_policy.space().concurrency();

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    buffer.resize(num_worker_threads, m_shared);
  }

  void execute_range(const int i) const {
    const int t = Kokkos::Experimental::HPX::impl_hardware_thread_id();
    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    const auto r =
        get_chunk_range(i, 0, m_policy.chunk_size(), m_policy.league_size());
    for (int league_rank = r.begin; league_rank < r.end; ++league_rank) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor(Member(m_policy, 0, league_rank, buffer.get(t), m_shared));
      } else {
        m_functor(WorkTag{},
                  Member(m_policy, 0, league_rank, buffer.get(t), m_shared));
      }
    }
  }

  void finalize() const {}

  void execute() const {
    const int num_chunks =
        get_num_chunks(0, m_policy.chunk_size(), m_policy.league_size());
    m_policy.space().impl_bulk_setup_finalize(
        false, is_light_weight_policy<Policy>(), *this, num_chunks,
        hpx::threads::thread_stacksize::nostack);
  }

  ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league(arg_policy.league_size()),
        m_shared(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     arg_functor, arg_policy.team_size())) {}
};

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>,
                     Kokkos::Experimental::HPX> {
 private:
  using Policy = TeamPolicyInternal<Kokkos::Experimental::HPX, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using Member  = typename Policy::member_type;
  using WorkTag = typename Policy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const int m_league;
  const Policy m_policy;
  pointer_type m_result_ptr;
  const std::size_t m_shared;
  const bool m_force_synchronous;

 public:
  void setup() const {
    const ReducerType &reducer   = m_functor_reducer.get_reducer();
    const std::size_t value_size = reducer.value_size();
    const int num_worker_threads = m_policy.space().concurrency();

    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    buffer.resize(num_worker_threads, value_size + m_shared);

    for (int t = 0; t < num_worker_threads; ++t) {
      reducer.init(reinterpret_cast<pointer_type>(buffer.get(t)));
    }
  }

  void execute_range(const int i) const {
    const ReducerType &reducer   = m_functor_reducer.get_reducer();
    const std::size_t value_size = reducer.value_size();
    std::size_t t = Kokkos::Experimental::HPX::impl_hardware_thread_id();
    hpx_thread_buffer &buffer = m_policy.space().impl_get_buffer();
    reference_type update =
        ReducerType::reference(reinterpret_cast<pointer_type>(buffer.get(t)));
    const auto r =
        get_chunk_range(i, 0, m_policy.chunk_size(), m_policy.league_size());
    char *local_buffer = static_cast<char *>(buffer.get(t)) + value_size;
    for (int league_rank = r.begin; league_rank < r.end; ++league_rank) {
      if constexpr (std::is_same_v<WorkTag, void>) {
        m_functor_reducer.get_functor()(
            Member(m_policy, 0, league_rank, local_buffer, m_shared), update);
      } else {
        m_functor_reducer.get_functor()(
            WorkTag{}, Member(m_policy, 0, league_rank, local_buffer, m_shared),
            update);
      }
    }
  }

  void finalize() const {
    hpx_thread_buffer &buffer    = m_policy.space().impl_get_buffer();
    const ReducerType &reducer   = m_functor_reducer.get_reducer();
    const int num_worker_threads = m_policy.space().concurrency();
    const pointer_type ptr = reinterpret_cast<pointer_type>(buffer.get(0));
    for (int t = 1; t < num_worker_threads; ++t) {
      reducer.join(ptr, reinterpret_cast<pointer_type>(buffer.get(t)));
    }

    reducer.final(ptr);

    if (m_result_ptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = ptr[j];
      }
    }
  }

  void execute() const {
    if (m_policy.league_size() * m_policy.team_size() == 0) {
      if (m_result_ptr) {
        const ReducerType &reducer = m_functor_reducer.get_reducer();
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
      return;
    }

    const int num_chunks =
        get_num_chunks(0, m_policy.chunk_size(), m_policy.league_size());
    m_policy.space().impl_bulk_setup_finalize(
        m_force_synchronous, is_light_weight_policy<Policy>(), *this,
        num_chunks, hpx::threads::thread_stacksize::nostack);
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType &arg_functor_reducer,
                 const Policy &arg_policy, const ViewType &arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_league(arg_policy.league_size()),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_shared(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     m_functor_reducer.get_functor(), arg_policy.team_size())),
        m_force_synchronous(!arg_result.impl_track().has_record()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "HPX reduce result must be a View accessible from HostSpace");
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>
    TeamThreadRange(const Impl::HPXTeamMember &thread, const iType &count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::HPXTeamMember>
TeamThreadRange(const Impl::HPXTeamMember &thread, const iType1 &i_begin,
                const iType2 &i_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, iType(i_begin), iType(i_end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>
    TeamVectorRange(const Impl::HPXTeamMember &thread, const iType &count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::HPXTeamMember>
TeamVectorRange(const Impl::HPXTeamMember &thread, const iType1 &i_begin,
                const iType2 &i_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, iType(i_begin), iType(i_end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
    ThreadVectorRange(const Impl::HPXTeamMember &thread, const iType &count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::HPXTeamMember>
ThreadVectorRange(const Impl::HPXTeamMember &thread, const iType1 &i_begin,
                  const iType2 &i_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>(
      thread, iType(i_begin), iType(i_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::HPXTeamMember> PerTeam(
    const Impl::HPXTeamMember &thread) {
  return Impl::ThreadSingleStruct<Impl::HPXTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::HPXTeamMember> PerThread(
    const Impl::HPXTeamMember &thread) {
  return Impl::VectorSingleStruct<Impl::HPXTeamMember>(thread);
}

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType,
          typename = std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda, ValueType &result) {
  result = ValueType();
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }
}

/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i);
  }
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType,
          typename = std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda, ValueType &result) {
  result = ValueType();
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }
}

template <typename iType, class Lambda, typename ReducerType,
          typename = std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda, const ReducerType &reducer) {
  reducer.init(reducer.reference());
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
}

template <typename iType, class Lambda, typename ReducerType,
          typename = std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const Lambda &lambda, const ReducerType &reducer) {
  reducer.init(reducer.reference());
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
}

template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HPXTeamMember> const
        &loop_boundaries,
    const FunctorType &lambda) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, FunctorType,
      void>::value_type;

  value_type scan_val = value_type();

  // Intra-member scan
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, false);
  }

  // 'scan_val' output is the exclusive prefix sum
  scan_val = loop_boundaries.thread.team_scan(scan_val);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, true);
  }
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes
 * lambda(iType i, ValueType & val, bool final) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan
 * operation is performed. Depending on the target execution space the operator
 * might be called twice: once with final=false and once with final=true. When
 * final==true val contains the prefix sum value. The contribution of this "i"
 * needs to be added to val no matter whether final==true or not. In a serial
 * execution (i.e. team_size==1) the operator is only called once with
 * final==true. Scan_val will be set to the final sum value over all vector
 */
template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const FunctorType &lambda) {
  using value_type =
      typename Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                     TeamPolicy<Experimental::HPX>, FunctorType,
                                     void>::value_type;

  value_type scan_val = value_type();

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, true);
  }
}

/** \brief  Intra-thread vector parallel scan with reducer
 *
 */
template <typename iType, class FunctorType, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HPXTeamMember>
        &loop_boundaries,
    const FunctorType &lambda, const ReducerType &reducer) {
  typename ReducerType::value_type scan_val;
  reducer.init(scan_val);

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, true);
  }
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::HPXTeamMember> &,
    const FunctorType &lambda) {
  lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::HPXTeamMember> &,
    const FunctorType &lambda) {
  lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::HPXTeamMember> &,
    const FunctorType &lambda, ValueType &val) {
  lambda(val);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::HPXTeamMember> &,
    const FunctorType &lambda, ValueType &val) {
  lambda(val);
}

}  // namespace Kokkos

#include <HPX/Kokkos_HPX_Task.hpp>

#endif /* #if defined( KOKKOS_ENABLE_HPX ) */
#endif /* #ifndef KOKKOS_HPX_HPP */
