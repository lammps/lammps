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

#ifndef KOKKOS_SYCL_TEAM_POLICY_HPP
#define KOKKOS_SYCL_TEAM_POLICY_HPP

#include <SYCL/Kokkos_SYCL_Team.hpp>

#include <vector>

template <typename... Properties>
class Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::SYCL,
                                       Properties...>
    : public PolicyTraits<Properties...> {
 public:
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  template <typename ExecSpace, typename... OtherProperties>
  friend class TeamPolicyInternal;

 private:
  typename traits::execution_space m_space;
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];
  int m_chunk_size;
  bool m_tune_team_size;
  bool m_tune_vector_length;

 public:
  using execution_space = Kokkos::Experimental::SYCL;

  template <class... OtherProperties>
  TeamPolicyInternal(TeamPolicyInternal<OtherProperties...> const& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    m_space                  = p.m_space;
    m_tune_team_size         = p.m_tune_team_size;
    m_tune_vector_length     = p.m_tune_vector_length;
  }

  template <typename FunctorType>
  int team_size_max(FunctorType const& f, ParallelForTag const&) const {
    return internal_team_size_max_for(f);
  }

  template <class FunctorType>
  inline int team_size_max(const FunctorType& f,
                           const ParallelReduceTag&) const {
    return internal_team_size_max_reduce<void>(f);
  }

  template <class FunctorType, class ReducerType>
  inline int team_size_max(const FunctorType& f, const ReducerType& /*r*/,
                           const ParallelReduceTag&) const {
    return internal_team_size_max_reduce<typename ReducerType::value_type>(f);
  }

  template <typename FunctorType>
  int team_size_recommended(FunctorType const& f, ParallelForTag const&) const {
    return internal_team_size_recommended_for(f);
  }

  template <typename FunctorType>
  inline int team_size_recommended(FunctorType const& f,
                                   ParallelReduceTag const&) const {
    return internal_team_size_recommended_reduce<void>(f);
  }

  template <class FunctorType, class ReducerType>
  int team_size_recommended(FunctorType const& f, ReducerType const&,
                            ParallelReduceTag const&) const {
    return internal_team_size_recommended_reduce<
        typename ReducerType::value_type>(f);
  }
  inline bool impl_auto_vector_length() const { return m_tune_vector_length; }
  inline bool impl_auto_team_size() const { return m_tune_team_size; }
  // FIXME_SYCL This is correct in most cases, but not necessarily in case a
  // custom sycl::queue is used to initialize the execution space.
  static int vector_length_max() {
    std::vector<size_t> sub_group_sizes =
        execution_space{}
            .impl_internal_space_instance()
            ->m_queue->get_device()
            .template get_info<sycl::info::device::sub_group_sizes>();
    return *std::max_element(sub_group_sizes.begin(), sub_group_sizes.end());
  }

 private:
  static int verify_requested_vector_length(int requested_vector_length) {
    int test_vector_length =
        std::min(requested_vector_length, vector_length_max());

    // Allow only power-of-two vector_length
    if (!(is_integral_power_of_two(test_vector_length))) {
      int test_pow2 = 1;
      while (test_pow2 < test_vector_length) test_pow2 <<= 1;
      test_vector_length = test_pow2 >> 1;
    }

    return test_vector_length;
  }

 public:
  static int scratch_size_max(int level) {
    return level == 0 ? 1024 * 32
                      :           // FIXME_SYCL arbitrarily setting this to 32kB
               20 * 1024 * 1024;  // FIXME_SYCL arbitrarily setting this to 20MB
  }
  inline void impl_set_vector_length(size_t size) { m_vector_length = size; }
  inline void impl_set_team_size(size_t size) { m_team_size = size; }
  int impl_vector_length() const { return m_vector_length; }

  int team_size() const { return m_team_size; }

  int league_size() const { return m_league_size; }

  size_t scratch_size(int level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  size_t team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }

  size_t thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  typename traits::execution_space space() const { return m_space; }

  TeamPolicyInternal()
      : m_space(typename traits::execution_space()),
        m_league_size(0),
        m_team_size(-1),
        m_vector_length(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(vector_length_max()),
        m_tune_team_size(false),
        m_tune_vector_length(false) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request, int vector_length_request = 1)
      : m_space(space_),
        m_league_size(league_size_),
        m_team_size(team_size_request),
        m_vector_length(
            (vector_length_request > 0)
                ? verify_requested_vector_length(vector_length_request)
                : (verify_requested_vector_length(1))),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(vector_length_max()),
        m_tune_team_size(bool(team_size_request <= 0)),
        m_tune_vector_length(bool(vector_length_request <= 0)) {
    // FIXME_SYCL Check that league size is permissible,
    // https://github.com/intel/llvm/pull/4064

    // Make sure total block size is permissible
    if (m_team_size * m_vector_length >
        static_cast<int>(
            m_space.impl_internal_space_instance()->m_maxWorkgroupSize)) {
      Impl::throw_runtime_exception(
          std::string("Kokkos::TeamPolicy<SYCL> the team size is too large. "
                      "Team size x vector length is " +
                      std::to_string(m_team_size * m_vector_length) +
                      " but must be smaller than ") +
          std::to_string(
              m_space.impl_internal_space_instance()->m_maxWorkgroupSize));
    }
  }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : TeamPolicyInternal(space_, league_size_, -1, vector_length_request) {}
  // FLAG
  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */
                     )
      : TeamPolicyInternal(space_, league_size_, team_size_request, -1)

  {}

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(space_, league_size_, -1, -1)

  {}

  TeamPolicyInternal(int league_size_, int team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request) {}

  TeamPolicyInternal(int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_, -1,
                           vector_length_request) {}

  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(int league_size_, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, -1)

  {}

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(typename traits::execution_space(), league_size_, -1,
                           -1) {}

  int chunk_size() const { return m_chunk_size; }

  TeamPolicyInternal& set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  TeamPolicyInternal& set_scratch_size(int level,
                                       PerTeamValue const& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  TeamPolicyInternal& set_scratch_size(int level,
                                       PerThreadValue const& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  TeamPolicyInternal& set_scratch_size(int level, PerTeamValue const& per_team,
                                       PerThreadValue const& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  using member_type = Kokkos::Impl::SYCLTeamMember;

 protected:
  template <class FunctorType>
  int internal_team_size_max_for(const FunctorType& /*f*/) const {
    // nested_reducer_memsize = (sizeof(double) * (m_team_size + 2)
    // custom: m_team_scratch_size[0] + m_thread_scratch_size[0] * m_team_size
    // total:
    // 2*sizeof(double)+m_team_scratch_size[0]
    // + m_team_size(sizeof(double)+m_thread_scratch_size[0])
    const int max_threads_for_memory =
        (space().impl_internal_space_instance()->m_maxShmemPerBlock -
         2 * sizeof(double) - m_team_scratch_size[0]) /
        (sizeof(double) + m_thread_scratch_size[0]);
    return std::min({
             int(m_space.impl_internal_space_instance()->m_maxWorkgroupSize),
      // FIXME_SYCL Avoid requesting too many registers on NVIDIA GPUs.
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
                 256,
#endif
                 max_threads_for_memory
           }) /
           impl_vector_length();
  }

  template <class ValueType, class FunctorType>
  int internal_team_size_max_reduce(const FunctorType& f) const {
    using Analysis =
        FunctorAnalysis<FunctorPatternInterface::REDUCE, TeamPolicyInternal,
                        FunctorType, ValueType>;
    using value_type      = typename Analysis::value_type;
    const int value_count = Analysis::value_count(f);

    // nested_reducer_memsize = (sizeof(double) * (m_team_size + 2)
    // reducer_memsize = sizeof(value_type) * m_team_size * value_count
    // custom: m_team_scratch_size[0] + m_thread_scratch_size[0] * m_team_size
    // total:
    // 2*sizeof(double)+m_team_scratch_size[0]
    // + m_team_size(sizeof(double)+sizeof(value_type)*value_count
    //               +m_thread_scratch_size[0])
    const int max_threads_for_memory =
        (space().impl_internal_space_instance()->m_maxShmemPerBlock -
         2 * sizeof(double) - m_team_scratch_size[0]) /
        (sizeof(double) + sizeof(value_type) * value_count +
         m_thread_scratch_size[0]);
    return std::min<int>({
             int(m_space.impl_internal_space_instance()->m_maxWorkgroupSize),
      // FIXME_SYCL Avoid requesting too many registers on NVIDIA GPUs.
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
                 256,
#endif
                 max_threads_for_memory
           }) /
           impl_vector_length();
  }

  template <class FunctorType>
  int internal_team_size_recommended_for(const FunctorType& f) const {
    // FIXME_SYCL improve
    return 1 << Kokkos::Impl::int_log2(internal_team_size_max_for(f));
  }

  template <class ValueType, class FunctorType>
  int internal_team_size_recommended_reduce(const FunctorType& f) const {
    // FIXME_SYCL improve
    return 1 << Kokkos::Impl::int_log2(
               internal_team_size_max_reduce<ValueType>(f));
  }
};

#endif
