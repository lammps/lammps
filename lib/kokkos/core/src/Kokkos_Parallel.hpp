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

/// \file Kokkos_Parallel.hpp
/// \brief Declaration of parallel operators

#ifndef KOKKOS_PARALLEL_HPP
#define KOKKOS_PARALLEL_HPP

#include <cstddef>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_ExecPolicy.hpp>

#include <impl/Kokkos_Tools.hpp>
#include <type_traits>
#include <typeinfo>

#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

#ifdef KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
#include <iostream>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class T, class = void>
struct is_detected_execution_space : std::false_type {
  using type = not_a_type;
};

template <class T>
struct is_detected_execution_space<T, void_t<typename T::execution_space>>
    : std::true_type {
  using type = typename T::execution_space;
};

template <class T>
using detected_execution_space_t =
    typename is_detected_execution_space<T>::type;

template <class T, class = void>
struct is_detected_device_type : std::false_type {
  using type = not_a_type;
};

template <class T>
struct is_detected_device_type<T, void_t<typename T::device_type>>
    : std::true_type {
  using type = typename T::device_type;
};

template <class T>
using detected_device_type_t = typename is_detected_device_type<T>::type;

//----------------------------------------------------------------------------
/** \brief  Given a Functor and Execution Policy query an execution space.
 *
 *  if       the Policy has an execution space use that
 *  else if  the Functor has an execution_space use that
 *  else if  the Functor has a device_type use that for backward compatibility
 *  else     use the default
 */

template <class Functor, class Policy>
struct FunctorPolicyExecutionSpace {
  using execution_space = std::conditional_t<
      is_detected_execution_space<Policy>::value,
      detected_execution_space_t<Policy>,
      std::conditional_t<
          is_detected_execution_space<Functor>::value,
          detected_execution_space_t<Functor>,
          std::conditional_t<
              is_detected_device_type<Functor>::value,
              detected_execution_space_t<detected_device_type_t<Functor>>,
              Kokkos::DefaultExecutionSpace>>>;
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief Execute \c functor in parallel according to the execution \c policy.
 *
 * A "functor" is a class containing the function to execute in parallel,
 * data needed for that execution, and an optional \c execution_space
 * alias.  Here is an example functor for parallel_for:
 *
 * \code
 *  class FunctorType {
 *  public:
 *    using execution_space = ...;
 *    void operator() ( WorkType iwork ) const ;
 *  };
 * \endcode
 *
 * In the above example, \c WorkType is any integer type for which a
 * valid conversion from \c size_t to \c IntType exists.  Its
 * <tt>operator()</tt> method defines the operation to parallelize,
 * over the range of integer indices <tt>iwork=[0,work_count-1]</tt>.
 * This compares to a single iteration \c iwork of a \c for loop.
 * If \c execution_space is not defined DefaultExecutionSpace will be used.
 */
template <class ExecPolicy, class FunctorType>
inline void parallel_for(
    const ExecPolicy& policy, const FunctorType& functor,
    const std::string& str = "",
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<ExecPolicy>::value>::type* =
        nullptr) {
  uint64_t kpID = 0;

  ExecPolicy inner_policy = policy;
  Kokkos::Tools::Impl::begin_parallel_for(inner_policy, functor, str, kpID);

  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelFor<FunctorType, ExecPolicy> closure(functor, inner_policy);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_for(inner_policy, functor, str, kpID);
}

template <class FunctorType>
inline void parallel_for(const size_t work_count, const FunctorType& functor,
                         const std::string& str = "") {
  using execution_space =
      typename Impl::FunctorPolicyExecutionSpace<FunctorType,
                                                 void>::execution_space;
  using policy = RangePolicy<execution_space>;

  uint64_t kpID = 0;

  policy execution_policy = policy(0, work_count);

  Kokkos::Tools::Impl::begin_parallel_for(execution_policy, functor, str, kpID);

  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelFor<FunctorType, policy> closure(functor, execution_policy);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_for(execution_policy, functor, str, kpID);
}

template <class ExecPolicy, class FunctorType>
inline void parallel_for(const std::string& str, const ExecPolicy& policy,
                         const FunctorType& functor) {
#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_for kernel: " << str << std::endl;
#endif

  ::Kokkos::parallel_for(policy, functor, str);

#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End   parallel_for kernel: " << str << std::endl;
#endif
  (void)str;
}

}  // namespace Kokkos

#include <Kokkos_Parallel_Reduce.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/// \fn parallel_scan
/// \tparam ExecutionPolicy The execution policy type.
/// \tparam FunctorType     The scan functor type.
///
/// \param policy  [in] The execution policy.
/// \param functor [in] The scan functor.
///
/// This function implements a parallel scan pattern.  The scan can
/// be either inclusive or exclusive, depending on how you implement
/// the scan functor.
///
/// A scan functor looks almost exactly like a reduce functor, except
/// that its operator() takes a third \c bool argument, \c final_pass,
/// which indicates whether this is the last pass of the scan
/// operation.  We will show below how to use the \c final_pass
/// argument to control whether the scan is inclusive or exclusive.
///
/// Here is the minimum required interface of a scan functor for a POD
/// (plain old data) value type \c PodType.  That is, the result is a
/// View of zero or more PodType.  It is also possible for the result
/// to be an array of (same-sized) arrays of PodType, but we do not
/// show the required interface for that here.
/// \code
/// template< class ExecPolicy , class FunctorType >
/// class ScanFunctor {
/// public:
///   // The Kokkos device type
///   using execution_space = ...;
///   // Type of an entry of the array containing the result;
///   // also the type of each of the entries combined using
///   // operator() or join().
///   using value_type = PodType;
///
///   void operator () (const ExecPolicy::member_type & i, value_type& update,
///   const bool final_pass) const; void init (value_type& update) const; void
///   join (volatile value_type& update, volatile const value_type& input) const
/// };
/// \endcode
///
/// Here is an example of a functor which computes an inclusive plus-scan
/// of an array of \c int, in place.  If given an array [1, 2, 3, 4], this
/// scan will overwrite that array with [1, 3, 6, 10].
///
/// \code
/// template<class SpaceType>
/// class InclScanFunctor {
/// public:
///   using execution_space = SpaceType;
///   using value_type = int;
///   using size_type = typename SpaceType::size_type;
///
///   InclScanFunctor( Kokkos::View<value_type*, execution_space> x
///                  , Kokkos::View<value_type*, execution_space> y ) : m_x(x),
///                  m_y(y) {}
///
///   void operator () (const size_type i, value_type& update, const bool
///   final_pass) const {
///     update += m_x(i);
///     if (final_pass) {
///       m_y(i) = update;
///     }
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input)
///   const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, execution_space> m_x;
///   Kokkos::View<value_type*, execution_space> m_y;
/// };
/// \endcode
///
/// Here is an example of a functor which computes an <i>exclusive</i>
/// scan of an array of \c int, in place.  In operator(), note both
/// that the final_pass test and the update have switched places, and
/// the use of a temporary.  If given an array [1, 2, 3, 4], this scan
/// will overwrite that array with [0, 1, 3, 6].
///
/// \code
/// template<class SpaceType>
/// class ExclScanFunctor {
/// public:
///   using execution_space = SpaceType;
///   using value_type = int;
///   using size_type = typename SpaceType::size_type;
///
///   ExclScanFunctor (Kokkos::View<value_type*, execution_space> x) : x_ (x) {}
///
///   void operator () (const size_type i, value_type& update, const bool
///   final_pass) const {
///     const value_type x_i = x_(i);
///     if (final_pass) {
///       x_(i) = update;
///     }
///     update += x_i;
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input)
///   const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, execution_space> x_;
/// };
/// \endcode
///
/// Here is an example of a functor which builds on the above
/// exclusive scan example, to compute an offsets array from a
/// population count array, in place.  We assume that the pop count
/// array has an extra entry at the end to store the final count.  If
/// given an array [1, 2, 3, 4, 0], this scan will overwrite that
/// array with [0, 1, 3, 6, 10].
///
/// \code
/// template<class SpaceType>
/// class OffsetScanFunctor {
/// public:
///   using execution_space = SpaceType;
///   using value_type = int;
///   using size_type = typename SpaceType::size_type;
///
///   // lastIndex_ is the last valid index (zero-based) of x.
///   // If x has length zero, then lastIndex_ won't be used anyway.
///   OffsetScanFunctor( Kokkos::View<value_type*, execution_space> x
///                    , Kokkos::View<value_type*, execution_space> y )
///      : m_x(x), m_y(y), last_index_ (x.dimension_0 () == 0 ? 0 :
///      x.dimension_0 () - 1)
///   {}
///
///   void operator () (const size_type i, int& update, const bool final_pass)
///   const {
///     if (final_pass) {
///       m_y(i) = update;
///     }
///     update += m_x(i);
///     // The last entry of m_y gets the final sum.
///     if (final_pass && i == last_index_) {
///       m_y(i+1) = update;
// i/     }
///   }
///   void init (value_type& update) const {
///     update = 0;
///   }
///   void join (volatile value_type& update, volatile const value_type& input)
///   const {
///     update += input;
///   }
///
/// private:
///   Kokkos::View<value_type*, execution_space> m_x;
///   Kokkos::View<value_type*, execution_space> m_y;
///   const size_type last_index_;
/// };
/// \endcode
///
template <class ExecutionPolicy, class FunctorType>
inline void parallel_scan(
    const ExecutionPolicy& policy, const FunctorType& functor,
    const std::string& str = "",
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<ExecutionPolicy>::value>::type* =
        nullptr) {
  uint64_t kpID                = 0;
  ExecutionPolicy inner_policy = policy;
  Kokkos::Tools::Impl::begin_parallel_scan(inner_policy, functor, str, kpID);

  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelScan<FunctorType, ExecutionPolicy> closure(functor,
                                                           inner_policy);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_scan(inner_policy, functor, str, kpID);
}

template <class FunctorType>
inline void parallel_scan(const size_t work_count, const FunctorType& functor,
                          const std::string& str = "") {
  using execution_space =
      typename Kokkos::Impl::FunctorPolicyExecutionSpace<FunctorType,
                                                         void>::execution_space;

  using policy = Kokkos::RangePolicy<execution_space>;

  uint64_t kpID = 0;
  policy execution_policy(0, work_count);
  Kokkos::Tools::Impl::begin_parallel_scan(execution_policy, functor, str,
                                           kpID);
  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelScan<FunctorType, policy> closure(functor, execution_policy);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_scan(execution_policy, functor, str, kpID);
}

template <class ExecutionPolicy, class FunctorType>
inline void parallel_scan(const std::string& str, const ExecutionPolicy& policy,
                          const FunctorType& functor) {
#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_scan kernel: " << str << std::endl;
#endif

  ::Kokkos::parallel_scan(policy, functor, str);

#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End parallel_scan kernel: " << str << std::endl;
#endif
  (void)str;
}

template <class ExecutionPolicy, class FunctorType, class ReturnType>
inline void parallel_scan(
    const ExecutionPolicy& policy, const FunctorType& functor,
    ReturnType& return_value, const std::string& str = "",
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<ExecutionPolicy>::value>::type* =
        nullptr) {
  uint64_t kpID                = 0;
  ExecutionPolicy inner_policy = policy;
  Kokkos::Tools::Impl::begin_parallel_scan(inner_policy, functor, str, kpID);

  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelScanWithTotal<FunctorType, ExecutionPolicy, ReturnType> closure(
      functor, inner_policy, return_value);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_scan(inner_policy, functor, str, kpID);

  policy.space().fence();
}

template <class FunctorType, class ReturnType>
inline void parallel_scan(const size_t work_count, const FunctorType& functor,
                          ReturnType& return_value,
                          const std::string& str = "") {
  using execution_space =
      typename Kokkos::Impl::FunctorPolicyExecutionSpace<FunctorType,
                                                         void>::execution_space;

  using policy = Kokkos::RangePolicy<execution_space>;

  policy execution_policy(0, work_count);
  uint64_t kpID = 0;
  Kokkos::Tools::Impl::begin_parallel_scan(execution_policy, functor, str,
                                           kpID);

  Kokkos::Impl::shared_allocation_tracking_disable();
  Impl::ParallelScanWithTotal<FunctorType, policy, ReturnType> closure(
      functor, execution_policy, return_value);
  Kokkos::Impl::shared_allocation_tracking_enable();

  closure.execute();

  Kokkos::Tools::Impl::end_parallel_scan(execution_policy, functor, str, kpID);

  execution_space().fence();
}

template <class ExecutionPolicy, class FunctorType, class ReturnType>
inline void parallel_scan(const std::string& str, const ExecutionPolicy& policy,
                          const FunctorType& functor,
                          ReturnType& return_value) {
#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG Start parallel_scan kernel: " << str << std::endl;
#endif

  ::Kokkos::parallel_scan(policy, functor, return_value, str);

#if KOKKOS_ENABLE_DEBUG_PRINT_KERNEL_NAMES
  Kokkos::fence();
  std::cout << "KOKKOS_DEBUG End parallel_scan kernel: " << str << std::endl;
#endif
  (void)str;
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType,
          bool HasTeamShmemSize =
              has_member_team_shmem_size<FunctorType>::value,
          bool HasShmemSize = has_member_shmem_size<FunctorType>::value>
struct FunctorTeamShmemSize {
  KOKKOS_INLINE_FUNCTION static size_t value(const FunctorType&, int) {
    return 0;
  }
};

template <class FunctorType>
struct FunctorTeamShmemSize<FunctorType, true, false> {
  static inline size_t value(const FunctorType& f, int team_size) {
    return f.team_shmem_size(team_size);
  }
};

template <class FunctorType>
struct FunctorTeamShmemSize<FunctorType, false, true> {
  static inline size_t value(const FunctorType& f, int team_size) {
    return f.shmem_size(team_size);
  }
};
template <class FunctorType>
struct FunctorTeamShmemSize<FunctorType, true, true> {
  static inline size_t value(const FunctorType& /*f*/, int /*team_size*/) {
    Kokkos::abort(
        "Functor with both team_shmem_size and shmem_size defined is "
        "not allowed");
    return 0;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_PARALLEL_HPP */
