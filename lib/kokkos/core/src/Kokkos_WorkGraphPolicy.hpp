/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_WORKGRAPHPOLICY_HPP
#define KOKKOS_WORKGRAPHPOLICY_HPP

namespace Kokkos {
namespace Impl {
namespace Experimental {

template< class functor_type , class execution_space, class ... policy_args >
class WorkGraphExec;

}}} // namespace Kokkos::Impl::Experimental

namespace Kokkos {
namespace Experimental {

template< class ... Properties >
class WorkGraphPolicy
{
public:

  using self_type = WorkGraphPolicy<Properties ... >;
  using traits = Kokkos::Impl::PolicyTraits<Properties ... >;
  using index_type = typename traits::index_type;
  using execution_space = typename traits::execution_space;
  using work_tag = typename traits::work_tag;
  using memory_space = typename execution_space::memory_space;
  using graph_type = Kokkos::Experimental::Crs<index_type, execution_space, void, index_type>;
  using member_type = index_type;

private:
   
  graph_type m_graph;

  using ints_type = Kokkos::View<std::int32_t*, memory_space>;
  using range_type = Kokkos::pair<std::int32_t, std::int32_t>;
  using ranges_type = Kokkos::View<range_type*, memory_space>;
  const std::int32_t m_total_work;
  ints_type m_counts;
  ints_type m_queue;
  ranges_type m_ranges;

public:

  struct TagZeroRanges {};
  KOKKOS_INLINE_FUNCTION
  void operator()(TagZeroRanges, std::int32_t i) const {
    m_ranges[i] = range_type(0, 0);
  }
  void zero_ranges() {
    using policy_type = RangePolicy<std::int32_t, execution_space, TagZeroRanges>;
    using closure_type = Kokkos::Impl::ParallelFor<self_type, policy_type>;
    const closure_type closure(*this, policy_type(0, 1));
    closure.execute();
    execution_space::fence();
  }

  struct TagFillQueue {};
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFillQueue, std::int32_t i) const {
    if (*((volatile std::int32_t*)(&m_counts(i))) == 0) push_work(i);
  }
  void fill_queue() {
    using policy_type = RangePolicy<std::int32_t, execution_space, TagFillQueue>;
    using closure_type = Kokkos::Impl::ParallelFor<self_type, policy_type>;
    const closure_type closure(*this, policy_type(0, m_total_work));
    closure.execute();
    execution_space::fence();
  }

private:

  inline
  void setup() {
    if (m_graph.numRows() > std::numeric_limits<std::int32_t>::max()) {
      Kokkos::abort("WorkGraphPolicy work must be indexable using int32_t");
    }
    get_crs_transpose_counts(m_counts, m_graph);
    m_queue = ints_type(ViewAllocateWithoutInitializing("queue"), m_total_work);
    deep_copy(m_queue, std::int32_t(-1));
    m_ranges = ranges_type("ranges", 1);
    fill_queue();
  }

  KOKKOS_INLINE_FUNCTION
  std::int32_t pop_work() const {
    range_type w(-1,-1);
    while (true) {
      const range_type w_new( w.first + 1 , w.second );
      w = atomic_compare_exchange( &m_ranges(0) , w , w_new );
      if ( w.first < w.second ) { // there was work in the queue
        if ( w_new.first == w.first + 1 && w_new.second == w.second ) {
          // we got a work item
          std::int32_t i;
          // the push_work function may have incremented the end counter
          // but not yet written the work index into the queue.
          // wait until the entry is valid.
          while ( -1 == ( i = *((volatile std::int32_t*)(&m_queue( w.first ))) ) );
          return i;
        } // we got a work item
      } else { // there was no work in the queue
#ifdef KOKKOS_DEBUG
        if ( w_new.first == w.first + 1 && w_new.second == w.second ) {
          Kokkos::abort("bug in pop_work");
        }
#endif
        if (w.first == m_total_work) { // all work is done
          return -1;
        } else { // need to wait for more work to be pushed
          // take a guess that one work item will be pushed
          // the key thing is we can't leave (w) alone, because
          // otherwise the next compare_exchange may succeed in
          // popping work from an empty queue
          w.second++;
        }
      } // there was no work in the queue
    } // while (true)
  }

  KOKKOS_INLINE_FUNCTION
  void push_work(std::int32_t i) const {
    range_type w(-1,-1);
    while (true) {
      const range_type w_new( w.first , w.second + 1 );
      // try to increment the end counter
      w = atomic_compare_exchange( &m_ranges(0) , w , w_new );
      // stop trying if the increment was successful
      if ( w.first == w_new.first && w.second + 1 == w_new.second ) break;
    }
    // write the work index into the claimed spot in the queue
    *((volatile std::int32_t*)(&m_queue( w.second ))) = i;
    // push this write out into the memory system
    memory_fence();
  }

  template< class functor_type , class execution_space, class ... policy_args >
  friend class Kokkos::Impl::Experimental::WorkGraphExec;

public:

  WorkGraphPolicy(graph_type arg_graph)
    : m_graph(arg_graph)
    , m_total_work( arg_graph.numRows() )
  {
    setup();
  }

};

}} // namespace Kokkos::Experimental

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
namespace Experimental {

template< class functor_type , class execution_space, class ... policy_args >
class WorkGraphExec
{
 public:

  using self_type = WorkGraphExec< functor_type, execution_space, policy_args ... >;
  using policy_type = Kokkos::Experimental::WorkGraphPolicy< policy_args ... >;
  using member_type = typename policy_type::member_type;
  using memory_space = typename execution_space::memory_space;

 protected:

  const functor_type m_functor;
  const policy_type  m_policy;

 protected:

  KOKKOS_INLINE_FUNCTION
  std::int32_t before_work() const {
    return m_policy.pop_work();
  }

  KOKKOS_INLINE_FUNCTION
  void after_work(std::int32_t i) const {
    /* fence any writes that were done by the work item itself
       (usually writing its result to global memory) */
    memory_fence();
    const std::int32_t begin = m_policy.m_graph.row_map( i );
    const std::int32_t end = m_policy.m_graph.row_map( i + 1 );
    for (std::int32_t j = begin; j < end; ++j) {
      const std::int32_t next = m_policy.m_graph.entries( j );
      const std::int32_t old_count = atomic_fetch_add( &(m_policy.m_counts(next)), -1 );
      if ( old_count == 1 )  m_policy.push_work( next );
    }
  }

  inline
  WorkGraphExec( const functor_type & arg_functor
               , const policy_type  & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
  {
  }
};

}}} // namespace Kokkos::Impl::Experimental

#ifdef KOKKOS_ENABLE_SERIAL
#include "impl/Kokkos_Serial_WorkGraphPolicy.hpp"
#endif

#ifdef KOKKOS_ENABLE_OPENMP
#include "OpenMP/Kokkos_OpenMP_WorkGraphPolicy.hpp"
#endif

#ifdef KOKKOS_ENABLE_CUDA
#include "Cuda/Kokkos_Cuda_WorkGraphPolicy.hpp"
#endif

#ifdef KOKKOS_ENABLE_THREADS
#include "Threads/Kokkos_Threads_WorkGraphPolicy.hpp"
#endif

#endif /* #define KOKKOS_WORKGRAPHPOLICY_HPP */
