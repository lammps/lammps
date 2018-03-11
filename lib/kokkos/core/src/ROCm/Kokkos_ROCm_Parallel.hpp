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

#include <algorithm>
#include <typeinfo>
#include <ROCm/Kokkos_ROCm_Reduce.hpp>
#include <ROCm/Kokkos_ROCm_Scan.hpp>
#include <ROCm/Kokkos_ROCm_Vectorization.hpp>


namespace Kokkos {
namespace Impl {

struct ROCmTeamMember ;

template< class ... Properties >
class TeamPolicyInternal< Kokkos::Experimental::ROCm, Properties ... >: public PolicyTraits<Properties ...> {
private:
  int m_league_size ;
  int m_team_size ;
  int m_vector_length ;
  int m_team_scratch_size[2] ;
  int m_thread_scratch_size[2] ;
  int m_chunk_size ;


public:

  using execution_policy = TeamPolicyInternal ;
  using execution_space  = Kokkos::Experimental::ROCm ;
  typedef PolicyTraits<Properties ... > traits;

  TeamPolicyInternal& operator = (const TeamPolicyInternal& p) {
    m_league_size = p.m_league_size;
    m_team_size = p.m_team_size;
    m_vector_length = p.m_vector_length;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    return *this;
  }

  TeamPolicyInternal()
    : m_league_size( 0 )
    , m_team_size( 0 )
    , m_vector_length( 0 )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 64 )
   {}

  TeamPolicyInternal( const int arg_league_size
            , const int arg_team_size )
    : m_league_size( arg_league_size ),
      m_team_size( arg_team_size )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 64 )
    {}

  TeamPolicyInternal( const int arg_league_size
            , const int arg_team_size
            , const int vector_length_request=1)
    : m_league_size( arg_league_size ),
      m_team_size( arg_team_size ),
      m_vector_length (vector_length_request)
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 64 )
    {}

  TeamPolicyInternal( const int arg_league_size
            , const Kokkos::AUTO_t )
    : m_league_size( arg_league_size ), m_team_size( -1 )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 64 )
    {}

  TeamPolicyInternal( const int arg_league_size
            , const Kokkos::AUTO_t
            , const int vector_length_request)
    : m_league_size( arg_league_size ),
      m_team_size( -1 ),
      m_vector_length (vector_length_request)
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 64 )
    {}

  inline int chunk_size() const { return m_chunk_size ; }

  /** \brief set chunk_size to a discrete value*/
  KOKKOS_INLINE_FUNCTION TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
    p.m_chunk_size = chunk_size_;
    return p;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    return p;
  };

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

// TODO:  evaluate proper team_size_max requirements
  template< class Functor_Type>
  KOKKOS_INLINE_FUNCTION static
  int team_size_max( const Functor_Type & functor)
  {
    typedef typename Kokkos::Impl::FunctorValueTraits<Functor_Type, void>::value_type value_type;
    return team_size_recommended(functor);
    // return std::min(Kokkos::Impl::get_max_tile_size() / sizeof(value_type), Kokkos::Impl::get_max_tile_thread());
  }

  template< class Functor_Type>
  KOKKOS_INLINE_FUNCTION static int team_size_recommended(const Functor_Type & functor)
  { return Kokkos::Impl::get_tile_size<typename Kokkos::Impl::FunctorValueTraits<Functor_Type, void>::value_type>(); }

  template< class Functor_Type >
  KOKKOS_INLINE_FUNCTION static int team_size_recommended(const Functor_Type &functor, const int vector_length)
 {
   int max = team_size_recommended( functor )/vector_length;
   if(max < 1) max = 1;
   return(max);
 }

  template<class F>
  KOKKOS_INLINE_FUNCTION int team_size(const F& f) const { return (m_team_size > 0) ? m_team_size : team_size_recommended(f); }
  KOKKOS_INLINE_FUNCTION int team_size() const { return (m_team_size > 0) ? m_team_size : Impl::get_max_tile_thread(); ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }


  inline int vector_length()   const { return m_vector_length ; }
  inline int scratch_size(int level, int team_size_ = -1) const {
    if(team_size_<0) team_size_ = m_team_size;
    return m_team_scratch_size[level] + team_size_*m_thread_scratch_size[level];
  }
  inline size_t team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }
  inline size_t thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  typedef Impl::ROCmTeamMember member_type;
};

  struct ROCmTeamMember {
    typedef Kokkos::Experimental::ROCm                             execution_space ;
    typedef Kokkos::ScratchMemorySpace<Kokkos::Experimental::ROCm> scratch_memory_space ;

    KOKKOS_INLINE_FUNCTION
    const scratch_memory_space & team_shmem() const 
      { return m_team_shared.set_team_thread_mode(0,1,0); }
    KOKKOS_INLINE_FUNCTION
    const execution_space::scratch_memory_space & team_scratch(const int& level) const
      { return m_team_shared.set_team_thread_mode(level,1,0) ; }
    KOKKOS_INLINE_FUNCTION
    const execution_space::scratch_memory_space & thread_scratch(const int& level) const
      { return m_team_shared.set_team_thread_mode(level,
                                             team_size(),
                                             team_rank()) ; }


    /* Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const { return m_idx.tile[0]; }
    /* Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
    /* Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const { return m_idx.local[0] / m_vector_length; }
    /* Rank of this thread within this thread */
    KOKKOS_INLINE_FUNCTION int vector_rank() const { return m_idx.local[0] % m_vector_length; }
    KOKKOS_INLINE_FUNCTION int lindex() const { return m_idx.local[0]; }
    KOKKOS_INLINE_FUNCTION int gindex() const { return m_idx.global[0]; }
    KOKKOS_INLINE_FUNCTION int tindex() const { return m_idx.tile[0]; }
    KOKKOS_INLINE_FUNCTION int tile_dim() const { return m_idx.tile_dim[0]; }
    KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }
    KOKKOS_INLINE_FUNCTION int vector_length() const { return m_vector_length; }


    KOKKOS_INLINE_FUNCTION
    ROCmTeamMember( const hc::tiled_index< 1 > & arg_idx, int league_size_,int team_size_ )
      : m_league_size( league_size_ )
      , m_team_size( team_size_ )
      , m_team_shared( nullptr, 0 )
      , m_vector_length( 1 )
      , m_idx( arg_idx )
      {}

    KOKKOS_INLINE_FUNCTION
    ROCmTeamMember( const hc::tiled_index< 1 > & arg_idx, int league_size_,int team_size_, char * shared,  std::size_t shsize, std::size_t scratch_size0, char * scratch_ptr, std::size_t scratch_size1, std::size_t vector_length)
      : m_league_size( league_size_ )
      , m_team_size( team_size_ )
      , m_team_shared( shared +  
                          arg_idx.tile[0]*(shsize+scratch_size0), 
                       (shsize+scratch_size0)*league_size_, 
                       scratch_ptr + arg_idx.tile[0]*scratch_size1, 
                       scratch_size1*league_size_)
      , m_vector_length( vector_length )
      , m_idx( arg_idx )
      {}

    KOKKOS_INLINE_FUNCTION
    void team_barrier() const {
      m_idx.barrier.wait();
    }

    template<class ValueType>
    KOKKOS_INLINE_FUNCTION
    void team_broadcast(const ValueType& value, const int& thread_id ) const 
    {
      static_assert(std::is_trivially_default_constructible<ValueType>(), "Only trivial constructible types can be broadcasted");
      tile_static ValueType local_value;
      zero_init(local_value);
      if (this->team_rank() == thread_id) {
        local_value = value;
      }
      this->team_barrier();
      value = local_value;
    }
// Reduce across a team of threads.
//
// Each thread has vector_length elements.
// This reduction is for TeamThreadRange operations, where the range
// is spread across threads.  Effectively, there are vector_length
// independent reduction operations.
// This is different from a reduction across the elements of a thread,
// which reduces every vector element.

    template< class ValueType, class JoinOp >
    KOKKOS_INLINE_FUNCTION
    ValueType team_reduce( const ValueType & value , const JoinOp & op_in) const
    {
      typedef JoinLambdaAdapter<ValueType,JoinOp> JoinOpFunctor ;
      const JoinOpFunctor op(op_in);

      tile_static ValueType buffer[512];
      const auto local = lindex();
      const auto team  = team_rank();
      auto vector_rank = local%m_vector_length;
      auto thread_base = team*m_vector_length;

      const std::size_t size = next_pow_2(m_team_size+1)/2;
#if defined(ROCM15)
      buffer[local] = value;
#else
        // ROCM 1.5 handles address spaces better, previous version didn't
      lds_for(buffer[local], [&](ValueType& x)
      {
          x = value;
      });
#endif
      m_idx.barrier.wait();

      for(std::size_t s = 1; s < size; s *= 2)
      {
          const std::size_t index = 2 * s * team;
          if (index < size)
          {
#if defined(ROCM15)
                op.join(buffer[vector_rank+index*m_vector_length],
                        buffer[vector_rank+(index+s)*m_vector_length]);
#else
              lds_for(buffer[vector_rank+index*m_vector_length], [&](ValueType& x)
              {
                  lds_for(buffer[vector_rank+(index+s)*m_vector_length],
                                [&](ValueType& y)
                  {
                      op.join(x, y);
                  });
              });
#endif
          }
          m_idx.barrier.wait();
      }

      if (local == 0)
      {
          for(int i=size*m_vector_length; i<m_team_size*m_vector_length; i+=m_vector_length)
#if defined(ROCM15)
              op.join(buffer[vector_rank], buffer[vector_rank+i]);
#else
              lds_for(buffer[vector_rank], [&](ValueType& x)
              {
                  lds_for(buffer[vector_rank+i],
                                [&](ValueType& y)
                  {
                      op.join(x, y);
                  });
              });
#endif
      }
      m_idx.barrier.wait();

      return buffer[0];
    }

// Reduce across a team of threads, with a reducer data type
//
// Each thread has vector_length elements.
// This reduction is for TeamThreadRange operations, where the range
// is spread across threads.  Effectively, there are vector_length
// independent reduction operations.
// This is different from a reduction across the elements of a thread,
// which reduces every vector element.

    template< class ReducerType >
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if< is_reducer< ReducerType >::value >::type
    team_reduce( const ReducerType & reducer) const
    {
      typedef typename ReducerType::value_type value_type ;

      tile_static value_type buffer[512];
      const auto local = lindex();
      const auto team  = team_rank();
      auto vector_rank = local%m_vector_length;
      auto thread_base = team*m_vector_length;

      const std::size_t size = next_pow_2(m_team_size+1)/2;
#if defined(ROCM15)
      buffer[local] = reducer.reference();
#else
        // ROCM 1.5 handles address spaces better, previous version didn't
      lds_for(buffer[local], [&](ValueType& x)
      {
          x = value;
      });
#endif
      m_idx.barrier.wait();

      for(std::size_t s = 1; s < size; s *= 2)
      {
          const std::size_t index = 2 * s * team;
          if (index < size)
          {
#if defined(ROCM15)
                reducer.join(buffer[vector_rank+index*m_vector_length],
                        buffer[vector_rank+(index+s)*m_vector_length]);
#else
              lds_for(buffer[vector_rank+index*m_vector_length], [&](ValueType& x)
              {
                  lds_for(buffer[vector_rank+(index+s)*m_vector_length],
                                [&](ValueType& y)
                  {
                      reducer.join(x, y);
                  });
              });
#endif
          }
          m_idx.barrier.wait();
      }

      if (local == 0)
      {
          for(int i=size*m_vector_length; i<m_team_size*m_vector_length; i+=m_vector_length)
#if defined(ROCM15)
              reducer.join(buffer[vector_rank], buffer[vector_rank+i]);
#else
              lds_for(buffer[vector_rank], [&](ValueType& x)
              {
                  lds_for(buffer[vector_rank+i],
                                [&](ValueType& y)
                  {
                      reducer.join(x, y);
                  });
              });
#endif
      }
      m_idx.barrier.wait();
    }

    /** \brief  Intra-team vector reduce 
     *          with intra-team non-deterministic ordering accumulation.
     *
     *  The intra-team accumulation value will, at the end of the
     *  league's parallel execution, be the reduction's total.
     *  Parallel execution ordering of the league's teams is non-deterministic.
     *  As such the base value for each team's vector reduce operation is
     *  similarly non-deterministic.
     */
    template< class ValueType, class JoinOp >
    KOKKOS_INLINE_FUNCTION
    ValueType thread_reduce( const ValueType & value , const JoinOp & op_in) const
    {
      typedef JoinLambdaAdapter<ValueType,JoinOp> JoinOpFunctor ;
      const JoinOpFunctor op(op_in);

      const auto local = m_idx.local[0];
      tile_static ValueType buffer[512];
      const std::size_t size = m_vector_length; //vector length must be power of 2
      auto vector_rank = local%m_vector_length;
      auto thread_base = team_rank()*m_vector_length;
      lds_for(buffer[local], [&](ValueType& x)
      {
          x = value;
      });
      m_idx.barrier.wait();
      for(std::size_t s = 1; s < size; s *= 2)
      {
          const std::size_t index = 2 * s * vector_rank;
          if (index < size)
          {
#if defined(ROCM15)
              op.join(buffer[thread_base+index], buffer[thread_base+index+s]);
#else

              lds_for(buffer[thread_base+index], [&](ValueType& x)
              {
                  lds_for(buffer[thread_base+index+s], [&](ValueType& y)
                  {
                      op.join(x, y);
                  });
              });
#endif
          }
          m_idx.barrier.wait();
      }

      m_idx.barrier.wait();
      return buffer[thread_base];
    }

  template< typename ReducerType >
  KOKKOS_INLINE_FUNCTION static
  typename std::enable_if< is_reducer< ReducerType >::value >::type
  vector_reduce( ReducerType const & reducer )
    {
      #ifdef __HCC_ACCELERATOR__
      if(blockDim_x == 1) return;

      // Intra vector lane shuffle reduction:
      typename ReducerType::value_type tmp ( reducer.reference() );

      for ( int i = blockDim_x ; ( i >>= 1 ) ; ) {
        shfl_down( reducer.reference() , i , blockDim_x );
        if ( (int)threadIdx_x < i ) { reducer.join( tmp , reducer.reference() ); }
      }

      // Broadcast from root lane to all other lanes.
      // Cannot use "butterfly" algorithm to avoid the broadcast
      // because floating point summation is not associative
      // and thus different threads could have different results.

      shfl( reducer.reference() , 0 , blockDim_x );
      #endif
    }



    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
     *          with intra-team non-deterministic ordering accumulation.
     *
     *  The global inter-team accumulation value will, at the end of the
     *  league's parallel execution, be the scan's total.
     *  Parallel execution ordering of the league's teams is non-deterministic.
     *  As such the base value for each team's scan operation is similarly
     *  non-deterministic.
     */
    template< typename Type >
    KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum = nullptr ) const
    {
  #if 0
      const auto local = m_idx.local[0];
      const auto last = m_team_size - 1;
      const auto init = 0;
      tile_static Type buffer[256];

      if (local == last) buffer[0] = init;
      else buffer[local] = value;

      m_idx.barrier.wait();

      for(std::size_t s = 1; s < m_team_size; s *= 2)
      {
          if (local >= s) buffer[local] += buffer[local - s];
          m_idx.barrier.wait();
      }

      if ( global_accum )
      { 
         if(local == last)
         {
            atomic_fetch_add(global_accum, buffer[local] + value);
         }
         m_idx.barrier.wait();
         buffer[local] += *global_accum;
      }
      m_idx.barrier.wait();
      return buffer[local];
#else
      tile_static Type sarray[2][256+1];
      int lid = m_idx.local[0];
      int lp1 = lid+1;

      int toggle = 1;
      int _toggle = 0;
      m_idx.barrier.wait();

      if(lid == 0) 
      {
         sarray[1][0] = 0;
         sarray[0][0] = 0;
      }
      sarray[1][lp1] = value;

      m_idx.barrier.wait();
      for(int stride = 1; stride < m_team_size; stride*=2)
      {
         if(lid >= stride)
         {
            sarray[_toggle][lp1] =
                          sarray[toggle][lp1]+sarray[toggle][lp1-stride];
         }
         else
         {
            sarray[_toggle][lp1] = sarray[toggle][lp1];
         }
         toggle = _toggle;
         _toggle = 1-toggle;
         m_idx.barrier.wait();
      }

      if ( global_accum )
      { 
         if(m_team_size == lp1)
         {
            sarray[toggle][m_team_size] = atomic_fetch_add(global_accum,sarray[toggle][m_team_size]);
         }
         m_idx.barrier.wait();
         sarray[toggle][lid] += sarray[toggle][m_team_size];
      }
      m_idx.barrier.wait();
      return sarray[toggle][lid];
#endif
    }

  private:
    int m_league_size ;
    int m_team_size ;
    const scratch_memory_space  m_team_shared;

  public:
    int m_vector_length;
    hc::tiled_index<1> m_idx;
  };
}
} // namespace Kokkos
#include <ROCm/Kokkos_ROCm_ReduceScan.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , class... Traits >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Traits... >, Kokkos::Experimental::ROCm >
{
private:

  typedef Kokkos::RangePolicy< Traits... > Policy ;

public:

  inline
  ParallelFor( const FunctorType & f
             , const Policy      & policy )
    {


      const auto len = policy.end()-policy.begin();
      const auto offset = policy.begin();
      if(len == 0) return;
// define a lambda to work around a compiler issue.  The compiler does not
// properly dereference f inside the pfe.
auto foo = [=](size_t i){rocm_invoke<typename Policy::work_tag>(f, i);};

#if __hcc_workweek__ > 16600
      hc::parallel_for_each(hc::extent<1>(len) , [=](const hc::index<1> & idx) [[hc]]  [[hc_max_workgroup_dim(1024,1,1)]]
#else
      hc::parallel_for_each(hc::extent<1>(len).tile(256) , [=](const hc::index<1> & idx) [[hc]]
#endif
      {
        if(idx[0]<len)  // workaround for Carrizo (and Fiji?)
          foo(idx[0] + offset);
      }).wait();

    }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

};

//----------------------------------------------------------------------------

template< class F , class... Traits >
class ParallelFor< F
                 , Kokkos::TeamPolicy< Traits... >
                 , Kokkos::Experimental::ROCm >
{
  using Policy = Kokkos::Impl::TeamPolicyInternal< Kokkos::Experimental::ROCm, Traits... >;
  typedef Kokkos::Impl::FunctorValueTraits<F, typename Policy::work_tag> ValueTraits;

public:
  inline
  ParallelFor( const F & f
             , const Policy      & policy )
    {
      const auto league_size  = policy.league_size();
      const auto team_size    = policy.team_size();
      const int vector_length = policy.vector_length();
      const auto total_size   = league_size * team_size * vector_length;
      const int scratch_size0 = policy.scratch_size(0,team_size);
      const int scratch_size1 = policy.scratch_size(1,team_size);

      if(total_size == 0) return;

      const auto shared_size = FunctorTeamShmemSize< F >::value( f , team_size );
      char * scratch = NULL;
      char * shared = (char *)rocm_device_allocate(shared_size * league_size +
                                                   scratch_size0*league_size);
      if(0<scratch_size1)
        scratch = (char *)rocm_device_allocate(scratch_size1*league_size);

      hc::extent< 1 > flat_extent( total_size );

      hc::tiled_extent< 1 > team_extent = flat_extent.tile(team_size*vector_length);
      hc::parallel_for_each( team_extent , [=](hc::tiled_index<1> idx) [[hc]]
      {
        rocm_invoke<typename Policy::work_tag>(f, typename Policy::member_type(idx, league_size, team_size, shared, shared_size, scratch_size0, scratch, scratch_size1,vector_length));
      }).wait();

      if(0<scratch_size1)
        rocm_device_free(scratch);
      rocm_device_free(shared);
    }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

};


//----------------------------------------------------------------------------

template< class FunctorType , class ReducerType, class... Traits >
class ParallelReduce<
  FunctorType , Kokkos::RangePolicy< Traits... >, ReducerType, Kokkos::Experimental::ROCm >
{
public:

  typedef Kokkos::RangePolicy< Traits... > Policy ;

  // TODO: Use generic lambdas instead
  struct invoke_fn
  {
    template<class F, class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(std::size_t size, F&& f, hc::tiled_index<1> idx, tile_desc td, Ts&&... xs) const
    {
      auto global = idx.global[0];
      if (global < size) f(idx.global[0], static_cast<Ts&&>(xs)...);
    }
  };

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType  & f,
                  const Policy       & policy,
                  const ViewType & result_view,
                  typename std::enable_if<
                               Kokkos::is_view< ViewType >::value &&
                              !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    {
      typedef typename Policy::work_tag Tag;
      typedef Kokkos::Impl::FunctorValueTraits< FunctorType , Tag > ValueTraits;
      typedef Kokkos::Impl::FunctorValueInit< FunctorType , Tag > ValueInit;
      typedef typename ValueTraits::reference_type reference_type;

      const auto total_size = policy.end() - policy.begin();

      if(total_size==0) {
        if (result_view.data()) {
           ValueInit::init( f , result_view.data() );
        }
        return;
      }

      Kokkos::Impl::reduce_enqueue< Tag >
        ( total_size 
        , f
        , InvalidType{}
        , rocm_capture(invoke_fn{}, total_size)
        , result_view.data()
        , result_view.extent(0)
        );
    }

  inline
  ParallelReduce( const FunctorType & f,
                  Policy       policy,
                  const ReducerType& reducer )
  {
      typedef typename Policy::work_tag Tag;

      typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value,                                   FunctorType, ReducerType> ReducerConditional;
      typedef typename ReducerConditional::type ReducerTypeFwd;
      typedef Kokkos::Impl::FunctorValueTraits< FunctorType , Tag > ValueTraits;
      typedef Kokkos::Impl::FunctorValueInit< ReducerType, Tag > ValueInit ;

      typedef typename ValueTraits::reference_type reference_type;

      const auto total_size = policy.end() - policy.begin();

      if(total_size==0) {
        if (reducer.view().data()) {
           ValueInit::init( ReducerConditional::select(f,reducer), 
                            reducer.view().data() );
        }
        return;
      }

      Kokkos::Impl::reduce_enqueue< Tag >
        ( total_size 
        , f
        , reducer
        , rocm_capture(invoke_fn{}, total_size)
        , reducer.view().data()
        , reducer.view().extent(0)
        );
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

};

template< class FunctorType, class ReducerType, class... Traits >
class ParallelReduce<
   FunctorType , Kokkos::TeamPolicy< Traits... >, ReducerType, Kokkos::Experimental::ROCm >
{
  using Policy = Kokkos::Impl::TeamPolicyInternal< Kokkos::Experimental::ROCm, Traits... >;
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType, typename Policy::work_tag> ValueTraits;

public:

  struct invoke_fn
  {
    template<class Create, class F, class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(Create&& create, F&& f, hc::tiled_index<1> idx, tile_desc td, Ts&&... xs) const
    {
      f(create(idx, td), static_cast<Ts&&>(xs)...);
    }
  };

  template< class ViewType >
  inline
  ParallelReduce( const FunctorType  & f,
                  const Policy       & policy,
                  const ViewType     & result_view,
                typename std::enable_if<
                  Kokkos::is_view< ViewType >::value &&
                  !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    {
      const int league_size = policy.league_size();
      const int team_size = policy.team_size(f);
      const int vector_length = policy.vector_length();
      const int scratch_size0 = policy.scratch_size(0,team_size);
      const int scratch_size1 = policy.scratch_size(1,team_size);
      const int total_size = league_size * team_size ;

      if(total_size == 0) return;

      const int reduce_size = ValueTraits::value_size( f );
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( f , team_size );

      char * shared;
      char * scratch = NULL;

      shared = (char *)rocm_device_allocate(league_size *
                             (shared_size + scratch_size0));
      if(0<scratch_size1)
        scratch = (char *)rocm_device_allocate(scratch_size1 * league_size);

      auto create_team_member = [=](hc::tiled_index<1> idx, tile_desc td) 
      { 

        return typename Policy::member_type(idx, league_size, td.team_size, 
                                          shared, shared_size, scratch_size0,
                                          scratch, scratch_size1, 
                                          vector_length); 
      };

      Kokkos::Impl::reduce_enqueue< typename Policy::work_tag >
      ( total_size*vector_length
        , f
        , InvalidType{}
        , rocm_capture(invoke_fn{}, create_team_member)
        , result_view.ptr_on_device()
        , result_view.dimension_0()
        , team_size 
        , vector_length 
        , shared_size
      );

      if(0<scratch_size1)
        rocm_device_free(scratch);
      rocm_device_free(shared);
    }

  inline
  ParallelReduce( const FunctorType & f,
                  Policy       policy,
                  const ReducerType& reducer )
  {
    const int league_size = policy.league_size();
      const int team_size = policy.team_size(f);
      const int vector_length = policy.vector_length();
      const int total_size = league_size * team_size;

      if(total_size == 0) return;

      const int reduce_size = ValueTraits::value_size( f );
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( f , team_size );
      const int scratch_size0 = policy.scratch_size(0,team_size);
      const int scratch_size1 = policy.scratch_size(1,team_size);

      char * shared;
      char * scratch = NULL;
      shared = (char *)rocm_device_allocate((shared_size + scratch_size0) *
                                            league_size);
      if(0<scratch_size1)
        scratch = (char *)rocm_device_allocate(scratch_size1 * league_size);

      auto create_team_member = [=](hc::tiled_index<1> idx, tile_desc td) 
      { 
        return typename Policy::member_type(idx, league_size, td.tile_size, shared, shared_size, scratch_size0, scratch, scratch_size1, vector_length); 
      };

      Kokkos::Impl::reduce_enqueue< typename Policy::work_tag >
      ( league_size
        , f
        , reducer
        , rocm_capture(invoke_fn{}, create_team_member)
        , reducer.view().data()
        , reducer.view().extent(0),team_size,vector_length
        , shared_size
     );

      if(0<scratch_size1)
        rocm_device_free(scratch);
      rocm_device_free(shared);
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

};


template< class FunctorType , class... Traits >
class ParallelScan< FunctorType , Kokkos::RangePolicy< Traits... >, Kokkos::Experimental::ROCm >
{
private:

  typedef Kokkos::RangePolicy< Traits... > Policy;
  typedef typename Policy::work_tag Tag;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, Tag>  ValueTraits;

public:

  //----------------------------------------

  inline
  ParallelScan( const FunctorType & f
              , const Policy      & policy )
  {
    const auto len = policy.end()-policy.begin();


    if(len==0) return;

    scan_enqueue<Tag>(len, f, [](hc::tiled_index<1> idx, int, int) { return idx.global[0]; });
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

  //----------------------------------------
};

template< class FunctorType , class... Traits>
class ParallelScan< FunctorType , Kokkos::TeamPolicy< Traits... >, Kokkos::Experimental::ROCm >
{
private:

  using Policy = Kokkos::Impl::TeamPolicyInternal< Kokkos::Experimental::ROCm, Traits... >;
  typedef typename Policy::work_tag Tag;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, Tag>  ValueTraits;

public:

  //----------------------------------------

  inline
  ParallelScan( const FunctorType & f
              , const Policy      & policy )
  {
    const auto league_size = policy.league_size();
    const auto team_size = policy.team_size(f);
    const auto len  = league_size * team_size;
      
    if(len == 0) return;

    scan_enqueue<Tag>(len, f, [&](hc::tiled_index<1> idx, int n_teams, int n_leagues) { return typename Policy::member_type(idx,n_leagues,n_teams); });
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

  //----------------------------------------
};

}
}

namespace Kokkos {
namespace Impl {
  template<typename iType>
  struct TeamThreadRangeBoundariesStruct<iType,ROCmTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;
    const ROCmTeamMember& thread;

#if defined( __HCC_ACCELERATOR__ )
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const ROCmTeamMember& thread_, const iType& count):
      start( thread_.team_rank() ),
      end( count ),
      increment( thread_.team_size() ),
      thread(thread_)
    {}
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const ROCmTeamMember& thread_,  const iType& begin_, const iType& end_):
      start( begin_ + thread_.team_rank() ),
      end( end_ ),
      increment( thread_.team_size() ),
      thread(thread_)
    {}
#else
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const ROCmTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 ),
      thread(thread_)
    {}
    KOKKOS_INLINE_FUNCTION
    TeamThreadRangeBoundariesStruct (const ROCmTeamMember& thread_,  const iType& begin_, const iType& end_):
      start( begin_ ),
      end( end_ ),
      increment( 1 ),
      thread(thread_)
    {}
#endif
  };
  template<typename iType>
  struct ThreadVectorRangeBoundariesStruct<iType,ROCmTeamMember> {
    typedef iType index_type;
    const iType start;
    const iType end;
    const iType increment;
    const ROCmTeamMember& thread;

#if defined( __HCC_ACCELERATOR__ )
    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const ROCmTeamMember& thread_, const iType& count):
      start( thread_.lindex()%thread_.vector_length() ),
      end( count ),
      increment( thread_.vector_length() ),
      thread(thread_)
    {}

//    KOKKOS_INLINE_FUNCTION
//    ThreadVectorRangeBoundariesStruct (const iType& count):
//      start( 0 ),
//      end( count ),
//      increment( 1 )
//    {}
#else
    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const ROCmTeamMember& thread_, const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 ),
      thread(thread_)
    {}
    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const iType& count):
      start( 0 ),
      end( count ),
      increment( 1 )
    {}
#endif
  };

}
}

namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>
  TeamThreadRange(const Impl::ROCmTeamMember& thread, const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>(thread,count);
}

template<typename iType1,typename iType2>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<typename std::common_type< iType1, iType2 >::type,Impl::ROCmTeamMember>
  TeamThreadRange(const Impl::ROCmTeamMember& thread, const iType1& begin, const iType2& end) {
  typedef typename std::common_type< iType1, iType2 >::type iType;
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>(thread,begin,end);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >
  ThreadVectorRange(const Impl::ROCmTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >(thread,count);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::ROCmTeamMember> PerTeam(const Impl::ROCmTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::ROCmTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::ROCmTeamMember> PerThread(const Impl::ROCmTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::ROCmTeamMember>(thread);
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::ROCmTeamMember>& single_struct, const FunctorType& lambda) {
  if(single_struct.team_member.vector_rank()==0) lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::ROCmTeamMember>& single_struct, const FunctorType& lambda) {
  if((single_struct.team_member.lindex()==0)) lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::ROCmTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
#if defined(ROCM15)
  // 1.5 needs this more proper restriction on which work units run
  if( single_struct.team_member.vector_rank()==0) lambda(val);
  val = shfl(val,0,single_struct.team_member.vector_length());
#else
  // but older compilers are fine with this (TestTeamVector::Test< Kokkos::Experimental::ROCm >(4))
  lambda(val);
#endif
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::ROCmTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
  if(single_struct.team_member.lindex()==0) lambda(val);
  single_struct.team_member.team_broadcast(val,0);
}

}

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread thread range parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }
  result = loop_boundaries.thread.team_reduce(result,
                                              Impl::JoinAdd<ValueType>());
//  Impl::rocm_intra_workgroup_reduction( loop_boundaries.thread, result,
//               Impl::JoinAdd<ValueType>());
//  Impl::rocm_inter_workgroup_reduction( loop_boundaries.thread, result,
//               Impl::JoinAdd<ValueType>());
}

/** \brief  Inter-thread thread range parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ReducerType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>& loop_boundaries,
                     const Lambda & lambda, ReducerType const & reducer) {
  reducer.init( reducer.reference() );

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,reducer.reference());
  }
  loop_boundaries.thread.team_reduce(reducer);
}

/** \brief  Intra-thread thread range parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::ROCmTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& result) {

#if defined(ROCM15)
  ValueType tmp = result;
  //  Simpler code works with ROCM1.5
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,tmp);
  }
  result = loop_boundaries.thread.team_reduce(tmp,join);
#else
  // this workaround freezes up with ROCM1.5, but needed for earlier compilers
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }
  result = loop_boundaries.thread.team_reduce(result,join);
#endif
//  Impl::rocm_intra_workgroup_reduction( loop_boundaries.thread, result,join);
//  Impl::rocm_inter_workgroup_reduction( loop_boundaries.thread, result,join);
}

} //namespace Kokkos


namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
    loop_boundaries, const Lambda& lambda) {

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
  result = ValueType();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }
  result = loop_boundaries.thread.thread_reduce(result,Impl::JoinAdd<ValueType>());
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ValueType& result) {

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);  
    loop_boundaries.thread.team_barrier();
  }
  result = loop_boundaries.thread.thread_reduce(result,join);
}


/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ReducerType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
      loop_boundaries, const Lambda & lambda, ReducerType const & reducer) {
  reducer.init( reducer.reference() );

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,reducer.reference());
  }
  loop_boundaries.thread.vector_reduce(reducer);
}
/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ReducerType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ReducerType const & reducer) {

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,reducer.reference());  
    loop_boundaries.thread.team_barrier();
  }
  reducer.reference() = loop_boundaries.thread.thread_reduce(reducer.reference(),join);
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes lambda(iType i, ValueType & val, bool final)
 *          for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan operation is performed.
 * Depending on the target execution space the operator might be called twice: once with final=false
 * and once with final=true. When final==true val contains the prefix sum value. The contribution of this
 * "i" needs to be added to val no matter whether final==true or not. In a serial execution
 * (i.e. team_size==1) the operator is only called once with final==true. Scan_val will be set
 * to the final sum value over all vector lanes.
 * This functionality requires C++11 support.*/
template< typename iType, class FunctorType >
KOKKOS_INLINE_FUNCTION
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::ROCmTeamMember >&
      loop_boundaries, const FunctorType & lambda) {

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void > ValueTraits ;
  typedef typename ValueTraits::value_type value_type ;

  value_type scan_val = value_type();
#if (__ROCM_ARCH__ >= 800)
// adopt the cuda vector shuffle method
  const int VectorLength = loop_boundaries.increment;
  int lid = loop_boundaries.thread.lindex();
  int vector_rank = lid%VectorLength;

  iType loop_bound = ((loop_boundaries.end+VectorLength-1)/VectorLength) * VectorLength;
  value_type val ;
  for(int _i = vector_rank; _i < loop_bound; _i += VectorLength) {
    val = value_type();
    if(_i<loop_boundaries.end)
      lambda(_i , val , false);

    value_type tmp = val;
    value_type result_i;

    if(vector_rank == 0)
      result_i = tmp;
    if (VectorLength > 1) {
      const value_type tmp2 = shfl_up(tmp, 1,VectorLength);
      if(vector_rank > 0)
        tmp+=tmp2;
    }
    if(vector_rank == 1)
      result_i = tmp;
    if (VectorLength > 3) {
      const value_type tmp2 = shfl_up(tmp, 2,VectorLength);
      if(vector_rank > 1)
        tmp+=tmp2;
    }
    if ((vector_rank >= 2) &&
        (vector_rank < 4))
      result_i = tmp;
    if (VectorLength > 7) {
      const value_type tmp2 = shfl_up(tmp, 4,VectorLength);
      if(vector_rank > 3)
        tmp+=tmp2;
    }
    if ((vector_rank >= 4) &&
        (vector_rank < 8))
      result_i = tmp;
    if (VectorLength > 15) {
      const value_type tmp2 = shfl_up(tmp, 8,VectorLength);
      if(vector_rank > 7)
        tmp+=tmp2;
    }
    if ((vector_rank >= 8) &&
        (vector_rank < 16))
      result_i = tmp;
    if (VectorLength > 31) {
      const value_type tmp2 = shfl_up(tmp, 16,VectorLength);
      if(vector_rank > 15)
        tmp+=tmp2;
    }
    if ((vector_rank >=16) &&
        (vector_rank < 32))
      result_i = tmp;
    if (VectorLength > 63) {
      const value_type tmp2 = shfl_up(tmp, 32,VectorLength);
      if(vector_rank > 31)
        tmp+=tmp2;
    }

    if (vector_rank >= 32)
      result_i = tmp;

    val = scan_val + result_i - val;
    scan_val += shfl(tmp,VectorLength-1,VectorLength);
    if(_i<loop_boundaries.end)
      lambda(_i , val , true);
  }
#else
// for kaveri, call the LDS based thread_scan routine
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,scan_val,true);
  }
  scan_val = loop_boundaries.thread.team_scan(scan_val);

#endif
}

} // namespace Kokkos

