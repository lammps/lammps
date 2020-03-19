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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDA_PARALLEL_HPP
#define KOKKOS_CUDA_PARALLEL_HPP

#include <Kokkos_Macros.hpp>
#if defined( __CUDACC__ ) && defined( KOKKOS_ENABLE_CUDA )

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdint>

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <Cuda/Kokkos_Cuda_Locks.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>
#endif

#include <KokkosExp_MDRangePolicy.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

extern bool show_warnings() noexcept;

namespace Impl {

template< class ... Properties >
class TeamPolicyInternal< Kokkos::Cuda , Properties ... >: public PolicyTraits<Properties ... >
{
public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicyInternal      execution_policy ;

  typedef PolicyTraits<Properties ... > traits;

  template< class ExecSpace, class ... OtherProperties >
  friend class TeamPolicyInternal;

private:

  enum { MAX_WARP = 8 };

  typename traits::execution_space m_space;
  int m_league_size ;
  int m_team_size ;
  int m_vector_length ;
  int m_team_scratch_size[2] ;
  int m_thread_scratch_size[2] ;
  int m_chunk_size;

public:

  //! Execution space of this execution policy
  typedef Kokkos::Cuda  execution_space ;

  template<class ... OtherProperties>
  TeamPolicyInternal( const TeamPolicyInternal<OtherProperties...>& p ) {
    m_league_size = p.m_league_size;
    m_team_size = p.m_team_size;
    m_vector_length = p.m_vector_length;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    m_space = p.m_space;
  }

  TeamPolicyInternal& operator = (const TeamPolicyInternal& p) {
    m_league_size = p.m_league_size;
    m_team_size = p.m_team_size;
    m_vector_length = p.m_vector_length;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    m_space = p.m_space;
    return *this;
  }

  //----------------------------------------

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  template< class FunctorType >
  static inline
  int team_size_max( const FunctorType & functor )
    {
      int n = MAX_WARP * Impl::CudaTraits::WarpSize ;

      for ( ; n ; n >>= 1 ) {
        const int shmem_size =
          /* for global reduce */ Impl::cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,typename traits::work_tag>( functor , n )
          /* for team   reduce */ + ( n + 2 ) * sizeof(double)
          /* for team   shared */ + Impl::FunctorTeamShmemSize< FunctorType >::value( functor , n );

        if ( shmem_size < typename traits::execution_space().impl_internal_space_instance()->m_maxShmemPerBlock ) break ;
      }

      return n ;
    }
#endif

  template<class FunctorType>
  int team_size_max( const FunctorType& f, const ParallelForTag& ) const {
    typedef Impl::ParallelFor< FunctorType , TeamPolicy<Properties...> > closure_type;
    cudaFuncAttributes attr = CudaParallelLaunch< closure_type, typename traits::launch_bounds >::
        get_cuda_func_attributes();
    int block_size = Kokkos::Impl::cuda_get_max_block_size< FunctorType, typename traits::launch_bounds >( 
        space().impl_internal_space_instance(),attr,f ,(size_t) vector_length(),
        (size_t) team_scratch_size(0) + 2*sizeof(double), (size_t) thread_scratch_size(0) + sizeof(double) );
    return block_size/vector_length();
  }

  template<class FunctorType>
  int team_size_max( const FunctorType& f, const ParallelReduceTag& ) const {
    typedef Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,TeamPolicyInternal,FunctorType> functor_analysis_type;
    typedef typename Impl::ParallelReduceReturnValue<void,typename functor_analysis_type::value_type,FunctorType>::reducer_type reducer_type;
    typedef Impl::ParallelReduce< FunctorType , TeamPolicy<Properties...>, reducer_type > closure_type;
    typedef Impl::FunctorValueTraits< FunctorType , typename traits::work_tag > functor_value_traits;

    cudaFuncAttributes attr = CudaParallelLaunch< closure_type, typename traits::launch_bounds >::
        get_cuda_func_attributes();
    int block_size = Kokkos::Impl::cuda_get_max_block_size< FunctorType, typename traits::launch_bounds >( 
        space().impl_internal_space_instance(),attr,f ,(size_t) vector_length(),
        (size_t) team_scratch_size(0) + 2*sizeof(double), (size_t) thread_scratch_size(0) + sizeof(double) +
                                                          ((functor_value_traits::StaticValueSize!=0)?0:functor_value_traits::value_size( f )));

    // Currently we require Power-of-2 team size for reductions.
    int p2 = 1;
    while(p2<=block_size) p2*=2;
    p2/=2;
    return p2/vector_length();
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  template< class FunctorType >
  static int team_size_recommended( const FunctorType & functor )
    { return team_size_max( functor ); }

  template< class FunctorType >
  static int team_size_recommended( const FunctorType & functor , const int vector_length)
    {
      int max = team_size_max( functor )/vector_length;
      if(max<1) max = 1;
      return max;
    }
#endif

  template<class FunctorType>
  int team_size_recommended( const FunctorType& f, const ParallelForTag& ) const {
    typedef Impl::ParallelFor< FunctorType , TeamPolicy<Properties...> > closure_type;
    cudaFuncAttributes attr = CudaParallelLaunch< closure_type, typename traits::launch_bounds >::
        get_cuda_func_attributes();
    const int block_size = Kokkos::Impl::cuda_get_opt_block_size< FunctorType, typename traits::launch_bounds>(
        space().impl_internal_space_instance(),
        attr, f , (size_t) vector_length(),
        (size_t) team_scratch_size(0) + 2*sizeof(double), (size_t) thread_scratch_size(0) + sizeof(double));
    return block_size/vector_length();
  }

  template<class FunctorType>
  int team_size_recommended( const FunctorType& f, const ParallelReduceTag& ) const {
    typedef Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,TeamPolicyInternal,FunctorType> functor_analysis_type;
    typedef typename Impl::ParallelReduceReturnValue<void,typename functor_analysis_type::value_type,FunctorType>::reducer_type reducer_type;
    typedef Impl::ParallelReduce< FunctorType , TeamPolicy<Properties...>, reducer_type > closure_type;
    typedef Impl::FunctorValueTraits< FunctorType , typename traits::work_tag > functor_value_traits;

    cudaFuncAttributes attr = CudaParallelLaunch< closure_type, typename traits::launch_bounds >::
        get_cuda_func_attributes();
    const int block_size = Kokkos::Impl::cuda_get_opt_block_size< FunctorType, typename traits::launch_bounds>(
        space().impl_internal_space_instance(),
        attr, f , (size_t) vector_length(),
        (size_t) team_scratch_size(0) + 2*sizeof(double), (size_t) thread_scratch_size(0) + sizeof(double) +
                                                          ((functor_value_traits::StaticValueSize!=0)?0:functor_value_traits::value_size( f )));
    // Currently we require Power-of-2 team size for reductions.
    int p2 = 1;
    while(p2<=block_size) p2*=2;
    p2/=2;
    return p2/vector_length();
  }


  inline static
  int vector_length_max()
    { return Impl::CudaTraits::WarpSize; }

  inline static
  int verify_requested_vector_length( int requested_vector_length ) {
      int test_vector_length = std::min( requested_vector_length, vector_length_max() );

      // Allow only power-of-two vector_length
      if ( !(is_integral_power_of_two( test_vector_length ) ) ) {
         int test_pow2 = 1;
         for (int i = 0; i < 5; i++) {
            test_pow2 = test_pow2 << 1;
            if (test_pow2 > test_vector_length) {
               break;
            }
         }
         test_vector_length = test_pow2 >> 1;
      }

      return test_vector_length;
  }

  inline static
  int scratch_size_max(int level)
    { return (level==0?
        1024*40:             // 48kB is the max for CUDA, but we need some for team_member.reduce etc.
        20*1024*1024);   // arbitrarily setting this to 20MB, for a Volta V100 that would give us about 3.2GB for 2 teams per SM
    }

  //----------------------------------------

  inline int vector_length()   const { return m_vector_length ; }
  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }
  inline int scratch_size(int level, int team_size_ = -1) const {
    if(team_size_<0) team_size_ = m_team_size;
    return m_team_scratch_size[level] + team_size_*m_thread_scratch_size[level];
  }
  inline int team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }
  inline int thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  inline typename traits::execution_space space() const {
    return m_space;
  }

  TeamPolicyInternal()
    : m_space(typename traits::execution_space())
    , m_league_size( 0 )
    , m_team_size( -1 )
    , m_vector_length( 0 )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 32 )
   {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal( const execution_space space_
            , int league_size_
            , int team_size_request
            , int vector_length_request = 1 )
    : m_space( space_ )
    , m_league_size( league_size_ )
    , m_team_size( team_size_request )
    , m_vector_length( verify_requested_vector_length(vector_length_request) )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 32 )
    {
      // Make sure league size is permissible
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");

      // Make sure total block size is permissible
      if ( m_team_size * m_vector_length > 1024 ) {
        Impl::throw_runtime_exception(std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. Team size x vector length must be smaller than 1024."));
      }
    }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal( const execution_space space_
            , int league_size_
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1 )
    : m_space( space_ )
    , m_league_size( league_size_ )
    , m_team_size( -1 )
    , m_vector_length( verify_requested_vector_length(vector_length_request) )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 32 )
    {
      // Make sure league size is permissible
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  TeamPolicyInternal( int league_size_
            , int team_size_request
            , int vector_length_request = 1 )
    : m_space( typename traits::execution_space() )
    , m_league_size( league_size_ )
    , m_team_size( team_size_request )
    , m_vector_length ( verify_requested_vector_length(vector_length_request) )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 32 )
    {
      // Make sure league size is permissible
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");

      // Make sure total block size is permissible
      if ( m_team_size * m_vector_length > 1024 ) {
        Impl::throw_runtime_exception(std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. Team size x vector length must be smaller than 1024."));
      }
    }

  TeamPolicyInternal( int league_size_
            , const Kokkos::AUTO_t & /* team_size_request */
            , int vector_length_request = 1 )
    : m_space( typename traits::execution_space() )
    , m_league_size( league_size_ )
    , m_team_size( -1 )
    , m_vector_length ( verify_requested_vector_length(vector_length_request) )
    , m_team_scratch_size {0,0}
    , m_thread_scratch_size {0,0}
    , m_chunk_size ( 32 )
    {
      // Make sure league size is permissible
      if(league_size_ >= int(Impl::cuda_internal_maximum_grid_count()))
        Impl::throw_runtime_exception( "Requested too large league_size for TeamPolicy on Cuda execution space.");
    }

  inline int chunk_size() const { return m_chunk_size ; }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
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
#else
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal& set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) {
    m_team_scratch_size[level] = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
#endif

  typedef Kokkos::Impl::CudaTeamMember member_type ;

protected:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal internal_set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) {
    m_team_scratch_size[level] = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
#endif
};

} // namspace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Traits ... >
                 , Kokkos::Cuda
                 >
{
public:
  typedef Kokkos::RangePolicy< Traits ... > Policy;
private:

  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::launch_bounds LaunchBounds ;

  const FunctorType  m_functor ;
  const Policy       m_policy ;

  ParallelFor() = delete ;
  ParallelFor & operator = ( const ParallelFor & ) = delete ;

  template< class TagType >
  inline __device__
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member i ) const
    { m_functor( i ); }

  template< class TagType >
  inline __device__
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member i ) const
    { m_functor( TagType() , i ); }

public:

  typedef FunctorType functor_type ;

  inline
  __device__
  void operator()(void) const
    {
      const Member work_stride = blockDim.y * gridDim.x ;
      const Member work_end    = m_policy.end();

      for ( Member
              iwork =  m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x ;
              iwork <  work_end ;
              iwork = iwork < work_end - work_stride ? iwork + work_stride : work_end) {
        this-> template exec_range< WorkTag >( iwork );
      }
    }

  inline
  void execute() const
    {
      const typename Policy::index_type nwork = m_policy.end() - m_policy.begin();

      cudaFuncAttributes attr = CudaParallelLaunch< ParallelFor, LaunchBounds >::
          get_cuda_func_attributes();
      const int block_size = Kokkos::Impl::cuda_get_opt_block_size< FunctorType, LaunchBounds>(
          m_policy.space().impl_internal_space_instance(),
          attr, m_functor , 1, 0 , 0 );
      const dim3 block(  1 , block_size , 1);
      const dim3 grid( std::min( typename Policy::index_type(( nwork + block.y - 1 ) / block.y) ,
                       typename Policy::index_type(cuda_internal_maximum_grid_count()) ) , 1 , 1);

      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_policy.space().impl_internal_space_instance() , false );
    }

  ParallelFor( const FunctorType  & arg_functor ,
               const Policy       & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    { }

};


// MDRangePolicy impl
template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType
                 , Kokkos::MDRangePolicy< Traits ... >
                 , Kokkos::Cuda
                 >
{
public:
  typedef Kokkos::MDRangePolicy< Traits ...  > Policy ;
private:
  using RP = Policy;
  typedef typename Policy::array_index_type array_index_type;
  typedef typename Policy::index_type index_type;
  typedef typename Policy::launch_bounds LaunchBounds;


  const FunctorType m_functor ;
  const Policy      m_rp ;

public:

  inline
  __device__
  void operator()(void) const
    {
      Kokkos::Impl::Refactor::DeviceIterateTile<Policy::rank,Policy,FunctorType,typename Policy::work_tag>(m_rp,m_functor).exec_range();
    }


  inline
  void execute() const
  {
    if(m_rp.m_num_tiles==0) return;
    const array_index_type maxblocks = static_cast<array_index_type>(m_rp.space().impl_internal_space_instance()->m_maxBlock);
    if ( RP::rank == 2 )
    {
      const dim3 block( m_rp.m_tile[0] , m_rp.m_tile[1] , 1);
      const dim3 grid(
            std::min( ( m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1 ) / block.x , maxblocks )
          , std::min( ( m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1 ) / block.y , maxblocks )
          , 1
          );
      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_rp.space().impl_internal_space_instance() , false );
    }
    else if ( RP::rank == 3 )
    {
      const dim3 block( m_rp.m_tile[0] , m_rp.m_tile[1] , m_rp.m_tile[2] );
      const dim3 grid(
          std::min( ( m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1 ) / block.x , maxblocks )
        , std::min( ( m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1 ) / block.y , maxblocks )
        , std::min( ( m_rp.m_upper[2] - m_rp.m_lower[2] + block.z - 1 ) / block.z , maxblocks )
        );
      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_rp.space().impl_internal_space_instance() , false );
    }
    else if ( RP::rank == 4 )
    {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to threadIdx.z
      const dim3 block( m_rp.m_tile[0]*m_rp.m_tile[1] , m_rp.m_tile[2] , m_rp.m_tile[3] );
      const dim3 grid(
          std::min( static_cast<index_type>( m_rp.m_tile_end[0] * m_rp.m_tile_end[1] )
                  , static_cast<index_type>(maxblocks) )
        , std::min( ( m_rp.m_upper[2] - m_rp.m_lower[2] + block.y - 1 ) / block.y , maxblocks )
        , std::min( ( m_rp.m_upper[3] - m_rp.m_lower[3] + block.z - 1 ) / block.z , maxblocks )
        );
      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_rp.space().impl_internal_space_instance() , false );
    }
    else if ( RP::rank == 5 )
    {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4 to threadIdx.z
      const dim3 block( m_rp.m_tile[0]*m_rp.m_tile[1] , m_rp.m_tile[2]*m_rp.m_tile[3] , m_rp.m_tile[4] );
      const dim3 grid(
          std::min( static_cast<index_type>( m_rp.m_tile_end[0] * m_rp.m_tile_end[1] )
                  , static_cast<index_type>(maxblocks) )
        , std::min( static_cast<index_type>( m_rp.m_tile_end[2] * m_rp.m_tile_end[3] )
                  , static_cast<index_type>(maxblocks) )
        , std::min( ( m_rp.m_upper[4] - m_rp.m_lower[4] + block.z - 1 ) / block.z , maxblocks )
        );
      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_rp.space().impl_internal_space_instance() , false );
    }
    else if ( RP::rank == 6 )
    {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4,id5 to threadIdx.z
      const dim3 block( m_rp.m_tile[0]*m_rp.m_tile[1] , m_rp.m_tile[2]*m_rp.m_tile[3] , m_rp.m_tile[4]*m_rp.m_tile[5] );
      const dim3 grid(
          std::min( static_cast<index_type>( m_rp.m_tile_end[0] * m_rp.m_tile_end[1] )
                  , static_cast<index_type>(maxblocks) )
        ,  std::min( static_cast<index_type>( m_rp.m_tile_end[2] * m_rp.m_tile_end[3] )
                  , static_cast<index_type>(maxblocks) )
        , std::min( static_cast<index_type>( m_rp.m_tile_end[4] * m_rp.m_tile_end[5] )
                  , static_cast<index_type>(maxblocks) )
        );
      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this , grid , block , 0 , m_rp.space().impl_internal_space_instance() , false );
    }
    else
    {
      printf("Kokkos::MDRange Error: Exceeded rank bounds with Cuda\n");
      Kokkos::abort("Aborting");
    }

  } //end execute

//  inline
  ParallelFor( const FunctorType & arg_functor
             , Policy arg_policy )
    : m_functor( arg_functor )
    , m_rp(  arg_policy )
    {}
};


template< class FunctorType , class ... Properties >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Properties ... >
                 , Kokkos::Cuda
                 >
{
public:
  typedef TeamPolicyInternal< Kokkos::Cuda , Properties ... >   Policy ;
private:

  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::launch_bounds  LaunchBounds ;

public:

  typedef FunctorType      functor_type ;
  typedef Cuda::size_type  size_type ;

private:

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType  m_functor ;
  const Policy       m_policy ;
  const size_type    m_league_size ;
  int    m_team_size ;
  const size_type    m_vector_size ;
  int m_shmem_begin ;
  int m_shmem_size ;
  void*              m_scratch_ptr[2] ;
  int m_scratch_size[2] ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const Member & member ) const
    { m_functor( member ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const Member & member ) const
    { m_functor( TagType() , member ); }

public:

  __device__ inline
  void operator()(void) const
  {
    // Iterate this block through the league
    int64_t threadid = 0;
    if ( m_scratch_size[1]>0 ) {
      __shared__ int64_t base_thread_id;
      if (threadIdx.x==0 && threadIdx.y==0 ) {
        threadid = (blockIdx.x*blockDim.z + threadIdx.z) %
          (Kokkos::Impl::g_device_cuda_lock_arrays.n / (blockDim.x * blockDim.y));
        threadid *= blockDim.x * blockDim.y;
        int done = 0;
        while (!done) {
          done = (0 == atomicCAS(&Kokkos::Impl::g_device_cuda_lock_arrays.scratch[threadid],0,1));
          if(!done) {
            threadid += blockDim.x * blockDim.y;
            if(int64_t(threadid+blockDim.x * blockDim.y) >= int64_t(Kokkos::Impl::g_device_cuda_lock_arrays.n)) threadid = 0;
          }
        }
        base_thread_id = threadid;
      }
      __syncthreads();
      threadid = base_thread_id;
    }


    const int int_league_size = (int)m_league_size;
    for ( int league_rank = blockIdx.x ; league_rank < int_league_size ; league_rank += gridDim.x ) {

      this-> template exec_team< WorkTag >(
        typename Policy::member_type( kokkos_impl_cuda_shared_memory<void>()
                                    , m_shmem_begin
                                    , m_shmem_size
                                    , (void*) ( ((char*)m_scratch_ptr[1]) + ptrdiff_t(threadid/(blockDim.x*blockDim.y)) * m_scratch_size[1])
                                    , m_scratch_size[1]
                                    , league_rank
                                    , m_league_size ) );
    }
    if ( m_scratch_size[1]>0 ) {
      __syncthreads();
      if (threadIdx.x==0 && threadIdx.y==0 )
        Kokkos::Impl::g_device_cuda_lock_arrays.scratch[threadid]=0;
    }
  }

  inline
  void execute() const
    {
      const int64_t shmem_size_total = m_shmem_begin + m_shmem_size ;
      const dim3 grid( int(m_league_size) , 1 , 1 );
      const dim3 block( int(m_vector_size) , int(m_team_size) , 1 );

      CudaParallelLaunch< ParallelFor, LaunchBounds >( *this, grid, block, shmem_size_total, m_policy.space().impl_internal_space_instance() , true ); // copy to device and execute

    }

  ParallelFor( const FunctorType  & arg_functor
             , const Policy       & arg_policy
             )
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_league_size( arg_policy.league_size() )
    , m_team_size( arg_policy.team_size() )
    , m_vector_size( arg_policy.vector_length() )
    {
      cudaFuncAttributes attr = CudaParallelLaunch< ParallelFor, LaunchBounds >::
          get_cuda_func_attributes();
      m_team_size = m_team_size>=0?m_team_size:Kokkos::Impl::cuda_get_opt_block_size< FunctorType, LaunchBounds>(
        m_policy.space().impl_internal_space_instance(),
        attr, m_functor , m_vector_size,
        m_policy.team_scratch_size(0), m_policy.thread_scratch_size(0) )/m_vector_size;

      m_shmem_begin = ( sizeof(double) * ( m_team_size + 2 ) );
      m_shmem_size = ( m_policy.scratch_size(0,m_team_size) + FunctorTeamShmemSize< FunctorType >::value( m_functor , m_team_size ) );
      m_scratch_size[0] = m_policy.scratch_size(0,m_team_size);
      m_scratch_size[1] = m_policy.scratch_size(1,m_team_size);

      // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.
      m_scratch_ptr[0] = NULL;
      m_scratch_ptr[1] = m_team_size<=0?NULL:cuda_resize_scratch_space(static_cast<ptrdiff_t>(m_scratch_size[1])*static_cast<ptrdiff_t>(Cuda::concurrency()/(m_team_size*m_vector_size)));

      const int shmem_size_total = m_shmem_begin + m_shmem_size ;
      if ( m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock < shmem_size_total ) {
        printf("%i %i\n",m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock,shmem_size_total);
        Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
      }

      if ( int(m_team_size) >
           int(Kokkos::Impl::cuda_get_max_block_size< FunctorType, LaunchBounds >
                 ( m_policy.space().impl_internal_space_instance(),
        attr, arg_functor , arg_policy.vector_length(), arg_policy.team_scratch_size(0),arg_policy.thread_scratch_size(0) ) / arg_policy.vector_length())) {
        Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelFor< Cuda > requested too large team size."));
      }
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ReducerType, class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Traits ... >
                    , ReducerType
                    , Kokkos::Cuda
                    >
{
public:
  typedef Kokkos::RangePolicy< Traits ... >         Policy ;
private:


  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::launch_bounds LaunchBounds ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd, WorkTagFwd > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd, WorkTagFwd > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd, WorkTagFwd > ValueJoin ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::value_type      value_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Kokkos::Cuda::size_type                  size_type ;
  typedef typename Policy::index_type             index_type ;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y == blockDim.z == 1

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;
  const bool          m_result_ptr_device_accessible ;
  size_type *         m_scratch_space ;
  size_type *         m_scratch_flags ;
  size_type *         m_unified_space ;

  // Shall we use the shfl based reduction or not (only use it for static sized types of more than 128bit)
  enum { UseShflReduction = false };//((sizeof(value_type)>2*sizeof(double)) && ValueTraits::StaticValueSize) };
  // Some crutch to do function overloading
private:
  typedef double DummyShflReductionType;
  typedef int DummySHMEMReductionType;

public:
  // Make the exec_range calls call to Reduce::DeviceIterateTile
  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update ) const
    { m_functor( i , update ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update ) const
    { m_functor( TagType() , i , update ); }

  __device__ inline
  void operator() () const {
/*    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType, DummySHMEMReductionType>::select(1,1.0) );
  }

  __device__ inline
  void run(const DummySHMEMReductionType& ) const
  {*/
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) / sizeof(size_type) );

    {
      reference_type value =
        ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      const WorkRange range( m_policy , blockIdx.x , gridDim.x );

      for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
            iwork < iwork_end ; iwork += blockDim.y ) {
        this-> template exec_range< WorkTag >( iwork , value );
      }
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,ReducerTypeFwd,WorkTagFwd>(
           ReducerConditional::select(m_functor , m_reducer) , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_result_ptr_device_accessible? reinterpret_cast<size_type*>(m_result_ptr) : 
                                 ( m_unified_space ? m_unified_space : m_scratch_space );

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }

/*  __device__ inline
   void run(const DummyShflReductionType&) const
   {
     value_type value;
     ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &value);
     // Number of blocks is bounded so that the reduction can be limited to two passes.
     // Each thread block is given an approximately equal amount of work to perform.
     // Accumulate the values for this block.
     // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

     const WorkRange range( m_policy , blockIdx.x , gridDim.x );

     for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
           iwork < iwork_end ; iwork += blockDim.y ) {
       this-> template exec_range< WorkTag >( iwork , value );
     }

     pointer_type const result = (pointer_type) (m_unified_space ? m_unified_space : m_scratch_space) ;

     int max_active_thread = range.end()-range.begin() < blockDim.y ? range.end() - range.begin():blockDim.y;

     max_active_thread = (max_active_thread == 0)?blockDim.y:max_active_thread;

    value_type init;
    ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &init);
     if(Impl::cuda_inter_block_reduction<ReducerTypeFwd,ValueJoin,WorkTagFwd>
            (value,init,ValueJoin(ReducerConditional::select(m_functor , m_reducer)),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
       const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
       if(id==0) {
         Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , (void*) &value );
         *result = value;
       }
     }
   }*/

  // Determine block size constrained by shared memory:
  inline
  unsigned local_block_size( const FunctorType & f )
    {
      unsigned n = CudaTraits::WarpSize * 8 ;
      int shmem_size = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n );
      while ( (n && (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock < shmem_size)) ||
          (n > static_cast<unsigned>(Kokkos::Impl::cuda_get_max_block_size< ParallelReduce, LaunchBounds>( f , 1, shmem_size , 0 )))) {
        n >>= 1 ;
        shmem_size = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n );
      }
      return n ;
    }

  inline
  void execute()
    {
      const index_type nwork = m_policy.end() - m_policy.begin();
      if ( nwork ) {
        const int block_size = local_block_size( m_functor );

        m_scratch_space = cuda_internal_scratch_space( m_policy.space(), ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) * block_size /* block_size == max block_count */ );
        m_scratch_flags = cuda_internal_scratch_flags( m_policy.space(), sizeof(size_type) );
        m_unified_space = cuda_internal_scratch_unified( m_policy.space(), ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) );

        // REQUIRED ( 1 , N , 1 )
        const dim3 block( 1 , block_size , 1 );
        // Required grid.x <= block.y
        const dim3 grid( std::min( int(block.y) , int( ( nwork + block.y - 1 ) / block.y ) ) , 1 , 1 );

      const int shmem = UseShflReduction?0:cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( m_functor , block.y );

      CudaParallelLaunch< ParallelReduce, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute

      if(!m_result_ptr_device_accessible) {
        Cuda().fence();

        if ( m_result_ptr ) {
          if ( m_unified_space ) {
            const int count = ValueTraits::value_count( ReducerConditional::select(m_functor , m_reducer)  );
            for ( int i = 0 ; i < count ; ++i ) { m_result_ptr[i] = pointer_type(m_unified_space)[i] ; }
          }
          else {
            const int size = ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer)  );
            DeepCopy<HostSpace,CudaSpace>( m_result_ptr , m_scratch_space , size );
          }
        }
      }
    }
    else {
      if (m_result_ptr) {
        ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , m_result_ptr );
      }
    }
  }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType & arg_result
                , typename std::enable_if<
                   Kokkos::is_view< ViewType >::value
                ,void*>::type = NULL)
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( InvalidType() )
  , m_result_ptr( arg_result.data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ViewType::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  { }

  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ReducerType & reducer)
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( reducer )
  , m_result_ptr( reducer.view().data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ReducerType::result_view_type::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  { }
};


// MDRangePolicy impl
template< class FunctorType , class ReducerType, class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::MDRangePolicy< Traits ... >
                    , ReducerType
                    , Kokkos::Cuda
                    >
{
public:
  typedef Kokkos::MDRangePolicy< Traits ... > Policy ;
private:

  typedef typename Policy::array_index_type                 array_index_type;
  typedef typename Policy::index_type                       index_type;

  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::launch_bounds LaunchBounds;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd, WorkTagFwd > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd, WorkTagFwd > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd, WorkTagFwd > ValueJoin ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::value_type      value_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y == blockDim.z == 1

  const FunctorType   m_functor ;
  const Policy        m_policy ; // used for workrange and nwork
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;
  const bool          m_result_ptr_device_accessible ;
  size_type *         m_scratch_space ;
  size_type *         m_scratch_flags ;
  size_type *         m_unified_space ;

  typedef typename Kokkos::Impl::Reduce::DeviceIterateTile<Policy::rank, Policy, FunctorType, typename Policy::work_tag, reference_type> DeviceIteratePattern;

  // Shall we use the shfl based reduction or not (only use it for static sized types of more than 128bit
  enum { UseShflReduction = ((sizeof(value_type)>2*sizeof(double)) && (ValueTraits::StaticValueSize!=0)) };
  // Some crutch to do function overloading
private:
  typedef double DummyShflReductionType;
  typedef int DummySHMEMReductionType;

public:
  inline
  __device__
  void
  exec_range( reference_type update ) const
  {
    Kokkos::Impl::Reduce::DeviceIterateTile<Policy::rank,Policy,FunctorType,typename Policy::work_tag, reference_type>(m_policy, m_functor, update).exec_range();
  }

  inline
  __device__
  void operator() (void) const {
/*    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType, DummySHMEMReductionType>::select(1,1.0) );
  }

  __device__ inline
  void run(const DummySHMEMReductionType& ) const
  {*/
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) / sizeof(size_type) );

    {
      reference_type value =
        ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

      // Number of blocks is bounded so that the reduction can be limited to two passes.
      // Each thread block is given an approximately equal amount of work to perform.
      // Accumulate the values for this block.
      // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

      this-> exec_range( value );
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Problem: non power-of-two blockDim
    if ( cuda_single_inter_block_reduce_scan<false,ReducerTypeFwd,WorkTagFwd>(
           ReducerConditional::select(m_functor , m_reducer) , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location
      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_result_ptr_device_accessible? reinterpret_cast<size_type*>(m_result_ptr) :
                                 ( m_unified_space ? m_unified_space : m_scratch_space );

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }
  }

/*  __device__ inline
   void run(const DummyShflReductionType&) const
   {

     value_type value;
     ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &value);
     // Number of blocks is bounded so that the reduction can be limited to two passes.
     // Each thread block is given an approximately equal amount of work to perform.
     // Accumulate the values for this block.
     // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

     const Member work_part =
       ( ( m_policy.m_num_tiles + ( gridDim.x - 1 ) ) / gridDim.x ); //portion of tiles handled by each block

     this-> exec_range( value );

     pointer_type const result = (pointer_type) (m_unified_space ? m_unified_space : m_scratch_space) ;

     int max_active_thread = work_part < blockDim.y ? work_part:blockDim.y;
     max_active_thread = (max_active_thread == 0)?blockDim.y:max_active_thread;

     value_type init;
     ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &init);
     if(Impl::cuda_inter_block_reduction<ReducerTypeFwd,ValueJoin,WorkTagFwd>
         (value,init,ValueJoin(ReducerConditional::select(m_functor , m_reducer)),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
       const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
       if(id==0) {
         Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , (void*) &value );
         *result = value;
       }
     }
   }
*/
  // Determine block size constrained by shared memory:
  inline
  unsigned local_block_size( const FunctorType & f )
    {
      unsigned n = CudaTraits::WarpSize * 8 ;
      int shmem_size = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n );
      while ( (n && (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock < shmem_size)) ||
          (n > static_cast<unsigned>(Kokkos::Impl::cuda_get_max_block_size< ParallelReduce, LaunchBounds>( f , 1, shmem_size , 0 )))) {
        n >>= 1 ;
        shmem_size = cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n );
      }
      return n ;
    }

  inline
  void execute()
    {
      const int nwork = m_policy.m_num_tiles;
      if ( nwork ) {
        int block_size = m_policy.m_prod_tile_dims;
        // CONSTRAINT: Algorithm requires block_size >= product of tile dimensions
        // Nearest power of two
        int exponent_pow_two = std::ceil( std::log2(block_size) );
        block_size = std::pow(2, exponent_pow_two);
        int suggested_blocksize = local_block_size( m_functor );

        block_size = (block_size > suggested_blocksize) ? block_size : suggested_blocksize ; //Note: block_size must be less than or equal to 512


        m_scratch_space = cuda_internal_scratch_space( m_policy.space(), ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) * block_size /* block_size == max block_count */ );
        m_scratch_flags = cuda_internal_scratch_flags( m_policy.space(), sizeof(size_type) );
        m_unified_space = cuda_internal_scratch_unified( m_policy.space(), ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) );

        // REQUIRED ( 1 , N , 1 )
        const dim3 block( 1 , block_size , 1 );
        // Required grid.x <= block.y
        const dim3 grid( std::min( int(block.y) , int( nwork ) ) , 1 , 1 );

      const int shmem = UseShflReduction?0:cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( m_functor , block.y );

      CudaParallelLaunch< ParallelReduce, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute

      if(!m_result_ptr_device_accessible) {
        Cuda().fence();

        if ( m_result_ptr ) {
          if ( m_unified_space ) {
            const int count = ValueTraits::value_count( ReducerConditional::select(m_functor , m_reducer)  );
            for ( int i = 0 ; i < count ; ++i ) { m_result_ptr[i] = pointer_type(m_unified_space)[i] ; }
          }
          else {
            const int size = ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer)  );
            DeepCopy<HostSpace,CudaSpace>( m_result_ptr , m_scratch_space , size );
          }
        }
      }
    }
    else {
      if (m_result_ptr) {
        ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , m_result_ptr );
      }
    }
  }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType & arg_result
                , typename std::enable_if<
                   Kokkos::is_view< ViewType >::value
                ,void*>::type = NULL)
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( InvalidType() )
  , m_result_ptr( arg_result.data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ViewType::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  {}

  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ReducerType & reducer)
  : m_functor( arg_functor )
  , m_policy(  arg_policy )
  , m_reducer( reducer )
  , m_result_ptr( reducer.view().data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ReducerType::result_view_type::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  {}
};


//----------------------------------------------------------------------------

template< class FunctorType , class ReducerType, class ... Properties >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Properties ... >
                    , ReducerType
                    , Kokkos::Cuda
                    >
{
public:
  typedef TeamPolicyInternal< Kokkos::Cuda, Properties ... >  Policy ;
private:

  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::launch_bounds     LaunchBounds ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits< ReducerTypeFwd, WorkTagFwd > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd, WorkTagFwd > ValueInit ;
  typedef Kokkos::Impl::FunctorValueJoin<   ReducerTypeFwd, WorkTagFwd > ValueJoin ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef typename ValueTraits::value_type      value_type ;

public:

  typedef FunctorType      functor_type ;
  typedef Cuda::size_type  size_type ;

  enum { UseShflReduction = (true && (ValueTraits::StaticValueSize!=0)) };

private:
  typedef double DummyShflReductionType;
  typedef int DummySHMEMReductionType;

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y == blockDim.z == 1
  // shared memory utilization:
  //
  //  [ global reduce space ]
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;
  const bool          m_result_ptr_device_accessible ;
  size_type *         m_scratch_space ;
  size_type *         m_scratch_flags ;
  size_type *         m_unified_space ;
  size_type           m_team_begin ;
  size_type           m_shmem_begin ;
  size_type           m_shmem_size ;
  void*               m_scratch_ptr[2] ;
  int                 m_scratch_size[2] ;
  const size_type     m_league_size ;
  int                 m_team_size ;
  const size_type     m_vector_size ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_team( const Member & member , reference_type update ) const
    { m_functor( member , update ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_team( const Member & member , reference_type update ) const
    { m_functor( TagType() , member , update ); }

public:

  __device__ inline
  void operator() () const {
    int64_t threadid = 0;
    if ( m_scratch_size[1]>0 ) {
      __shared__ int64_t base_thread_id;
      if (threadIdx.x==0 && threadIdx.y==0 ) {
        threadid = (blockIdx.x*blockDim.z + threadIdx.z) %
          (Kokkos::Impl::g_device_cuda_lock_arrays.n / (blockDim.x * blockDim.y));
        threadid *= blockDim.x * blockDim.y;
        int done = 0;
        while (!done) {
          done = (0 == atomicCAS(&Kokkos::Impl::g_device_cuda_lock_arrays.scratch[threadid],0,1));
          if(!done) {
            threadid += blockDim.x * blockDim.y;
            if(int64_t(threadid + blockDim.x * blockDim.y) >= int64_t(Kokkos::Impl::g_device_cuda_lock_arrays.n)) threadid = 0;
          }
        }
        base_thread_id = threadid;
      }
      __syncthreads();
      threadid = base_thread_id;
    }

    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType, DummySHMEMReductionType>::select(1,1.0), threadid );
    if ( m_scratch_size[1]>0 ) {
      __syncthreads();
      if (threadIdx.x==0 && threadIdx.y==0 )
        Kokkos::Impl::g_device_cuda_lock_arrays.scratch[threadid]=0;
    }
  }

  __device__ inline
  void run(const DummySHMEMReductionType&, const int& threadid) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) / sizeof(size_type) );

    reference_type value =
      ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , kokkos_impl_cuda_shared_memory<size_type>() + threadIdx.y * word_count.value );

    // Iterate this block through the league
    const int int_league_size = (int)m_league_size;
    for ( int league_rank = blockIdx.x ; league_rank < int_league_size ; league_rank += gridDim.x ) {
      this-> template exec_team< WorkTag >
        ( Member( kokkos_impl_cuda_shared_memory<char>() + m_team_begin
                                        , m_shmem_begin
                                        , m_shmem_size
                                        , (void*) ( ((char*)m_scratch_ptr[1]) + ptrdiff_t(threadid/(blockDim.x*blockDim.y)) * m_scratch_size[1])
                                        , m_scratch_size[1]
                                        , league_rank
                                        , m_league_size )
        , value );
    }

    // Reduce with final value at blockDim.y - 1 location.
    if ( cuda_single_inter_block_reduce_scan<false,FunctorType,WorkTag>(
           ReducerConditional::select(m_functor , m_reducer) , blockIdx.x , gridDim.x ,
           kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags ) ) {

      // This is the final block with the final result at the final threads' location

      size_type * const shared = kokkos_impl_cuda_shared_memory<size_type>() + ( blockDim.y - 1 ) * word_count.value ;
      size_type * const global = m_result_ptr_device_accessible? reinterpret_cast<size_type*>(m_result_ptr) :
                                 ( m_unified_space ? m_unified_space : m_scratch_space );

      if ( threadIdx.y == 0 ) {
        Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , shared );
      }

      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); }

      for ( unsigned i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i]; }
    }

  }

  __device__ inline
  void run(const DummyShflReductionType&, const int& threadid) const
  {
    value_type value;
    ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &value);

    // Iterate this block through the league
    const int int_league_size = (int)m_league_size;
    for ( int league_rank = blockIdx.x ; league_rank < int_league_size ; league_rank += gridDim.x ) {
      this-> template exec_team< WorkTag >
        ( Member( kokkos_impl_cuda_shared_memory<char>() + m_team_begin
                                        , m_shmem_begin
                                        , m_shmem_size
                                        , (void*) ( ((char*)m_scratch_ptr[1]) + ptrdiff_t(threadid/(blockDim.x*blockDim.y)) * m_scratch_size[1])
                                        , m_scratch_size[1]
                                        , league_rank
                                        , m_league_size )
        , value );
    }

    pointer_type const result = m_result_ptr_device_accessible? m_result_ptr :
                                (pointer_type) ( m_unified_space ? m_unified_space : m_scratch_space );

    value_type init;
    ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , &init);
    if(
        Impl::cuda_inter_block_reduction<FunctorType,ValueJoin,WorkTag>
           (value,init,ValueJoin(ReducerConditional::select(m_functor , m_reducer)),m_scratch_space,result,m_scratch_flags,blockDim.y)
        //This breaks a test
        //   Kokkos::Impl::CudaReductionsFunctor<FunctorType,WorkTag,false,true>::scalar_inter_block_reduction(ReducerConditional::select(m_functor , m_reducer) , blockIdx.x , gridDim.x ,
        //              kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags)
    ) {
      const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
      if(id==0) {
        Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final( ReducerConditional::select(m_functor , m_reducer) , (void*) &value );
        *result = value;
      }
    }
  }

  inline
  void execute()
    {
      const int nwork = m_league_size * m_team_size ;
      if ( nwork ) {
        const int block_count = UseShflReduction? std::min( m_league_size , size_type(1024*32) )
          :std::min( int(m_league_size) , m_team_size );

        m_scratch_space = cuda_internal_scratch_space(m_policy.space(), ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) * block_count );
        m_scratch_flags = cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type) );
        m_unified_space = cuda_internal_scratch_unified( m_policy.space(),ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) ) );

        const dim3 block( m_vector_size , m_team_size , 1 );
        const dim3 grid( block_count , 1 , 1 );
        const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;

        CudaParallelLaunch< ParallelReduce, LaunchBounds >( *this, grid, block, shmem_size_total , m_policy.space().impl_internal_space_instance() , true ); // copy to device and execute

        if(!m_result_ptr_device_accessible) {
          Cuda().fence();

          if ( m_result_ptr ) {
            if ( m_unified_space ) {
              const int count = ValueTraits::value_count( ReducerConditional::select(m_functor , m_reducer) );
              for ( int i = 0 ; i < count ; ++i ) { m_result_ptr[i] = pointer_type(m_unified_space)[i] ; }
            }
            else {
              const int size = ValueTraits::value_size( ReducerConditional::select(m_functor , m_reducer) );
              DeepCopy<HostSpace,CudaSpace>( m_result_ptr, m_scratch_space, size );
            }
          }
        }
      }
      else {
        if (m_result_ptr) {
          ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , m_result_ptr );
        }
      }
    }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType & arg_result
                , typename std::enable_if<
                                   Kokkos::is_view< ViewType >::value
                                ,void*>::type = NULL)
  : m_functor( arg_functor )
  , m_policy ( arg_policy )
  , m_reducer( InvalidType() )
  , m_result_ptr( arg_result.data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ViewType::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_team_begin( 0 )
  , m_shmem_begin( 0 )
  , m_shmem_size( 0 )
  , m_scratch_ptr{NULL,NULL}
  , m_league_size( arg_policy.league_size() )
  , m_team_size( arg_policy.team_size() )
  , m_vector_size( arg_policy.vector_length() )
  {
    cudaFuncAttributes attr = CudaParallelLaunch< ParallelReduce, LaunchBounds >::
        get_cuda_func_attributes();
    m_team_size = m_team_size>=0?m_team_size:
        Kokkos::Impl::cuda_get_opt_block_size< FunctorType, LaunchBounds>(
      m_policy.space().impl_internal_space_instance(),
      attr, m_functor , m_vector_size,
      m_policy.team_scratch_size(0), m_policy.thread_scratch_size(0) )/m_vector_size;

    // Return Init value if the number of worksets is zero
    if( m_league_size*m_team_size == 0) {
      ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , arg_result.data() );
      return ;
    }

    m_team_begin = UseShflReduction?0:cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( arg_functor , m_team_size );
    m_shmem_begin = sizeof(double) * ( m_team_size + 2 );
    m_shmem_size = m_policy.scratch_size(0,m_team_size) + FunctorTeamShmemSize< FunctorType >::value( arg_functor , m_team_size );
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1,m_team_size);
    m_scratch_ptr[1] = m_team_size<=0?NULL:cuda_resize_scratch_space(static_cast<std::int64_t>(m_scratch_size[1])*(static_cast<std::int64_t>(Cuda::concurrency()/(m_team_size*m_vector_size))));

    // The global parallel_reduce does not support vector_length other than 1 at the moment
    if( (arg_policy.vector_length() > 1) && !UseShflReduction )
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a vector length of greater than 1 is not currently supported for CUDA for dynamic sized reduction types.");

    if( (m_team_size < 32) && !UseShflReduction )
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller than 32 is not currently supported with CUDA for dynamic sized reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;

    if (! Kokkos::Impl::is_integral_power_of_two( m_team_size )  && !UseShflReduction ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    if ( m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too much L0 scratch memory"));
    }

    if ( int(m_team_size) > arg_policy.team_size_max(m_functor,ParallelReduceTag()) ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too large team size."));
    }

  }

  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ReducerType & reducer)
  : m_functor( arg_functor )
  , m_policy( arg_policy )
  , m_reducer( reducer )
  , m_result_ptr( reducer.view().data() )
  , m_result_ptr_device_accessible(MemorySpaceAccess< Kokkos::CudaSpace , typename ReducerType::result_view_type::memory_space>::accessible )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_unified_space( 0 )
  , m_team_begin( 0 )
  , m_shmem_begin( 0 )
  , m_shmem_size( 0 )
  , m_scratch_ptr{NULL,NULL}
  , m_league_size( arg_policy.league_size() )
  , m_team_size( arg_policy.team_size() )
  , m_vector_size( arg_policy.vector_length() )
  {
    cudaFuncAttributes attr = CudaParallelLaunch< ParallelReduce, LaunchBounds >::
        get_cuda_func_attributes();
    m_team_size = m_team_size>=0?m_team_size:
        Kokkos::Impl::cuda_get_opt_block_size< FunctorType, LaunchBounds>(
      m_policy.space().impl_internal_space_instance(),
      attr, m_functor , m_vector_size,
      m_policy.team_scratch_size(0), m_policy.thread_scratch_size(0) )/m_vector_size;

    // Return Init value if the number of worksets is zero
    if( arg_policy.league_size() == 0) {
      ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , m_result_ptr );
      return ;
    }

    m_team_begin = UseShflReduction?0:cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( arg_functor , m_team_size );
    m_shmem_begin = sizeof(double) * ( m_team_size + 2 );
    m_shmem_size = m_policy.scratch_size(0,m_team_size) + FunctorTeamShmemSize< FunctorType >::value( arg_functor , m_team_size );
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1,m_team_size);
    m_scratch_ptr[1] = m_team_size<=0?NULL:cuda_resize_scratch_space(static_cast<ptrdiff_t>(m_scratch_size[1])*static_cast<ptrdiff_t>(Cuda::concurrency()/(m_team_size*m_vector_size)));

    // The global parallel_reduce does not support vector_length other than 1 at the moment
    if( (arg_policy.vector_length() > 1) && !UseShflReduction )
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a vector length of greater than 1 is not currently supported for CUDA for dynamic sized reduction types.");

    if( (m_team_size < 32) && !UseShflReduction )
      Impl::throw_runtime_exception( "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller than 32 is not currently supported with CUDA for dynamic sized reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size ;

    if ( (! Kokkos::Impl::is_integral_power_of_two( m_team_size )  && !UseShflReduction ) ||
         m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock < shmem_size_total ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }
    if ( int(m_team_size) > arg_policy.team_size_max(m_functor,ParallelReduceTag()) ) {
      Kokkos::Impl::throw_runtime_exception(std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too large team size."));
    }

  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Traits ... >
                  , Kokkos::Cuda
                  >
{
public:
  typedef Kokkos::RangePolicy< Traits ... >  Policy ;
private:
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::launch_bounds  LaunchBounds ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType, WorkTag > ValueOps ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

private:

  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type         m_final ;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final_result ) const
    { m_functor( i , update , final_result ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final_result ) const
    { m_functor( TagType() , i , update , final_result ); }

  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    size_type * const shared_value = kokkos_impl_cuda_shared_memory<size_type>() + word_count.value * threadIdx.y ;

    ValueInit::init( m_functor , shared_value );

    // Number of blocks is bounded so that the reduction can be limited to two passes.
    // Each thread block is given an approximately equal amount of work to perform.
    // Accumulate the values for this block.
    // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
          iwork < iwork_end ; iwork += blockDim.y ) {
      this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , blockIdx.x , gridDim.x , kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags );
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_data   = kokkos_impl_cuda_shared_memory<size_type>();
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.y ;
    size_type * const shared_accum  = shared_data + word_count.value * ( blockDim.y + 1 );

    // Starting value for this thread block is the previous block's total.
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.y ) {
      ValueInit::init( m_functor , shared_accum );
    }

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork_base = range.begin(); iwork_base < range.end() ; iwork_base += blockDim.y ) {
      #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
      unsigned MASK=KOKKOS_IMPL_CUDA_ACTIVEMASK;
      #endif
      const typename Policy::member_type iwork = iwork_base + threadIdx.y ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      ValueInit::init( m_functor , shared_prefix + word_count.value );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }
      #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
      KOKKOS_IMPL_CUDA_SYNCWARP_MASK(MASK);
      #else
      KOKKOS_IMPL_CUDA_SYNCWARP;
      #endif
      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix + word_count.value ) , false );
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , typename ValueTraits::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.y ;
        for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix ) , true );
      }
    }
  }

public:

  //----------------------------------------

  __device__ inline
  void operator()(void) const
  {
    if ( ! m_final ) {
      initial();
    }
    else {
      final();
    }
  }

  // Determine block size constrained by shared memory:
  inline
  unsigned local_block_size( const FunctorType & f )
    {
      // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
      // gridDim.x <= blockDim.y * blockDim.y
      //
      // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

      unsigned n = CudaTraits::WarpSize * 4 ;
      while ( n && unsigned(m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock) < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  inline
  void execute()
    {
      const int nwork    = m_policy.end() - m_policy.begin();
      if ( nwork ) {
        enum { GridMaxComputeCapability_2x = 0x0ffff };

        const int block_size = local_block_size( m_functor );

        const int grid_max =
          ( block_size * block_size ) < GridMaxComputeCapability_2x ?
          ( block_size * block_size ) : GridMaxComputeCapability_2x ;

        // At most 'max_grid' blocks:
        const int max_grid = std::min( int(grid_max) , int(( nwork + block_size - 1 ) / block_size ));

        // How much work per block:
        const int work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

        // How many block are really needed for this much work:
        const int grid_x = ( nwork + work_per_block - 1 ) / work_per_block ;

        m_scratch_space = cuda_internal_scratch_space( m_policy.space(), ValueTraits::value_size( m_functor ) * grid_x );
        m_scratch_flags = cuda_internal_scratch_flags( m_policy.space(), sizeof(size_type) * 1 );

        const dim3 grid( grid_x , 1 , 1 );
        const dim3 block( 1 , block_size , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )
        const int shmem = ValueTraits::value_size( m_functor ) * ( block_size + 2 );

        m_final = false ;
        CudaParallelLaunch< ParallelScan, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute

        m_final = true ;
        CudaParallelLaunch< ParallelScan, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute
      }
    }

  ParallelScan( const FunctorType  & arg_functor ,
                const Policy       & arg_policy )
  : m_functor( arg_functor )
  , m_policy( arg_policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_final( false )
  { }
};

//----------------------------------------------------------------------------
template< class FunctorType, class ReturnType, class ... Traits >
class ParallelScanWithTotal< FunctorType
                           , Kokkos::RangePolicy< Traits ... >
                           , ReturnType
                           , Kokkos::Cuda
                           >
{
public:
  typedef Kokkos::RangePolicy< Traits ... >  Policy ;

private:
  typedef typename Policy::member_type  Member ;
  typedef typename Policy::work_tag     WorkTag ;
  typedef typename Policy::WorkRange    WorkRange ;
  typedef typename Policy::launch_bounds  LaunchBounds ;

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType, WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType, WorkTag > ValueInit ;
  typedef Kokkos::Impl::FunctorValueOps<    FunctorType, WorkTag > ValueOps ;

public:

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;
  typedef FunctorType                           functor_type ;
  typedef Cuda::size_type                       size_type ;

private:

  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const FunctorType m_functor ;
  const Policy      m_policy ;
  size_type *       m_scratch_space ;
  size_type *       m_scratch_flags ;
  size_type         m_final ;
  ReturnType      & m_returnvalue;

  template< class TagType >
  __device__ inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final_result ) const
    { m_functor( i , update , final_result ); }

  template< class TagType >
  __device__ inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_range( const Member & i , reference_type update , const bool final_result ) const
    { m_functor( TagType() , i , update , final_result ); }

  //----------------------------------------

  __device__ inline
  void initial(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    size_type * const shared_value = kokkos_impl_cuda_shared_memory<size_type>() + word_count.value * threadIdx.y ;

    ValueInit::init( m_functor , shared_value );

    // Number of blocks is bounded so that the reduction can be limited to two passes.
    // Each thread block is given an approximately equal amount of work to perform.
    // Accumulate the values for this block.
    // The accumulation ordering does not match the final pass, but is arithmatically equivalent.

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( Member iwork = range.begin() + threadIdx.y , iwork_end = range.end() ;
          iwork < iwork_end ; iwork += blockDim.y ) {
      this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_value ) , false );
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups' totals.
    // Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i < gridDim.x
    cuda_single_inter_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , blockIdx.x , gridDim.x , kokkos_impl_cuda_shared_memory<size_type>() , m_scratch_space , m_scratch_flags );
  }

  //----------------------------------------

  __device__ inline
  void final(void) const
  {
    const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
      word_count( ValueTraits::value_size( m_functor ) / sizeof(size_type) );

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] , value[2] , ... }
    size_type * const shared_data   = kokkos_impl_cuda_shared_memory<size_type>();
    size_type * const shared_prefix = shared_data + word_count.value * threadIdx.y ;
    size_type * const shared_accum  = shared_data + word_count.value * ( blockDim.y + 1 );

    // Starting value for this thread block is the previous block's total.
    if ( blockIdx.x ) {
      size_type * const block_total = m_scratch_space + word_count.value * ( blockIdx.x - 1 );
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i] ; }
    }
    else if ( 0 == threadIdx.y ) {
      ValueInit::init( m_functor , shared_accum );
    }

    const WorkRange range( m_policy , blockIdx.x , gridDim.x );

    for ( typename Policy::member_type iwork_base = range.begin(); iwork_base < range.end() ; iwork_base += blockDim.y ) {
      #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
      unsigned MASK=KOKKOS_IMPL_CUDA_ACTIVEMASK;
      #endif

      const typename Policy::member_type iwork = iwork_base + threadIdx.y ;

      __syncthreads(); // Don't overwrite previous iteration values until they are used

      ValueInit::init( m_functor , shared_prefix + word_count.value );

      // Copy previous block's accumulation total into thread[0] prefix and inclusive scan value of this block
      for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i] ;
      }

      #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
      KOKKOS_IMPL_CUDA_SYNCWARP_MASK(MASK);
      #else
      KOKKOS_IMPL_CUDA_SYNCWARP;
      #endif
      if ( CudaTraits::WarpSize < word_count.value ) { __syncthreads(); } // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix + word_count.value ) , false );
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true,FunctorType,WorkTag>( m_functor , typename ValueTraits::pointer_type(shared_data+word_count.value) );

      {
        size_type * const block_total = shared_data + word_count.value * blockDim.y ;
        for ( unsigned i = threadIdx.y ; i < word_count.value ; ++i ) { shared_accum[i] = block_total[i]; }
      }

      // Call functor with exclusive scan value
      if ( iwork < range.end() ) {
        this-> template exec_range< WorkTag >( iwork , ValueOps::reference( shared_prefix ) , true );
      }
    }
  }

public:

  //----------------------------------------

  __device__ inline
  void operator()(void) const
  {
    if ( ! m_final ) {
      initial();
    }
    else {
      final();
    }
  }

  // Determine block size constrained by shared memory:
  inline
  unsigned local_block_size( const FunctorType & f )
    {
      // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512 (16 warps)
      // gridDim.x <= blockDim.y * blockDim.y
      //
      // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit testing

      unsigned n = CudaTraits::WarpSize * 4 ;
      while ( n && unsigned(m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock) < cuda_single_inter_block_reduce_scan_shmem<false,FunctorType,WorkTag>( f , n ) ) { n >>= 1 ; }
      return n ;
    }

  inline
  void execute()
    {
      const int nwork    = m_policy.end() - m_policy.begin();
      if ( nwork ) {
        enum { GridMaxComputeCapability_2x = 0x0ffff };

        const int block_size = local_block_size( m_functor );

        const int grid_max =
          ( block_size * block_size ) < GridMaxComputeCapability_2x ?
          ( block_size * block_size ) : GridMaxComputeCapability_2x ;

        // At most 'max_grid' blocks:
        const int max_grid = std::min( int(grid_max) , int(( nwork + block_size - 1 ) / block_size ));

        // How much work per block:
        const int work_per_block = ( nwork + max_grid - 1 ) / max_grid ;

        // How many block are really needed for this much work:
        const int grid_x = ( nwork + work_per_block - 1 ) / work_per_block ;

        m_scratch_space = cuda_internal_scratch_space( m_policy.space(), ValueTraits::value_size( m_functor ) * grid_x );
        m_scratch_flags = cuda_internal_scratch_flags( m_policy.space(), sizeof(size_type) * 1 );

        const dim3 grid( grid_x , 1 , 1 );
        const dim3 block( 1 , block_size , 1 ); // REQUIRED DIMENSIONS ( 1 , N , 1 )
        const int shmem = ValueTraits::value_size( m_functor ) * ( block_size + 2 );

        m_final = false ;
        CudaParallelLaunch< ParallelScanWithTotal, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute

        m_final = true ;
        CudaParallelLaunch< ParallelScanWithTotal, LaunchBounds >( *this, grid, block, shmem , m_policy.space().impl_internal_space_instance() , false ); // copy to device and execute

        const int size = ValueTraits::value_size( m_functor );
        DeepCopy<HostSpace,CudaSpace>( &m_returnvalue, m_scratch_space + (grid_x - 1)*size/sizeof(int), size );
      }
    }

  ParallelScanWithTotal( const FunctorType  & arg_functor ,
                         const Policy       & arg_policy ,   
                         ReturnType         & arg_returnvalue )
  : m_functor( arg_functor )
  , m_policy( arg_policy )
  , m_scratch_space( 0 )
  , m_scratch_flags( 0 )
  , m_final( false )
  , m_returnvalue( arg_returnvalue )
  { }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {
  template< class FunctorType, class ExecPolicy, class ValueType , class Tag = typename ExecPolicy::work_tag>
  struct CudaFunctorAdapter {
    const FunctorType f;
    typedef ValueType value_type;
    CudaFunctorAdapter(const FunctorType& f_):f(f_) {}

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, j, val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, j, k, val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, j, k, l, val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, const typename ExecPolicy::member_type& m, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, j, k, l, m, val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::work_tag, const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, const typename ExecPolicy::member_type& m, const typename ExecPolicy::member_type& n, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals third argument type of FunctorType::operator()
      f(typename ExecPolicy::work_tag(), i, j, k, l, m, n, val);
    }

  };

  template< class FunctorType, class ExecPolicy, class ValueType >
  struct CudaFunctorAdapter<FunctorType,ExecPolicy,ValueType,void> {
    const FunctorType f;
    typedef ValueType value_type;
    CudaFunctorAdapter(const FunctorType& f_):f(f_) {}

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,val);
    }

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,val);
    }

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,val);
    }

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,val);
    }

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, const typename ExecPolicy::member_type& m, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,m,val);
    }

    __device__ inline
    void operator() (const typename ExecPolicy::member_type& i, const typename ExecPolicy::member_type& j, const typename ExecPolicy::member_type& k, const typename ExecPolicy::member_type& l, const typename ExecPolicy::member_type& m, const typename ExecPolicy::member_type& n, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,m,n,val);
    }


    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, typename ExecPolicy::member_type& j, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, typename ExecPolicy::member_type& j, typename ExecPolicy::member_type& k, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, typename ExecPolicy::member_type& j, typename ExecPolicy::member_type& k, typename ExecPolicy::member_type& l, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, typename ExecPolicy::member_type& j, typename ExecPolicy::member_type& k, typename ExecPolicy::member_type& l, typename ExecPolicy::member_type& m, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,m,val);
    }

    __device__ inline
    void operator() (typename ExecPolicy::member_type& i, typename ExecPolicy::member_type& j, typename ExecPolicy::member_type& k, typename ExecPolicy::member_type& l, typename ExecPolicy::member_type& m, typename ExecPolicy::member_type& n, ValueType& val) const {
      //Insert Static Assert with decltype on ValueType equals second argument type of FunctorType::operator()
      f(i,j,k,l,m,n,val);
    }

  };

  template<class FunctorType, class ResultType, class Tag, bool Enable = IsNonTrivialReduceFunctor<FunctorType>::value >
  struct FunctorReferenceType {
    typedef ResultType& reference_type;
  };

  template<class FunctorType, class ResultType, class Tag>
  struct FunctorReferenceType<FunctorType, ResultType, Tag, true> {
    typedef typename Kokkos::Impl::FunctorValueTraits< FunctorType ,Tag >::reference_type reference_type;
  };

  template< class FunctorTypeIn, class ExecPolicy, class ValueType>
  struct ParallelReduceFunctorType<FunctorTypeIn,ExecPolicy,ValueType,Cuda> {

    enum {FunctorHasValueType = IsNonTrivialReduceFunctor<FunctorTypeIn>::value };
    typedef typename Kokkos::Impl::if_c<FunctorHasValueType, FunctorTypeIn, Impl::CudaFunctorAdapter<FunctorTypeIn,ExecPolicy,ValueType> >::type functor_type;
    static functor_type functor(const FunctorTypeIn& functor_in) {
      return Impl::if_c<FunctorHasValueType,FunctorTypeIn,functor_type>::select(functor_in,functor_type(functor_in));
    }
  };

}

} // namespace Kokkos

#endif /* defined( __CUDACC__ ) */
#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */

