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

#ifndef KOKKOS_ROCM_EXP_ITERATE_TILE_REFACTOR_HPP
#define KOKKOS_ROCM_EXP_ITERATE_TILE_REFACTOR_HPP

#include <Kokkos_Macros.hpp>
#if defined( __HCC__ ) && defined( KOKKOS_ENABLE_ROCM )

#include <iostream>
#include <algorithm>
#include <cstdio>

#include <utility>

// #include<ROCm/Kokkos_ROCmExec.hpp>
// Including the file above leads to following type of errors:
// /home/ndellin/kokkos/core/src/ROCm/Kokkos_ROCmExec.hpp(84): error: incomplete type is not allowed
// use existing Kokkos functionality, e.g. max blocks, once resolved

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>
#endif


#define threadIdx_x (hc_get_workitem_id(0))
#define threadIdx_y (hc_get_workitem_id(1))
#define threadIdx_z (hc_get_workitem_id(2))

#define blockIdx_x  (hc_get_group_id(0))
#define blockIdx_y  (hc_get_group_id(1))
#define blockIdx_z  (hc_get_group_id(2))

#define blockDim_x  (hc_get_group_size(0))
#define blockDim_y  (hc_get_group_size(1))
#define blockDim_z  (hc_get_group_size(2))

#define gridDim_x   (hc_get_num_groups(0))
#define gridDim_y   (hc_get_num_groups(1))
#define gridDim_z   (hc_get_num_groups(2))


namespace Kokkos { namespace Impl {

namespace Refactor {

// ------------------------------------------------------------------ //
// ParallelFor iteration pattern
template< int N , typename RP , typename Functor , typename Tag >
struct DeviceIterateTile;

//Rank 2
// Specializations for void tag type
template< typename RP , typename Functor >
struct DeviceIterateTile<2,RP,Functor,void >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    // LL
    if (RP::inner_direction == RP::Left) {
      for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
        const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
        if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

          for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
            const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
            if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {
              m_func(offset_0 , offset_1);
            }
          }
        }
      }
    }
    // LR
    else {
      for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
        const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {
              m_func(offset_0 , offset_1);
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

// Specializations for tag type
template< typename RP , typename Functor , typename Tag >
struct DeviceIterateTile<2,RP,Functor,Tag>
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if (RP::inner_direction == RP::Left) {
      // Loop over size maxnumblocks until full range covered
      for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
        const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
        if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

          for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
            const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
            if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {
              m_func(Tag(), offset_0 , offset_1);
            }
          }
        }
      }
    }
    else {
      for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
        const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {
              m_func(Tag(), offset_0 , offset_1);
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};


//Rank 3
// Specializations for void tag type
template< typename RP , typename Functor >
struct DeviceIterateTile<3,RP,Functor,void >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    // LL
    if (RP::inner_direction == RP::Left) {
      for ( index_type tile_id2 = (index_type)blockIdx_z; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_z ) {
        const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[2];
        if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_z < m_rp.m_tile[2] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

              for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
                const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
                if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {
                  m_func(offset_0 , offset_1 , offset_2);
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
        const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

              for ( index_type tile_id2 = (index_type)blockIdx_z; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_z ) {
                const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_z < m_rp.m_tile[2] ) {
                  m_func(offset_0 , offset_1 , offset_2);
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

// Specializations for void tag type
template< typename RP , typename Functor , typename Tag >
struct DeviceIterateTile<3,RP,Functor,Tag>
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if (RP::inner_direction == RP::Left) {
      for ( index_type tile_id2 = (index_type)blockIdx_z; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_z ) {
        const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[2];
        if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_z < m_rp.m_tile[2] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

              for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
                const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
                if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {
                  m_func(Tag(), offset_0 , offset_1 , offset_2);
                }
              }
            }
          }
        }
      }
    }
    else {
      for ( index_type tile_id0 = (index_type)blockIdx_x; tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim_x ) {
        const index_type offset_0 = tile_id0*m_rp.m_tile[0] + (index_type)threadIdx_x + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && (index_type)threadIdx_x < m_rp.m_tile[0] ) {

          for ( index_type tile_id1 = (index_type)blockIdx_y; tile_id1 < m_rp.m_tile_end[1]; tile_id1 += gridDim_y ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && (index_type)threadIdx_y < m_rp.m_tile[1] ) {

              for ( index_type tile_id2 = (index_type)blockIdx_z; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_z ) {
                const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_z < m_rp.m_tile[2] ) {
                  m_func(Tag(), offset_0 , offset_1 , offset_2);
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};


//Rank 4
// Specializations for void tag type
template< typename RP , typename Functor >
struct DeviceIterateTile<4,RP,Functor,void >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    // LL
    if (RP::inner_direction == RP::Left) {
      const index_type temp0  =  m_rp.m_tile_end[0];
      const index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      for ( index_type tile_id3 = (index_type)blockIdx_z; tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim_z ) {
        const index_type offset_3 = tile_id3*m_rp.m_tile[3] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[3];
        if ( offset_3 < m_rp.m_upper[3] && (index_type)threadIdx_z < m_rp.m_tile[3] ) {

          for ( index_type tile_id2 = (index_type)blockIdx_y; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_y ) {
            const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[2];
            if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_y < m_rp.m_tile[2] ) {

              for ( index_type j = tile_id1 ; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                  for ( index_type i = tile_id0 ; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                    const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                    if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                      m_func(offset_0 , offset_1 , offset_2 , offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      const index_type temp0  =  m_rp.m_tile_end[0];
      const index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type tile_id2 = (index_type)blockIdx_y; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_y ) {
                const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_y < m_rp.m_tile[2] ) {

                  for ( index_type tile_id3 = (index_type)blockIdx_z; tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim_z ) {
                    const index_type offset_3 = tile_id3*m_rp.m_tile[3] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && (index_type)threadIdx_z < m_rp.m_tile[3] ) {
                      m_func(offset_0 , offset_1 , offset_2 , offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

// Specializations for void tag type
template< typename RP , typename Functor , typename Tag >
struct DeviceIterateTile<4,RP,Functor,Tag>
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if (RP::inner_direction == RP::Left) {
      const index_type temp0  =  m_rp.m_tile_end[0];
      const index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      for ( index_type tile_id3 = (index_type)blockIdx_z; tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim_z ) {
        const index_type offset_3 = tile_id3*m_rp.m_tile[3] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[3];
        if ( offset_3 < m_rp.m_upper[3] && (index_type)threadIdx_z < m_rp.m_tile[3] ) {

          for ( index_type tile_id2 = (index_type)blockIdx_y; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_y ) {
            const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[2];
            if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_y < m_rp.m_tile[2] ) {

              for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                  for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                    const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                    if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                      m_func(Tag(), offset_0 , offset_1 , offset_2 , offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else {
      const index_type temp0  =  m_rp.m_tile_end[0];
      const index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = tile_id1*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type tile_id2 = (index_type)blockIdx_y; tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim_y ) {
                const index_type offset_2 = tile_id2*m_rp.m_tile[2] + (index_type)threadIdx_y + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && (index_type)threadIdx_y < m_rp.m_tile[2] ) {

                  for ( index_type tile_id3 = (index_type)blockIdx_z; tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim_z ) {
                    const index_type offset_3 = tile_id3*m_rp.m_tile[3] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && (index_type)threadIdx_z < m_rp.m_tile[3] ) {
                      m_func(Tag() , offset_0 , offset_1 , offset_2 , offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};


//Rank 5
// Specializations for void tag type
template< typename RP , typename Functor >
struct DeviceIterateTile<5,RP,Functor,void >
{
  using index_type = typename RP::index_type;

  
  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    // LL
    if (RP::inner_direction == RP::Left) {

      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl2 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl3 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl2 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y % numbl2;
      const index_type tile_id3 = (index_type)blockIdx_y / numbl2;
      const index_type thr_id2 = (index_type)threadIdx_y % m_rp.m_tile[2];
      const index_type thr_id3 = (index_type)threadIdx_y / m_rp.m_tile[2];

      for ( index_type tile_id4 = (index_type)blockIdx_z; tile_id4 < m_rp.m_tile_end[4]; tile_id4 += gridDim_z ) {
        const index_type offset_4 = tile_id4*m_rp.m_tile[4] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[4];
        if ( offset_4 < m_rp.m_upper[4] && (index_type)threadIdx_z < m_rp.m_tile[4] ) {

          for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
            const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
            if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type j = tile_id1 ; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                    const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                    if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                      for ( index_type i = tile_id0 ; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                          m_func(offset_0 , offset_1 , offset_2 , offset_3, offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl3 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl2 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl3 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y / numbl3;
      const index_type tile_id3 = (index_type)blockIdx_y % numbl3;
      const index_type thr_id2 = (index_type)threadIdx_y / m_rp.m_tile[3];
      const index_type thr_id3 = (index_type)threadIdx_y % m_rp.m_tile[3];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                    const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                      for ( index_type tile_id4 = (index_type)blockIdx_z; tile_id4 < m_rp.m_tile_end[4]; tile_id4 += gridDim_z ) {
                        const index_type offset_4 = tile_id4*m_rp.m_tile[4] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[4];
                        if ( offset_4 < m_rp.m_upper[4] && (index_type)threadIdx_z < m_rp.m_tile[4] ) {
                          m_func(offset_0 , offset_1 , offset_2 , offset_3 , offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

// Specializations for tag type
template< typename RP , typename Functor , typename Tag >
struct DeviceIterateTile<5,RP,Functor,Tag>
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl2 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl3 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl2 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y % numbl2;
      const index_type tile_id3 = (index_type)blockIdx_y / numbl2;
      const index_type thr_id2 = (index_type)threadIdx_y % m_rp.m_tile[2];
      const index_type thr_id3 = (index_type)threadIdx_y / m_rp.m_tile[2];

      for ( index_type tile_id4 = (index_type)blockIdx_z; tile_id4 < m_rp.m_tile_end[4]; tile_id4 += gridDim_z ) {
        const index_type offset_4 = tile_id4*m_rp.m_tile[4] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[4];
        if ( offset_4 < m_rp.m_upper[4] && (index_type)threadIdx_z < m_rp.m_tile[4] ) {

          for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
            const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
            if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type j = tile_id1 ; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                    const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                    if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                      for ( index_type i = tile_id0 ; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                          m_func(Tag() , offset_0 , offset_1 , offset_2 , offset_3, offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl3 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl2 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl3 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y / numbl3;
      const index_type tile_id3 = (index_type)blockIdx_y % numbl3;
      const index_type thr_id2 = (index_type)threadIdx_y / m_rp.m_tile[3];
      const index_type thr_id3 = (index_type)threadIdx_y % m_rp.m_tile[3];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                    const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                      for ( index_type tile_id4 = (index_type)blockIdx_z; tile_id4 < m_rp.m_tile_end[4]; tile_id4 += gridDim_z ) {
                        const index_type offset_4 = tile_id4*m_rp.m_tile[4] + (index_type)threadIdx_z + (index_type)m_rp.m_lower[4];
                        if ( offset_4 < m_rp.m_upper[4] && (index_type)threadIdx_z < m_rp.m_tile[4] ) {
                          m_func(Tag() , offset_0 , offset_1 , offset_2 , offset_3 , offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};


//Rank 6
// Specializations for void tag type
template< typename RP , typename Functor >
struct DeviceIterateTile<6,RP,Functor,void >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl2 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl3 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl2 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y % numbl2;
      const index_type tile_id3 = (index_type)blockIdx_y / numbl2;
      const index_type thr_id2 = (index_type)threadIdx_y % m_rp.m_tile[2];
      const index_type thr_id3 = (index_type)threadIdx_y / m_rp.m_tile[2];

      temp0  =  m_rp.m_tile_end[4];
      temp1  =  m_rp.m_tile_end[5];
      const index_type numbl4 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl5 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl4 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id4 = (index_type)blockIdx_z % numbl4;
      const index_type tile_id5 = (index_type)blockIdx_z / numbl4;
      const index_type thr_id4 = (index_type)threadIdx_z % m_rp.m_tile[4];
      const index_type thr_id5 = (index_type)threadIdx_z / m_rp.m_tile[4];

      for ( index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5 ) {
        const index_type offset_5 = n*m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
        if ( offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5] ) {

          for ( index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4 ) {
            const index_type offset_4 = m*m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
            if ( offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4] ) {

              for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                  for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                    const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                    if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                      for ( index_type j = tile_id1 ; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                        const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                        if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                          for ( index_type i = tile_id0 ; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                            const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                            if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                              m_func(offset_0 , offset_1 , offset_2 , offset_3, offset_4, offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl3 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl2 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl3 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y / numbl3;
      const index_type tile_id3 = (index_type)blockIdx_y % numbl3;
      const index_type thr_id2 = (index_type)threadIdx_y / m_rp.m_tile[3];
      const index_type thr_id3 = (index_type)threadIdx_y % m_rp.m_tile[3];

      temp0  =  m_rp.m_tile_end[4];
      temp1  =  m_rp.m_tile_end[5];
      const index_type numbl5 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl4 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl5 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id4 = (index_type)blockIdx_z / numbl5;
      const index_type tile_id5 = (index_type)blockIdx_z % numbl5;
      const index_type thr_id4 = (index_type)threadIdx_z / m_rp.m_tile[5];
      const index_type thr_id5 = (index_type)threadIdx_z % m_rp.m_tile[5];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                    const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                      for ( index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4 ) {
                        const index_type offset_4 = m*m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
                        if ( offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4] ) {

                          for ( index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5 ) {
                            const index_type offset_5 = n*m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
                            if ( offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5] ) {
                              m_func(offset_0 , offset_1 , offset_2 , offset_3 , offset_4 , offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

// Specializations for tag type
template< typename RP , typename Functor , typename Tag >
struct DeviceIterateTile<6,RP,Functor,Tag>
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ )
  : m_rp(rp_)
  , m_func(f_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl0 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl1 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl0 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x % numbl0;
      const index_type tile_id1 = (index_type)blockIdx_x / numbl0;
      const index_type thr_id0 = (index_type)threadIdx_x % m_rp.m_tile[0];
      const index_type thr_id1 = (index_type)threadIdx_x / m_rp.m_tile[0];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl2 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl3 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl2 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y % numbl2;
      const index_type tile_id3 = (index_type)blockIdx_y / numbl2;
      const index_type thr_id2 = (index_type)threadIdx_y % m_rp.m_tile[2];
      const index_type thr_id3 = (index_type)threadIdx_y / m_rp.m_tile[2];

      temp0  =  m_rp.m_tile_end[4];
      temp1  =  m_rp.m_tile_end[5];
      const index_type numbl4 = ( temp0 <= max_blocks ? temp0 : max_blocks ) ;
      const index_type numbl5 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl4 ) :
          (  temp1 <= max_blocks ? temp1 : max_blocks ) );

      const index_type tile_id4 = (index_type)blockIdx_z % numbl4;
      const index_type tile_id5 = (index_type)blockIdx_z / numbl4;
      const index_type thr_id4 = (index_type)threadIdx_z % m_rp.m_tile[4];
      const index_type thr_id5 = (index_type)threadIdx_z / m_rp.m_tile[4];

      for ( index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5 ) {
        const index_type offset_5 = n*m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
        if ( offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5] ) {

          for ( index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4 ) {
            const index_type offset_4 = m*m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
            if ( offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4] ) {

              for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                  for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                    const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                    if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                      for ( index_type j = tile_id1 ; j < m_rp.m_tile_end[1]; j += numbl1 ) {
                        const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                        if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

                          for ( index_type i = tile_id0 ; i < m_rp.m_tile_end[0]; i += numbl0 ) {
                            const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
                            if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {
                              m_func(Tag() , offset_0 , offset_1 , offset_2 , offset_3, offset_4, offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0  =  m_rp.m_tile_end[0];
      index_type temp1  =  m_rp.m_tile_end[1];
      const index_type numbl1 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl0 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl1 ) :
          ( temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id0 = (index_type)blockIdx_x / numbl1;
      const index_type tile_id1 = (index_type)blockIdx_x % numbl1;
      const index_type thr_id0 = (index_type)threadIdx_x / m_rp.m_tile[1];
      const index_type thr_id1 = (index_type)threadIdx_x % m_rp.m_tile[1];

      temp0  =  m_rp.m_tile_end[2];
      temp1  =  m_rp.m_tile_end[3];
      const index_type numbl3 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl2 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl3 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id2 = (index_type)blockIdx_y / numbl3;
      const index_type tile_id3 = (index_type)blockIdx_y % numbl3;
      const index_type thr_id2 = (index_type)threadIdx_y / m_rp.m_tile[3];
      const index_type thr_id3 = (index_type)threadIdx_y % m_rp.m_tile[3];

      temp0  =  m_rp.m_tile_end[4];
      temp1  =  m_rp.m_tile_end[5];
      const index_type numbl5 = ( temp1 <= max_blocks ? temp1 : max_blocks ) ;
      const index_type numbl4 = ( temp0*temp1 > max_blocks ? index_type( max_blocks / numbl5 ) :
          (  temp0 <= max_blocks ? temp0 : max_blocks ) );

      const index_type tile_id4 = (index_type)blockIdx_z / numbl5;
      const index_type tile_id5 = (index_type)blockIdx_z % numbl5;
      const index_type thr_id4 = (index_type)threadIdx_z / m_rp.m_tile[5];
      const index_type thr_id5 = (index_type)threadIdx_z % m_rp.m_tile[5];

      for ( index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0 ) {
        const index_type offset_0 = i*m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if ( offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0] ) {

          for ( index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1 ) {
            const index_type offset_1 = j*m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if ( offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1] ) {

              for ( index_type k = tile_id2; k < m_rp.m_tile_end[2]; k += numbl2 ) {
                const index_type offset_2 = k*m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if ( offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2] ) {

                  for ( index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3 ) {
                    const index_type offset_3 = l*m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                    if ( offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3] ) {

                      for ( index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4 ) {
                        const index_type offset_4 = m*m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
                        if ( offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4] ) {

                          for ( index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5 ) {
                            const index_type offset_5 = n*m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
                            if ( offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5] ) {
                              m_func(Tag() , offset_0 , offset_1 , offset_2 , offset_3 , offset_4 , offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
};

} // Refactor

// ----------------------------------------------------------------------------------

namespace Reduce {

template < typename T >
using is_void = std::is_same< T, void >;

template < typename T >
struct is_array_type : std::false_type
{
  using value_type = T;
};

template < typename T >
struct is_array_type< T* > : std::true_type
{
  using value_type = T;
};

template < typename T >
struct is_array_type< T[] > : std::true_type
{
  using value_type = T;
};

// ------------------------------------------------------------------ //
template< int N , typename RP , typename Functor , typename Tag , typename ValueType , typename Enable = void >
struct DeviceIterateTile;

// ParallelReduce iteration pattern
// Scalar reductions

// num_blocks = min( num_tiles, max_num_blocks ); //i.e. determined by number of tiles and reduction algorithm constraints
// extract n-dim tile offsets (i.e. tile's global starting mulit-index) from the tileid = blockid using tile dimensions
// local indices within a tile extracted from (index_type)threadIdx_x using tile dims, constrained by blocksize
// combine tile and local id info for multi-dim global ids

// Pattern:
// Each block+thread is responsible for a tile+local_id combo (additional when striding by num_blocks)
// 1. create offset arrays
// 2. loop over number of tiles, striding by griddim (equal to num tiles, or max num blocks)
// 3. temps set for tile_idx and thrd_idx, which will be modified
// 4. if LL vs LR:
//      determine tile starting point offsets (multidim)
//      determine local index offsets (multidim)
//      concatentate tile offset + local offset for global multi-dim index
//    if offset withinin range bounds AND local offset within tile bounds, call functor

// ValueType = T
//Rank 2
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<2,RP,Functor,void,ValueType, typename std::enable_if< !is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            // Deduce this blocks tile_id
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          {
            m_func( m_offset[0], m_offset[1], m_v );
          }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_v ); }
        }
      }
    }

  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<2,RP,Functor,Tag, ValueType, typename std::enable_if< !is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


//Rank 3
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<3,RP,Functor,void,ValueType , typename std::enable_if< !is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// Specializations for void tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<3,RP,Functor,Tag, ValueType, typename std::enable_if< !is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


//Rank 4
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<4,RP,Functor,void,ValueType , typename std::enable_if< !is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// Specializations for void tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<4,RP,Functor,Tag,ValueType, typename std::enable_if< !is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


//Rank 5
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<5,RP,Functor,void,ValueType , typename std::enable_if< !is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<5,RP,Functor,Tag,ValueType, typename std::enable_if< !is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


//Rank 6
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<6,RP,Functor,void,ValueType , typename std::enable_if< !is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<6,RP,Functor,Tag,ValueType, typename std::enable_if< !is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , ValueType & v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  ValueType & m_v;
};


// ValueType = T[], T*
//Rank 2
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<2,RP,Functor,void,ValueType, typename std::enable_if< is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  [[hc]]
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          {
            m_func( m_offset[0], m_offset[1], m_v );
          }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<2,RP,Functor,Tag, ValueType, typename std::enable_if< is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_v ); }
        }
      } //end for loop over num_tiles - product of tiles in each direction
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


//Rank 3
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<3,RP,Functor,void,ValueType , typename std::enable_if< is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]); // Move this to first computation, add to m_offset right away
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


// Specializations for void tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<3,RP,Functor,Tag, ValueType, typename std::enable_if< is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


//Rank 4
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<4,RP,Functor,void,ValueType , typename std::enable_if< is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


// Specializations for void tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<4,RP,Functor,Tag,ValueType, typename std::enable_if< is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


//Rank 5
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<5,RP,Functor,void,ValueType , typename std::enable_if< is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<5,RP,Functor,Tag,ValueType, typename std::enable_if< is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


//Rank 6
// Specializations for void tag type
template< typename RP , typename Functor , typename ValueType >
struct DeviceIterateTile<6,RP,Functor,void,ValueType , typename std::enable_if< is_array_type<ValueType>::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};


// Specializations for tag type
template< typename RP , typename Functor , typename Tag, typename ValueType >
struct DeviceIterateTile<6,RP,Functor,Tag,ValueType, typename std::enable_if< is_array_type<ValueType>::value && !is_void< Tag >::value >::type >
{
  using index_type = typename RP::index_type;
  using value_type = typename is_array_type< ValueType >::value_type;

  KOKKOS_INLINE_FUNCTION
  DeviceIterateTile( const RP & rp_ , const Functor & f_ , value_type* v_)
  : m_rp(rp_)
  , m_func(f_)
  , m_v(v_)
  {}

  static constexpr index_type max_blocks = 65535;
  //static constexpr index_type max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount);

  KOKKOS_INLINE_FUNCTION
  void exec_range() const
  {
    //enum { max_blocks = static_cast<index_type>(Kokkos::Impl::ROCmTraits::UpperBoundGridCount) };
    //const index_type max_blocks = static_cast<index_type>( Kokkos::Impl::rocm_internal_maximum_grid_count() );
    if ( (index_type)blockIdx_x < m_rp.m_num_tiles && (index_type)threadIdx_y < m_rp.m_prod_tile_dims ) {
      index_type m_offset[RP::rank]; // tile starting global id offset
      index_type m_local_offset[RP::rank]; // tile starting global id offset

      for ( index_type tileidx = (index_type)blockIdx_x; tileidx < m_rp.m_num_tiles; tileidx += gridDim_x ) {
        index_type tile_idx = tileidx; // temp because tile_idx will be modified while determining tile starting point offsets
        index_type thrd_idx = (index_type)threadIdx_y;
        bool in_bounds = true;

        // LL
        if (RP::inner_direction == RP::Left) {
          for (int i=0; i<RP::rank; ++i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
        // LR
        else {
          for (int i=RP::rank-1; i>=0; --i) {
            m_offset[i] = (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i] ;
            tile_idx /= m_rp.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_rp.m_tile[i]);
            thrd_idx /= m_rp.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if ( !(m_offset[i] < m_rp.m_upper[i] && m_local_offset[i] < m_rp.m_tile[i]) ) {
              in_bounds &= false;
            }
          }
          if ( in_bounds )
          { m_func( Tag(), m_offset[0], m_offset[1], m_offset[2], m_offset[3], m_offset[4], m_offset[5], m_v ); }
        }
      }
    }
  } //end exec_range

private:
  const RP & m_rp;
  const Functor & m_func;
  value_type* m_v;
};

} // Reduce

// ----------------------------------------------------------------------------------

} } //end namespace Kokkos::Impl
#undef threadIdx_x
#undef threadIdx_y
#undef threadIdx_z

#undef blockIdx_x
#undef blockIdx_y
#undef blockIdx_z

#undef blockDim_x
#undef blockDim_y
#undef blockDim_z

#undef gridDim_x
#undef gridDim_y
#undef gridDim_z


#endif
#endif
