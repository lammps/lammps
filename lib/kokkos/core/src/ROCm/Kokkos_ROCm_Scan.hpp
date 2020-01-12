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

#include <ROCm/Kokkos_ROCm_Invoke.hpp>
#include <ROCm/Kokkos_ROCm_Join.hpp>

namespace Kokkos {
namespace Impl {

template< class Tag, class F, class TransformIndex>
void scan_enqueue(
  const int len,
  const F & f,
  TransformIndex transform_index)
{
    typedef Kokkos::Impl::FunctorValueTraits< F, Tag>  ValueTraits;
    typedef Kokkos::Impl::FunctorValueInit<   F, Tag>  ValueInit;
    typedef Kokkos::Impl::FunctorValueJoin<   F, Tag>  ValueJoin;
    typedef Kokkos::Impl::FunctorValueOps<    F, Tag>  ValueOps;

    typedef typename ValueTraits::value_type    value_type;
    typedef typename ValueTraits::pointer_type    pointer_type;
    typedef typename ValueTraits::reference_type  reference_type;

    const auto td = get_tile_desc<value_type>(len);
    std::vector<value_type> result_cpu(td.num_tiles);
    hc::array<value_type> result(td.num_tiles);
    hc::array<value_type> scratch(len);

    tile_for<value_type>(td, [&,f,len,td](hc::tiled_index<1> t_idx, tile_buffer<value_type> buffer) [[hc]] 
    {
        const auto local = t_idx.local[0];
        const auto global = t_idx.global[0];
        const auto tile = t_idx.tile[0];

        // Join tile buffer elements
        const auto join = [&](std::size_t i, std::size_t j)
        {
            buffer.action_at(i, j, [&](value_type& x, const value_type& y)
            {
                ValueJoin::join(f, &x, &y);
            });
        };

        // Copy into tile
        buffer.action_at(local, [&](value_type& state)
        {
            ValueInit::init(f, &state);
            if (global < len) rocm_invoke<Tag>(f, transform_index(t_idx, td.tile_size, td.num_tiles), state, false);
        });
        t_idx.barrier.wait();
        // Up sweep phase
        for(std::size_t d=1;d<buffer.size();d*=2)
        {
            auto d2 = 2*d;
            auto i = local*d2;
            if(i<len)
            {
               auto j = i + d - 1;
               auto k = i + d2 - 1;

               ValueJoin::join(f, &buffer[k], &buffer[j]);
            }
        }
        t_idx.barrier.wait();

        result[tile] = buffer[buffer.size()-1];
        buffer[buffer.size()-1] = 0;
        // Down sweep phase
        for(std::size_t d=buffer.size()/2;d>0;d/=2)
        {
            auto d2 = 2*d;
            auto i = local*d2;
            if(i<len)
            {
               auto j = i + d - 1;
               auto k = i + d2 - 1;
               auto t = buffer[k];

               ValueJoin::join(f, &buffer[k], &buffer[j]);
               buffer[j] = t;
            }
            t_idx.barrier.wait();
        }
        // Copy tiles into global memory
        if (global < len) scratch[global] = buffer[local];
    }).wait();
    copy(result,result_cpu.data());

   for(int i=1; i<td.num_tiles; i++)
      ValueJoin::join(f, &result_cpu[i], &result_cpu[i-1]);

    copy(result_cpu.data(),result);
    size_t launch_len = (((len - 1) / td.tile_size) + 1) * td.tile_size;
    hc::parallel_for_each(hc::extent<1>(launch_len).tile(td.tile_size), [&,f,len,td](hc::tiled_index<1> t_idx) [[hc]] 
    {
        const auto global = t_idx.global[0];
        const auto tile = t_idx.tile[0];

        if (global < len) 
        {
            auto final_state = scratch[global];

            if (tile != 0) ValueJoin::join(f, &final_state, &result[tile-1]);
            rocm_invoke<Tag>(f, transform_index(t_idx, td.tile_size, td.num_tiles), final_state, true);
        }
    }).wait();
}

template< class Tag, class ReturnType, class F, class TransformIndex>
void scan_enqueue(
  const int len,
  const F & f,
  ReturnType & return_val,
  TransformIndex transform_index)
{
    typedef Kokkos::Impl::FunctorValueTraits< F, Tag>  ValueTraits;
    typedef Kokkos::Impl::FunctorValueInit<   F, Tag>  ValueInit;
    typedef Kokkos::Impl::FunctorValueJoin<   F, Tag>  ValueJoin;
    typedef Kokkos::Impl::FunctorValueOps<    F, Tag>  ValueOps;

    typedef typename ValueTraits::value_type    value_type;
    typedef typename ValueTraits::pointer_type    pointer_type;
    typedef typename ValueTraits::reference_type  reference_type;

    const auto td = get_tile_desc<value_type>(len);
    std::vector<value_type> result_cpu(td.num_tiles);
    hc::array<value_type> result(td.num_tiles);
    hc::array<value_type> scratch(len);
    std::vector<ReturnType> total_cpu(1);
    hc::array<ReturnType> total(1);

    tile_for<value_type>(td, [&,f,len,td](hc::tiled_index<1> t_idx, tile_buffer<value_type> buffer) [[hc]] 
    {
        const auto local = t_idx.local[0];
        const auto global = t_idx.global[0];
        const auto tile = t_idx.tile[0];

        // Join tile buffer elements
        const auto join = [&](std::size_t i, std::size_t j)
        {
            buffer.action_at(i, j, [&](value_type& x, const value_type& y)
            {
                ValueJoin::join(f, &x, &y);
            });
        };

        // Copy into tile
        buffer.action_at(local, [&](value_type& state)
        {
            ValueInit::init(f, &state);
            if (global < len) rocm_invoke<Tag>(f, transform_index(t_idx, td.tile_size, td.num_tiles), state, false);
        });
        t_idx.barrier.wait();
        // Up sweep phase
        for(std::size_t d=1;d<buffer.size();d*=2)
        {
            auto d2 = 2*d;
            auto i = local*d2;
            if(i<len)
            {
               auto j = i + d - 1;
               auto k = i + d2 - 1;
               ValueJoin::join(f, &buffer[k], &buffer[j]);
            }
        }
        t_idx.barrier.wait();

        result[tile] = buffer[buffer.size()-1];
        buffer[buffer.size()-1] = 0;
        // Down sweep phase
        for(std::size_t d=buffer.size()/2;d>0;d/=2)
        {
            auto d2 = 2*d;
            auto i = local*d2;
            if(i<len)
            {
               auto j = i + d - 1;
               auto k = i + d2 - 1;
               auto t = buffer[k];
               ValueJoin::join(f, &buffer[k], &buffer[j]);
               buffer[j] = t;
            }
            t_idx.barrier.wait();
        }
        // Copy tiles into global memory
        if (global < len) scratch[global] = buffer[local];
    }).wait();
    copy(result,result_cpu.data());

   for(int i=1; i<td.num_tiles; i++)
      ValueJoin::join(f, &result_cpu[i], &result_cpu[i-1]);

    copy(result_cpu.data(),result);
    size_t launch_len = (((len - 1) / td.tile_size) + 1) * td.tile_size;
    hc::parallel_for_each(hc::extent<1>(launch_len).tile(td.tile_size), [&,f,len,td](hc::tiled_index<1> t_idx) [[hc]] 
    {
        const auto global = t_idx.global[0];
        const auto tile = t_idx.tile[0];

        if (global < len) 
        {
            auto final_state = scratch[global];

            if (tile != 0) ValueJoin::join(f, &final_state, &result[tile-1]);
            rocm_invoke<Tag>(f, transform_index(t_idx, td.tile_size, td.num_tiles), final_state, true);
            if(global==(len-1))  total[0] = final_state;
        }
    }).wait();
    copy(total,total_cpu.data());
    return_val = total_cpu[0];
}

} // namespace Impl
} // namespace Kokkos
