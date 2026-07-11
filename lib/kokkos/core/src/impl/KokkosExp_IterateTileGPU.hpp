// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_EXP_ITERATE_TILE_GPU_HPP
#define KOKKOS_EXP_ITERATE_TILE_GPU_HPP

#include <Kokkos_Macros.hpp>

#include <algorithm>

#include <utility>

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_SYCL
template <typename index_type>
struct EmulateCUDADim3 {
  index_type x;
  index_type y;
  index_type z;
};
#endif

template <class Tag, class Functor, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION std::enable_if_t<std::is_void_v<Tag>>
_tag_invoke(Functor const& f, Args&&... args) {
  f((Args&&)args...);
}

template <class Tag, class Functor, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION std::enable_if_t<!std::is_void_v<Tag>>
_tag_invoke(Functor const& f, Args&&... args) {
  f(Tag{}, (Args&&)args...);
}

template <class Tag, class Functor, class T, size_t N, size_t... Idxs,
          class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION void _tag_invoke_array_helper(
    Functor const& f, T (&vals)[N], std::integer_sequence<size_t, Idxs...>,
    Args&&... args) {
  _tag_invoke<Tag>(f, vals[Idxs]..., (Args&&)args...);
}

template <class Tag, class Functor, class T, size_t N, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION void _tag_invoke_array(Functor const& f,
                                                        T (&vals)[N],
                                                        Args&&... args) {
  _tag_invoke_array_helper<Tag>(f, vals, std::make_index_sequence<N>{},
                                (Args&&)args...);
}

// ------------------------------------------------------------------------- //
// Compute GPU launch parameters (grid/block dimensions) for MDRangePolicy
//
// Ranks 2-3: Direct mapping - each policy dimension maps to one GPU dimension.
// Ranks 4-6: Dimension packing - pairs of policy dimensions are packed
//            into single GPU dimensions to fit the 3D hardware limit.
//
// Returns: CUDA/HIP: std::pair<dim3 grid, dim3 block>
//          SYCL:     sycl::nd_range<3>{global, local}
//
template <typename... Traits, typename MaxGridSize>
auto compute_device_launch_params(
    const Kokkos::MDRangePolicy<Traits...>& policy,
    const MaxGridSize& max_grid_size) {
  using Policy           = Kokkos::MDRangePolicy<Traits...>;
  using array_index_type = typename Policy::array_index_type;

#ifdef KOKKOS_ENABLE_SYCL
  EmulateCUDADim3 block{1, 1, 1};
#else
  dim3 block{1, 1, 1};
#endif

  array_index_type grid_0 = 1;
  array_index_type grid_1 = 1;
  array_index_type grid_2 = 1;

  if constexpr (Policy::inner_direction == Iterate::Left) {
    if constexpr (Policy::rank == 2) {
      block.x = policy.m_tile[0];
      block.y = policy.m_tile[1];
      grid_0  = policy.m_tile_end[0];
      grid_1  = policy.m_tile_end[1];
    } else if constexpr (Policy::rank == 3) {
      block.x = policy.m_tile[0];
      block.y = policy.m_tile[1];
      block.z = policy.m_tile[2];
      grid_0  = policy.m_tile_end[0];
      grid_1  = policy.m_tile_end[1];
      grid_2  = policy.m_tile_end[2];
    } else if constexpr (Policy::rank == 4) {
      block.x = policy.m_tile[0] * policy.m_tile[1];
      block.y = policy.m_tile[2];
      block.z = policy.m_tile[3];
      grid_0  = policy.m_tile_end[0] * policy.m_tile_end[1];
      grid_1  = policy.m_tile_end[2];
      grid_2  = policy.m_tile_end[3];
    } else if constexpr (Policy::rank == 5) {
      block.x = policy.m_tile[0] * policy.m_tile[1];
      block.y = policy.m_tile[2] * policy.m_tile[3];
      block.z = policy.m_tile[4];
      grid_0  = policy.m_tile_end[0] * policy.m_tile_end[1];
      grid_1  = policy.m_tile_end[2] * policy.m_tile_end[3];
      grid_2  = policy.m_tile_end[4];
    } else if constexpr (Policy::rank == 6) {
      block.x = policy.m_tile[0] * policy.m_tile[1];
      block.y = policy.m_tile[2] * policy.m_tile[3];
      block.z = policy.m_tile[4] * policy.m_tile[5];
      grid_0  = policy.m_tile_end[0] * policy.m_tile_end[1];
      grid_1  = policy.m_tile_end[2] * policy.m_tile_end[3];
      grid_2  = policy.m_tile_end[4] * policy.m_tile_end[5];
    }
  } else {  // InnerDirection == Right
    if constexpr (Policy::rank == 2) {
      block.x = policy.m_tile[1];
      block.y = policy.m_tile[0];
      grid_0  = policy.m_tile_end[1];
      grid_1  = policy.m_tile_end[0];
    } else if constexpr (Policy::rank == 3) {
      block.x = policy.m_tile[2];
      block.y = policy.m_tile[1];
      block.z = policy.m_tile[0];
      grid_0  = policy.m_tile_end[2];
      grid_1  = policy.m_tile_end[1];
      grid_2  = policy.m_tile_end[0];
    } else if constexpr (Policy::rank == 4) {
      block.x = policy.m_tile[3] * policy.m_tile[2];
      block.y = policy.m_tile[1];
      block.z = policy.m_tile[0];
      grid_0  = policy.m_tile_end[3] * policy.m_tile_end[2];
      grid_1  = policy.m_tile_end[1];
      grid_2  = policy.m_tile_end[0];
    } else if constexpr (Policy::rank == 5) {
      block.x = policy.m_tile[4] * policy.m_tile[3];
      block.y = policy.m_tile[2] * policy.m_tile[1];
      block.z = policy.m_tile[0];
      grid_0  = policy.m_tile_end[4] * policy.m_tile_end[3];
      grid_1  = policy.m_tile_end[2] * policy.m_tile_end[1];
      grid_2  = policy.m_tile_end[0];
    } else if constexpr (Policy::rank == 6) {
      block.x = policy.m_tile[5] * policy.m_tile[4];
      block.y = policy.m_tile[3] * policy.m_tile[2];
      block.z = policy.m_tile[1] * policy.m_tile[0];
      grid_0  = policy.m_tile_end[5] * policy.m_tile_end[4];
      grid_1  = policy.m_tile_end[3] * policy.m_tile_end[2];
      grid_2  = policy.m_tile_end[1] * policy.m_tile_end[0];
    }
  }

#ifdef KOKKOS_ENABLE_SYCL
  // SYCL uses nd_range with global = grid * local sizes
  sycl::range<3> local_sizes(block.x, block.y, block.z);
  sycl::range<3> global_sizes(
      std::min<array_index_type>(grid_0, max_grid_size[0]) * local_sizes[0],
      std::min<array_index_type>(grid_1, max_grid_size[1]) * local_sizes[1],
      std::min<array_index_type>(grid_2, max_grid_size[2]) * local_sizes[2]);
  return sycl::nd_range<3>(global_sizes, local_sizes);
#else
  dim3 grid(std::min<array_index_type>(grid_0, max_grid_size[0]),
            std::min<array_index_type>(grid_1, max_grid_size[1]),
            std::min<array_index_type>(grid_2, max_grid_size[2]));
  return std::pair(grid, block);
#endif
}

// ------------------------------------------------------------------------- //
// ParallelFor iteration pattern - maps GPU threads to N-D iteration space
//
// For ranks 2-3: Direct mapping of hardware threads to iteration space
// dimensions.
// For ranks 4-6: Multiple logical indices are packed into single
// hardware dimensions.
//
// 1. Start iterating at the hardware thread identifier.
// 2. Extend the iteration space range with stride loops using grid dimensions.
// 3. Bounds check against m_upper to filter out-of-bounds iterations.
//
template <int Rank, typename array_index_type, typename index_type,
          typename Functor, Kokkos::Iterate Layout, typename Tag>
struct DeviceIterate;

template <int Rank, typename array_index_type, typename index_type,
          typename Functor, Kokkos::Iterate Layout, typename Tag>
struct DeviceIterate {
  using array_type = Kokkos::Array<array_index_type, Rank>;

 private:
  const array_type m_lower;
  const array_type m_upper;
  const array_type m_extent;  // tile_size * num_tiles
  const Functor& m_functor;

#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif

 public:
#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterate(
      const array_type& lower, const array_type& upper,
      const array_type& extent, const Functor& functor,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_lower(lower),
        m_upper(upper),
        m_extent(extent),
        m_functor(functor),
        gridDim(gridDim_),
        blockDim(blockDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else

  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterate(const array_type& lower,
                                            const array_type& upper,
                                            const array_type& extent,
                                            const Functor& functor)
      : m_lower(lower), m_upper(upper), m_extent(extent), m_functor(functor) {}
#endif

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const { iterate(std::integral_constant<unsigned, Rank>()); }

 private:
  // Runtime expression to determine if Dim is part of a packed pair
  // Packing occurs on consecutive dimension pairs for rank > 3
  template <unsigned Dim>
  KOKKOS_IMPL_DEVICE_FUNCTION static consteval bool is_packed_index() {
    return ((Dim == 0 || Dim == 1) && Rank > 3) ||
           ((Dim == 2 || Dim == 3) && Rank > 4) ||
           ((Dim == 4 || Dim == 5) && Rank > 5);
  }

  // Packed: returns flat hardware thread ID (unpacking happens in iterate())
  // Unpacked: returns global index (lower + blockIdx * blockDim + threadIdx)
  template <unsigned R>
  KOKKOS_IMPL_DEVICE_FUNCTION KOKKOS_IMPL_FORCEINLINE constexpr index_type
  my_begin() const noexcept {
    static_assert(R < 6);
    if constexpr (is_packed_index<R>()) {
      if constexpr (R == 0 || R == 1) {
        return blockIdx.x * blockDim.x + threadIdx.x;
      } else if constexpr (R == 2 || R == 3) {
        return blockIdx.y * blockDim.y + threadIdx.y;
      } else if constexpr (R == 4 || R == 5) {
        return blockIdx.z * blockDim.z + threadIdx.z;
      }
    } else {
      // No packed index
      if constexpr (Rank < 4) {
        if constexpr (R == 0) {
          return m_lower[R] + blockIdx.x * blockDim.x + threadIdx.x;
        } else if constexpr (R == 1) {
          return m_lower[R] + blockIdx.y * blockDim.y + threadIdx.y;
        } else if constexpr (R == 2) {
          return m_lower[R] + blockIdx.z * blockDim.z + threadIdx.z;
        }
      } else {
        // Mix of packed and unpacked for Rank 4 and 5
        if constexpr (R == 2) {
          return m_lower[R] + blockIdx.y * blockDim.y + threadIdx.y;
        } else if constexpr (R == 3 || R == 4) {
          return m_lower[R] + blockIdx.z * blockDim.z + threadIdx.z;
        }
      }
    }
    return m_lower[R];
  }

  // Packed: end at the product of two consecutive extents
  // Unpacked: directly use m_upper
  template <unsigned R>
  KOKKOS_IMPL_DEVICE_FUNCTION KOKKOS_IMPL_FORCEINLINE constexpr index_type
  my_end() const noexcept {
    static_assert(R < 6);
    if constexpr (is_packed_index<R>()) {
      if constexpr (R % 2 == 0) {
        return m_extent[R] * m_extent[R + 1];
      } else {
        return m_extent[R] * m_extent[R - 1];
      }
    } else {
      return m_upper[R];
    }
  }

  // Stride by the total number of threads in the GPU dimension
  template <unsigned R>
  KOKKOS_IMPL_DEVICE_FUNCTION KOKKOS_IMPL_FORCEINLINE constexpr index_type
  my_stride() const noexcept {
    static_assert(R < 6);
    if constexpr (is_packed_index<R>()) {
      if constexpr (R == 0 || R == 1) {
        return static_cast<index_type>(blockDim.x) *
               static_cast<index_type>(gridDim.x);
      } else if constexpr (R == 2 || R == 3) {
        return static_cast<index_type>(blockDim.y) *
               static_cast<index_type>(gridDim.y);
      } else if constexpr (R == 4 || R == 5) {
        return static_cast<index_type>(blockDim.z) *
               static_cast<index_type>(gridDim.z);
      }
    } else {
      // No packed index for all ranks
      if constexpr (Rank < 4) {
        if constexpr (R == 0) {
          return static_cast<index_type>(blockDim.x) *
                 static_cast<index_type>(gridDim.x);
        } else if constexpr (R == 1) {
          return static_cast<index_type>(blockDim.y) *
                 static_cast<index_type>(gridDim.y);
        } else if constexpr (R == 2) {
          return static_cast<index_type>(blockDim.z) *
                 static_cast<index_type>(gridDim.z);
        }
      } else {
        // Mix of packed and unpacked for Rank 4 and 5
        if constexpr (R == 2) {
          return static_cast<index_type>(blockDim.y) *
                 static_cast<index_type>(gridDim.y);
        } else if constexpr (R == 3 || R == 4) {
          return static_cast<index_type>(blockDim.z) *
                 static_cast<index_type>(gridDim.z);
        }
      }
    }
    return index_type{1};
  }

  // ----------------------------------------------------------------------- //
  // Nested loops with recursive template instantiation
  //
  // Accumulates indices in parameter pack Idxs...
  // The fastest changing index is always i0 (innermost loop).
  //
  // Functor call order depends on the Layout:
  //  Layout::Left:
  //    functor(i0, i1, i2, ..., iR)
  //  Layout::Right:
  //    functor(iR, ..., i2, i1, i0)
  //
  // For Layout::Right, bounds were previously swapped during ParallelFor
  // construction, so i0 correctly iterates over the range of iR while
  // remaining the "fastest-changing" index.
  //
  template <unsigned R, typename... Idxs>
  KOKKOS_IMPL_DEVICE_FUNCTION inline void iterate(
      std::integral_constant<unsigned, R>, Idxs... idxs) const {
    constexpr unsigned rankIdx = R - 1;
    const index_type start     = my_begin<rankIdx>();
    const index_type end       = my_end<rankIdx>();
    const index_type stride    = my_stride<rankIdx>();

    for (index_type idx = start; idx < end; idx += stride) {
      if constexpr (is_packed_index<rankIdx>()) {
        static_assert(R >= 2);
        // Unpack two consecutive indices
        constexpr unsigned idx1 = (rankIdx % 2 == 0) ? rankIdx : (rankIdx - 1);
        constexpr unsigned idx2 = (rankIdx % 2 == 0) ? (rankIdx + 1) : rankIdx;

        const index_type id_1 = idx % m_extent[idx1] + m_lower[idx1];
        const index_type id_2 = idx / m_extent[idx1] + m_lower[idx2];

        if (id_1 < m_upper[idx1] && id_2 < m_upper[idx2]) {
          if constexpr (Layout == Iterate::Left) {
            iterate(std::integral_constant<unsigned, R - 2>(), id_1, id_2,
                    idxs...);
          } else {
            iterate(std::integral_constant<unsigned, R - 2>(), idxs..., id_2,
                    id_1);
          }
        }
      } else {
        if constexpr (Layout == Iterate::Left) {
          iterate(std::integral_constant<unsigned, R - 1>(), idx, idxs...);
        } else {
          iterate(std::integral_constant<unsigned, R - 1>(), idxs..., idx);
        }
      }
    }
  }

  template <typename... Idxs>
  KOKKOS_IMPL_DEVICE_FUNCTION inline void iterate(
      std::integral_constant<unsigned, 0u>, Idxs... idxs) const {
    Impl::_tag_invoke<Tag>(m_functor, idxs...);
  }
};

// ----------------------------------------------------------------------------------

namespace Reduce {

template <typename T>
struct is_array_type : std::false_type {
  using value_type = T;
};

template <typename T>
struct is_array_type<T*> : std::true_type {
  using value_type = T;
};

template <typename T>
struct is_array_type<T[]> : std::true_type {
  using value_type = T;
};

// ------------------------------------------------------------------ //

template <typename T>
using value_type_storage_t =
    std::conditional_t<is_array_type<T>::value, std::decay_t<T>,
                       std::add_lvalue_reference_t<T>>;

// ParallelReduce iteration pattern
// Scalar reductions

// num_blocks = min( num_tiles, max_num_blocks ); //i.e. determined by number of
// tiles and reduction algorithm constraints extract n-dim tile offsets (i.e.
// tile's global starting mulit-index) from the tileid = blockid using tile
// dimensions local indices within a tile extracted from (index_type)threadIdx_x
// using tile dims, constrained by blocksize combine tile and local id info for
// multi-dim global ids

// Pattern:
// Each block+thread is responsible for a tile+local_id combo (additional when
// striding by num_blocks)
// 1. create offset arrays
// 2. loop over number of tiles, striding by griddim (equal to num tiles, or max
// num blocks)
// 3. temps set for tile_idx and thrd_idx, which will be modified
// 4. if LL vs LR:
//      determine tile starting point offsets (multidim)
//      determine local index offsets (multidim)
//      concatentate tile offset + local offset for global multi-dim index
//    if offset withinin range bounds AND local offset within tile bounds, call
//    functor

template <int N, typename PolicyType, typename Functor, typename Tag,
          typename ValueType, typename Enable = void>
struct DeviceIterateTile {
  using index_type         = typename PolicyType::index_type;
  using value_type_storage = value_type_storage_t<ValueType>;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_, value_type_storage v_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        m_v(v_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_,
                                                value_type_storage v_)
      : m_policy(policy_), m_func(f_), m_v(v_) {}
#endif

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    if (static_cast<index_type>(blockIdx.x) < m_policy.m_num_tiles &&
        static_cast<index_type>(threadIdx.y) < m_policy.m_prod_tile_dims) {
      index_type m_offset[PolicyType::rank];  // tile starting global id offset
      index_type
          m_local_offset[PolicyType::rank];  // tile starting global id offset

      for (index_type tileidx = static_cast<index_type>(blockIdx.x);
           tileidx < m_policy.m_num_tiles; tileidx += gridDim.x) {
        index_type tile_idx =
            tileidx;  // temp because tile_idx will be modified while
                      // determining tile starting point offsets
        index_type thrd_idx = static_cast<index_type>(threadIdx.y);
        bool in_bounds      = true;

        // LL
        if (PolicyType::inner_direction == Iterate::Left) {
          for (int i = 0; i < PolicyType::rank; ++i) {
            m_offset[i] =
                (tile_idx % m_policy.m_tile_end[i]) * m_policy.m_tile[i] +
                m_policy.m_lower[i];
            tile_idx /= m_policy.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_policy.m_tile[i]);
            thrd_idx /= m_policy.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if (!(m_offset[i] < m_policy.m_upper[i] &&
                  m_local_offset[i] < m_policy.m_tile[i])) {
              in_bounds = false;
            }
          }
          if (in_bounds) {
            Impl::_tag_invoke_array<Tag>(m_func, m_offset, m_v);
          }
        }
        // LR
        else {
          for (int i = PolicyType::rank - 1; i >= 0; --i) {
            m_offset[i] =
                (tile_idx % m_policy.m_tile_end[i]) * m_policy.m_tile[i] +
                m_policy.m_lower[i];
            tile_idx /= m_policy.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] =
                (thrd_idx %
                 m_policy.m_tile[i]);  // Move this to first computation,
                                       // add to m_offset right away
            thrd_idx /= m_policy.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if (!(m_offset[i] < m_policy.m_upper[i] &&
                  m_local_offset[i] < m_policy.m_tile[i])) {
              in_bounds = false;
            }
          }
          if (in_bounds) {
            Impl::_tag_invoke_array<Tag>(m_func, m_offset, m_v);
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
  value_type_storage m_v;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

}  // namespace Reduce
}  // namespace Impl
}  // namespace Kokkos
#endif
