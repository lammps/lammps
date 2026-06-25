/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   ILVES constraint solver: append-only growing memory pool and the STL
   allocator that draws from it.  Ported from GROMACS 2021 ILVES (LGPL-2.1),
   src/gromacs/mdlib/growing_mem_pool.h + growing_allocator.h, combined into
   one header.  See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_MEMPOOL_H
#define LMP_ILVES_MEMPOOL_H

#include <algorithm>
#include <cstddef>
#include <list>
#include <utility>

namespace LAMMPS_NS {
namespace ILVES {

/**
 * A memory pool that only creates new chunks, it does not allow deallocation.
 */
class GrowingMemPool {
 public:
  /**
   * Construct an empty memory pool. When requested, the pool will
   * allocate at least min_chunk_elements elements.
   *
   * @param min_chunk_elements The minimum number of elements that
   * the pool will allocate when requested.
   */
  GrowingMemPool(std::size_t min_chunk_elements) noexcept :
      min_chunk_elements(min_chunk_elements){};

  /**
   * Deallocate all the memory allocated.
   */
  ~GrowingMemPool() noexcept
  {
    for (auto &chunk : chunks) { ::operator delete(static_cast<void *>(chunk.data)); }
  };

  /**
   * Allocate memory for num_elements elements of size element_bytes.
   * Use previously allocated unused bytes if possible.
   *
   * @param num_elements The number of elements to allocate.
   * @param element_bytes The size of each element.
   * @return void* A pointer to the allocated memory.
   */
  void *allocate(std::size_t num_elements, std::size_t element_bytes)
  {
    const std::size_t new_bytes = num_elements * element_bytes;

    const bool new_chunk =
        chunks.empty() || chunks.back().total_bytes - chunks.back().used_bytes < new_bytes;

    void *ret = nullptr;

    if (new_chunk) {
      const std::size_t min_chunk_bytes = min_chunk_elements * element_bytes;
      const std::size_t chunk_bytes = std::max(min_chunk_bytes, new_bytes);
      ret = ::operator new(chunk_bytes);
      chunks.push_back({chunk_bytes, new_bytes, static_cast<std::byte *>(ret)});
    } else {
      ret = static_cast<void *>(chunks.back().data + chunks.back().used_bytes);
      chunks.back().used_bytes += new_bytes;
    }

    return ret;
  }

 private:
  struct Chunk {
    std::size_t total_bytes;
    std::size_t used_bytes;
    std::byte *data;
  };
  std::list<Chunk> chunks;

  std::size_t min_chunk_elements;
};

/**
 * STL allocator that uses a GrowingMemPool.
 *
 * An STL allocator must have no state, so we cannot store the memory
 * pool in the allocator. Instead, we store a pointer to the memory pool.
 */
template <class T> class GrowingAllocator {
 public:
  typedef T value_type;

  /**
   * Construct an allocator that uses the given memory pool.
   */
  GrowingAllocator(GrowingMemPool *mem_pool) noexcept : mem_pool(mem_pool){};

  GrowingAllocator(const GrowingAllocator &other) noexcept = default;

  GrowingAllocator &operator=(const GrowingAllocator &other) noexcept
  {
    if (this == &other) { return *this; }

    mem_pool = other.mem_pool;
    return *this;
  };

  /**
   * Create an allocator of type U that uses the same memory pool.
   */
  template <class U> GrowingAllocator(const GrowingAllocator<U> &other) noexcept :
      mem_pool(other.mem_pool){};

  GrowingAllocator(GrowingAllocator &&other) noexcept = default;

  GrowingAllocator &operator=(GrowingAllocator &&other) noexcept = default;

  /**
   * Move operator from an allocator of type U.
   */
  template <class U> GrowingAllocator(GrowingAllocator<U> &&other) noexcept :
      mem_pool(std::move(other.mem_pool))
  {
    other.mem_pool = nullptr;
  };

  /**
   * Allocate memory for num elements of type value_type.
   */
  value_type *allocate(std::size_t num)
  {
    return static_cast<value_type *>(mem_pool->allocate(num, sizeof(value_type)));
  }

  /**
   * We do not do anything, since a GrowingAllocator just grows.
   */
  void deallocate(value_type *, std::size_t)
  {
    // We don't do anything, we just grow.
  }

  template <class U> bool operator==(const GrowingAllocator<U> &) const { return true; }

  template <class U> bool operator!=(const GrowingAllocator<U> &) const { return false; }

 private:
  GrowingMemPool *mem_pool;

  // Add a friend declaration for the rebind constructor
  template <class U> friend class GrowingAllocator;
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif
