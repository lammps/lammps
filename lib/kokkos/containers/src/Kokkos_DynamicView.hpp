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

#ifndef KOKKOS_DYNAMIC_VIEW_HPP
#define KOKKOS_DYNAMIC_VIEW_HPP

#include <cstdio>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Experimental {

// Simple metafunction for choosing memory space
// In the current implementation, if memory_space == CudaSpace,
// use CudaUVMSpace for the chunk 'array' allocation, which
// contains will contain pointers to chunks of memory allocated
// in CudaSpace
namespace Impl {
template <class MemSpace>
struct ChunkArraySpace {
  using memory_space = MemSpace;
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ChunkArraySpace<Kokkos::CudaSpace> {
  using memory_space = typename Kokkos::CudaUVMSpace;
};
#endif
#ifdef KOKKOS_ENABLE_HIP
template <>
struct ChunkArraySpace<Kokkos::Experimental::HIPSpace> {
  using memory_space = typename Kokkos::Experimental::HIPHostPinnedSpace;
};
#endif
}  // end namespace Impl

/** \brief Dynamic views are restricted to rank-one and no layout.
 *         Resize only occurs on host outside of parallel_regions.
 *         Subviews are not allowed.
 */
template <typename DataType, typename... P>
class DynamicView : public Kokkos::ViewTraits<DataType, P...> {
 public:
  using traits = Kokkos::ViewTraits<DataType, P...>;

 private:
  template <class, class...>
  friend class DynamicView;

  using track_type = Kokkos::Impl::SharedAllocationTracker;

  static_assert(traits::rank == 1 && traits::rank_dynamic == 1,
                "DynamicView must be rank-one");

  // It is assumed that the value_type is trivially copyable;
  // when this is not the case, potential problems can occur.
  static_assert(std::is_same<typename traits::specialize, void>::value,
                "DynamicView only implemented for non-specialized View type");

  template <class Space, bool = Kokkos::Impl::MemorySpaceAccess<
                             Space, typename traits::memory_space>::accessible>
  struct verify_space {
    KOKKOS_FORCEINLINE_FUNCTION static void check() {}
  };

  template <class Space>
  struct verify_space<Space, false> {
    KOKKOS_FORCEINLINE_FUNCTION static void check() {
      Kokkos::abort(
          "Kokkos::DynamicView ERROR: attempt to access inaccessible memory "
          "space");
    };
  };

 private:
  track_type m_track;
  typename traits::value_type** m_chunks =
      nullptr;             // array of pointers to 'chunks' of memory
  unsigned m_chunk_shift;  // ceil(log2(m_chunk_size))
  unsigned m_chunk_mask;   // m_chunk_size - 1
  unsigned m_chunk_max;  // number of entries in the chunk array - each pointing
                         // to a chunk of extent == m_chunk_size entries
  unsigned m_chunk_size;  // 2 << (m_chunk_shift - 1)

 public:
  //----------------------------------------------------------------------

  /** \brief  Compatible view of array of scalar types */
  using array_type =
      DynamicView<typename traits::data_type, typename traits::device_type>;

  /** \brief  Compatible view of const data type */
  using const_type = DynamicView<typename traits::const_data_type,
                                 typename traits::device_type>;

  /** \brief  Compatible view of non-const data type */
  using non_const_type = DynamicView<typename traits::non_const_data_type,
                                     typename traits::device_type>;

  /** \brief  Must be accessible everywhere */
  using HostMirror = DynamicView;

  /** \brief Unified types */
  using uniform_device =
      Kokkos::Device<typename traits::device_type::execution_space,
                     Kokkos::AnonymousSpace>;
  using uniform_type               = array_type;
  using uniform_const_type         = const_type;
  using uniform_runtime_type       = array_type;
  using uniform_runtime_const_type = const_type;
  using uniform_nomemspace_type =
      DynamicView<typename traits::data_type, uniform_device>;
  using uniform_const_nomemspace_type =
      DynamicView<typename traits::const_data_type, uniform_device>;
  using uniform_runtime_nomemspace_type =
      DynamicView<typename traits::data_type, uniform_device>;
  using uniform_runtime_const_nomemspace_type =
      DynamicView<typename traits::const_data_type, uniform_device>;

  //----------------------------------------------------------------------

  enum { Rank = 1 };

  KOKKOS_INLINE_FUNCTION
  size_t allocation_extent() const noexcept {
    uintptr_t n = *reinterpret_cast<const uintptr_t*>(m_chunks + m_chunk_max);
    return (n << m_chunk_shift);
  }

  KOKKOS_INLINE_FUNCTION
  size_t chunk_size() const noexcept { return m_chunk_size; }

  KOKKOS_INLINE_FUNCTION
  size_t size() const noexcept {
    size_t extent_0 =
        *reinterpret_cast<const size_t*>(m_chunks + m_chunk_max + 1);
    return extent_0;
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION size_t extent(const iType& r) const {
    return r == 0 ? size() : 1;
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION size_t extent_int(const iType& r) const {
    return r == 0 ? size() : 1;
  }

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return 0; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    *s = 0;
  }

  //----------------------------------------
  // Allocation tracking properties

  KOKKOS_INLINE_FUNCTION
  int use_count() const { return m_track.use_count(); }

  inline const std::string label() const {
    return m_track.template get_label<typename traits::memory_space>();
  }

  //----------------------------------------------------------------------
  // Range span is the span which contains all members.

  using reference_type = typename traits::value_type&;
  using pointer_type   = typename traits::value_type*;

  enum {
    reference_type_is_lvalue_reference =
        std::is_lvalue_reference<reference_type>::value
  };

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return false;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const { return 0; }

  //----------------------------------------

  template <typename I0, class... Args>
  KOKKOS_INLINE_FUNCTION reference_type
  operator()(const I0& i0, const Args&... /*args*/) const {
    static_assert(Kokkos::Impl::are_integral<I0, Args...>::value,
                  "Indices must be integral type");

    DynamicView::template verify_space<
        Kokkos::Impl::ActiveExecutionMemorySpace>::check();

    // Which chunk is being indexed.
    const uintptr_t ic = uintptr_t(i0 >> m_chunk_shift);

    typename traits::value_type* volatile* const ch = m_chunks + ic;

    // Do bounds checking if enabled or if the chunk pointer is zero.
    // If not bounds checking then we assume a non-zero pointer is valid.

#if !defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    if (nullptr == *ch)
#endif
    {
      // Verify that allocation of the requested chunk in in progress.

      // The allocated chunk counter is m_chunks[ m_chunk_max ]
      const uintptr_t n =
          *reinterpret_cast<uintptr_t volatile*>(m_chunks + m_chunk_max);

      if (n <= ic) {
        Kokkos::abort("Kokkos::DynamicView array bounds error");
      }

      // Allocation of this chunk is in progress
      // so wait for allocation to complete.
      while (nullptr == *ch)
        ;
    }

    return (*ch)[i0 & m_chunk_mask];
  }

  //----------------------------------------
  /** \brief  Resizing in serial can grow or shrink the array size
   *          up to the maximum number of chunks
   * */
  template <typename IntType>
  inline typename std::enable_if<
      std::is_integral<IntType>::value &&
      Kokkos::Impl::MemorySpaceAccess<
          Kokkos::HostSpace,
          typename Impl::ChunkArraySpace<
              typename traits::memory_space>::memory_space>::accessible>::type
  resize_serial(IntType const& n) {
    using local_value_type   = typename traits::value_type;
    using value_pointer_type = local_value_type*;

    const uintptr_t NC =
        (n + m_chunk_mask) >>
        m_chunk_shift;  // New total number of chunks needed for resize

    if (m_chunk_max < NC) {
      Kokkos::abort("DynamicView::resize_serial exceeded maximum size");
    }

    // *m_chunks[m_chunk_max] stores the current number of chunks being used
    uintptr_t* const pc = reinterpret_cast<uintptr_t*>(m_chunks + m_chunk_max);
    std::string _label =
        m_track.template get_label<typename traits::memory_space>();
    if (*pc < NC) {
      while (*pc < NC) {
        m_chunks[*pc] = reinterpret_cast<value_pointer_type>(
            typename traits::memory_space().allocate(
                _label.c_str(), sizeof(local_value_type) << m_chunk_shift));
        ++*pc;
      }
    } else {
      while (NC + 1 <= *pc) {
        --*pc;
        typename traits::memory_space().deallocate(
            _label.c_str(), m_chunks[*pc],
            sizeof(local_value_type) << m_chunk_shift);
        m_chunks[*pc] = nullptr;
      }
    }
    // *m_chunks[m_chunk_max+1] stores the 'extent' requested by resize
    *(pc + 1) = n;
  }

  KOKKOS_INLINE_FUNCTION bool is_allocated() const {
    if (m_chunks == nullptr) {
      return false;
    } else {
      // *m_chunks[m_chunk_max] stores the current number of chunks being used
      uintptr_t* const pc =
          reinterpret_cast<uintptr_t*>(m_chunks + m_chunk_max);
      return (*(pc + 1) > 0);
    }
  }

  //----------------------------------------------------------------------

  ~DynamicView()                  = default;
  DynamicView()                   = default;
  DynamicView(DynamicView&&)      = default;
  DynamicView(const DynamicView&) = default;
  DynamicView& operator=(DynamicView&&) = default;
  DynamicView& operator=(const DynamicView&) = default;

  template <class RT, class... RP>
  DynamicView(const DynamicView<RT, RP...>& rhs)
      : m_track(rhs.m_track),
        m_chunks((typename traits::value_type**)rhs.m_chunks),
        m_chunk_shift(rhs.m_chunk_shift),
        m_chunk_mask(rhs.m_chunk_mask),
        m_chunk_max(rhs.m_chunk_max),
        m_chunk_size(rhs.m_chunk_size) {
    using SrcTraits = typename DynamicView<RT, RP...>::traits;
    using Mapping   = Kokkos::Impl::ViewMapping<traits, SrcTraits, void>;
    static_assert(Mapping::is_assignable,
                  "Incompatible DynamicView copy construction");
  }

  //----------------------------------------------------------------------

  struct Destroy {
    using local_value_type = typename traits::value_type;
    std::string m_label;
    local_value_type** m_chunks;
    unsigned m_chunk_max;
    bool m_destroy;
    unsigned m_chunk_size;

    // Initialize or destroy array of chunk pointers.
    // Two entries beyond the max chunks are allocation counters.
    inline void operator()(unsigned i) const {
      if (m_destroy && i < m_chunk_max && nullptr != m_chunks[i]) {
        typename traits::memory_space().deallocate(
            m_label.c_str(), m_chunks[i],
            sizeof(local_value_type) * m_chunk_size);
      }
      m_chunks[i] = nullptr;
    }

    void execute(bool arg_destroy) {
      using Range = Kokkos::RangePolicy<typename HostSpace::execution_space>;

      m_destroy = arg_destroy;

      Kokkos::Impl::ParallelFor<Destroy, Range> closure(
          *this,
          Range(0, m_chunk_max + 2));  // Add 2 to 'destroy' extra slots storing
                                       // num_chunks and extent; previously + 1

      closure.execute();

      typename traits::execution_space().fence();
      // Impl::ChunkArraySpace< typename traits::memory_space
      // >::memory_space::execution_space().fence();
    }

    void construct_shared_allocation() { execute(false); }

    void destroy_shared_allocation() { execute(true); }

    Destroy()               = default;
    Destroy(Destroy&&)      = default;
    Destroy(const Destroy&) = default;
    Destroy& operator=(Destroy&&) = default;
    Destroy& operator=(const Destroy&) = default;

    Destroy(std::string label, typename traits::value_type** arg_chunk,
            const unsigned arg_chunk_max, const unsigned arg_chunk_size)
        : m_label(label),
          m_chunks(arg_chunk),
          m_chunk_max(arg_chunk_max),
          m_destroy(false),
          m_chunk_size(arg_chunk_size) {}
  };

  /**\brief  Allocation constructor
   *
   *  Memory is allocated in chunks
   *  A maximum size is required in order to allocate a
   *  chunk-pointer array.
   */
  explicit inline DynamicView(const std::string& arg_label,
                              const unsigned min_chunk_size,
                              const unsigned max_extent)
      : m_track(),
        m_chunks(nullptr)
        // The chunk size is guaranteed to be a power of two
        ,
        m_chunk_shift(Kokkos::Impl::integral_power_of_two_that_contains(
            min_chunk_size))  // div ceil(log2(min_chunk_size))
        ,
        m_chunk_mask((1 << m_chunk_shift) - 1)  // mod
        ,
        m_chunk_max((max_extent + m_chunk_mask) >>
                    m_chunk_shift)  // max num pointers-to-chunks in array
        ,
        m_chunk_size(2 << (m_chunk_shift - 1)) {
    using chunk_array_memory_space = typename Impl::ChunkArraySpace<
        typename traits::memory_space>::memory_space;
    // A functor to deallocate all of the chunks upon final destruction
    using record_type =
        Kokkos::Impl::SharedAllocationRecord<chunk_array_memory_space, Destroy>;

    // Allocate chunk pointers and allocation counter
    record_type* const record =
        record_type::allocate(chunk_array_memory_space(), arg_label,
                              (sizeof(pointer_type) * (m_chunk_max + 2)));
    // Allocate + 2 extra slots so that *m_chunk[m_chunk_max] ==
    // num_chunks_alloc and *m_chunk[m_chunk_max+1] == extent This must match in
    // Destroy's execute(...) method

    m_chunks = reinterpret_cast<pointer_type*>(record->data());

    record->m_destroy = Destroy(arg_label, m_chunks, m_chunk_max, m_chunk_size);

    // Initialize to zero
    record->m_destroy.construct_shared_allocation();

    m_track.assign_allocated_record_to_uninitialized(record);
  }
};

}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {

template <class T, class... P>
inline typename Kokkos::Experimental::DynamicView<T, P...>::HostMirror
create_mirror_view(const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return src;
}

template <class T, class... DP, class... SP>
inline void deep_copy(const View<T, DP...>& dst,
                      const Kokkos::Experimental::DynamicView<T, SP...>& src) {
  using dst_type = View<T, DP...>;
  using src_type = Kokkos::Experimental::DynamicView<T, SP...>;

  using dst_execution_space = typename ViewTraits<T, DP...>::execution_space;
  using src_memory_space    = typename ViewTraits<T, SP...>::memory_space;

  enum {
    DstExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<dst_execution_space,
                                         src_memory_space>::accessible
  };

  if (DstExecCanAccessSrc) {
    // Copying data between views in accessible memory spaces and either
    // non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap<dst_type, src_type>(dst, src);
  } else {
    Kokkos::Impl::throw_runtime_exception(
        "deep_copy given views that would require a temporary allocation");
  }
}

template <class T, class... DP, class... SP>
inline void deep_copy(const Kokkos::Experimental::DynamicView<T, DP...>& dst,
                      const View<T, SP...>& src) {
  using dst_type = Kokkos::Experimental::DynamicView<T, SP...>;
  using src_type = View<T, DP...>;

  using dst_execution_space = typename ViewTraits<T, DP...>::execution_space;
  using src_memory_space    = typename ViewTraits<T, SP...>::memory_space;

  enum {
    DstExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<dst_execution_space,
                                         src_memory_space>::accessible
  };

  if (DstExecCanAccessSrc) {
    // Copying data between views in accessible memory spaces and either
    // non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap<dst_type, src_type>(dst, src);
  } else {
    Kokkos::Impl::throw_runtime_exception(
        "deep_copy given views that would require a temporary allocation");
  }
}

namespace Impl {
template <class Arg0, class... DP, class... SP>
struct CommonSubview<Kokkos::Experimental::DynamicView<DP...>,
                     Kokkos::Experimental::DynamicView<SP...>, 1, Arg0> {
  using DstType          = Kokkos::Experimental::DynamicView<DP...>;
  using SrcType          = Kokkos::Experimental::DynamicView<SP...>;
  using dst_subview_type = DstType;
  using src_subview_type = SrcType;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& /*arg0*/)
      : dst_sub(dst), src_sub(src) {}
};

template <class... DP, class SrcType, class Arg0>
struct CommonSubview<Kokkos::Experimental::DynamicView<DP...>, SrcType, 1,
                     Arg0> {
  using DstType          = Kokkos::Experimental::DynamicView<DP...>;
  using dst_subview_type = DstType;
  using src_subview_type = typename Kokkos::Subview<SrcType, Arg0>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0)
      : dst_sub(dst), src_sub(src, arg0) {}
};

template <class DstType, class... SP, class Arg0>
struct CommonSubview<DstType, Kokkos::Experimental::DynamicView<SP...>, 1,
                     Arg0> {
  using SrcType          = Kokkos::Experimental::DynamicView<SP...>;
  using dst_subview_type = typename Kokkos::Subview<DstType, Arg0>;
  using src_subview_type = SrcType;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0)
      : dst_sub(dst, arg0), src_sub(src) {}
};

template <class... DP, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<Kokkos::Experimental::DynamicView<DP...>, ViewTypeB, Layout,
                ExecSpace, 1, iType> {
  Kokkos::Experimental::DynamicView<DP...> a;
  ViewTypeB b;

  using policy_type = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<iType>>;

  ViewCopy(const Kokkos::Experimental::DynamicView<DP...>& a_,
           const ViewTypeB& b_)
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-1D", policy_type(0, b.extent(0)),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0) const { a(i0) = b(i0); };
};

template <class... DP, class... SP, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<Kokkos::Experimental::DynamicView<DP...>,
                Kokkos::Experimental::DynamicView<SP...>, Layout, ExecSpace, 1,
                iType> {
  Kokkos::Experimental::DynamicView<DP...> a;
  Kokkos::Experimental::DynamicView<SP...> b;

  using policy_type = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<iType>>;

  ViewCopy(const Kokkos::Experimental::DynamicView<DP...>& a_,
           const Kokkos::Experimental::DynamicView<SP...>& b_)
      : a(a_), b(b_) {
    const iType n = std::min(a.extent(0), b.extent(0));
    Kokkos::parallel_for("Kokkos::ViewCopy-1D", policy_type(0, n), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0) const { a(i0) = b(i0); };
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_DYNAMIC_VIEW_HPP */
