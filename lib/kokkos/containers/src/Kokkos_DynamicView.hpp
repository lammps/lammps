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
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNAMICVIEW
#endif

#include <cstdio>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Experimental {

namespace Impl {

/// Utility class to manage memory for chunked arrays on the host and
/// device. Allocates/deallocates memory on both the host and device along with
/// providing utilities for creating mirrors and deep copying between them.
template <typename MemorySpace, typename ValueType>
struct ChunkedArrayManager {
  using value_type   = ValueType;
  using pointer_type = ValueType*;
  using track_type   = Kokkos::Impl::SharedAllocationTracker;

  ChunkedArrayManager()                           = default;
  ChunkedArrayManager(ChunkedArrayManager const&) = default;
  ChunkedArrayManager(ChunkedArrayManager&&)      = default;
  ChunkedArrayManager& operator=(ChunkedArrayManager&&) = default;
  ChunkedArrayManager& operator=(const ChunkedArrayManager&) = default;

  template <typename Space, typename Value>
  friend struct ChunkedArrayManager;

  template <typename Space, typename Value>
  inline ChunkedArrayManager(const ChunkedArrayManager<Space, Value>& rhs)
      : m_valid(rhs.m_valid),
        m_chunk_max(rhs.m_chunk_max),
        m_chunks((ValueType**)(rhs.m_chunks)),
        m_track(rhs.m_track),
        m_chunk_size(rhs.m_chunk_size) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<MemorySpace, Space>::assignable,
        "Incompatible ChunkedArrayManager copy construction");
  }

  ChunkedArrayManager(const unsigned arg_chunk_max,
                      const unsigned arg_chunk_size)
      : m_chunk_max(arg_chunk_max), m_chunk_size(arg_chunk_size) {}

 private:
  struct ACCESSIBLE_TAG {};
  struct INACCESSIBLE_TAG {};

  ChunkedArrayManager(ACCESSIBLE_TAG, pointer_type* arg_chunks,
                      const unsigned arg_chunk_max)
      : m_valid(true), m_chunk_max(arg_chunk_max), m_chunks(arg_chunks) {}

  ChunkedArrayManager(INACCESSIBLE_TAG, const unsigned arg_chunk_max,
                      const unsigned arg_chunk_size)
      : m_chunk_max(arg_chunk_max), m_chunk_size(arg_chunk_size) {}

 public:
  template <typename Space, typename Enable_ = void>
  struct IsAccessibleFrom;

  template <typename Space>
  struct IsAccessibleFrom<
      Space, typename std::enable_if_t<Kokkos::Impl::MemorySpaceAccess<
                 MemorySpace, Space>::accessible>> : std::true_type {};

  template <typename Space>
  struct IsAccessibleFrom<
      Space, typename std::enable_if_t<!Kokkos::Impl::MemorySpaceAccess<
                 MemorySpace, Space>::accessible>> : std::false_type {};

  template <typename Space>
  static ChunkedArrayManager<Space, ValueType> create_mirror(
      ChunkedArrayManager<MemorySpace, ValueType> const& other,
      std::enable_if_t<IsAccessibleFrom<Space>::value>* = nullptr) {
    return ChunkedArrayManager<Space, ValueType>{
        ACCESSIBLE_TAG{}, other.m_chunks, other.m_chunk_max};
  }

  template <typename Space>
  static ChunkedArrayManager<Space, ValueType> create_mirror(
      ChunkedArrayManager<MemorySpace, ValueType> const& other,
      std::enable_if_t<!IsAccessibleFrom<Space>::value>* = nullptr) {
    using tag_type =
        typename ChunkedArrayManager<Space, ValueType>::INACCESSIBLE_TAG;
    return ChunkedArrayManager<Space, ValueType>{tag_type{}, other.m_chunk_max,
                                                 other.m_chunk_size};
  }

 public:
  void allocate_device(const std::string& label) {
    if (m_chunks == nullptr) {
      m_chunks = reinterpret_cast<pointer_type*>(MemorySpace().allocate(
          label.c_str(), (sizeof(pointer_type) * (m_chunk_max + 2))));
    }
  }

  void initialize() {
    for (unsigned i = 0; i < m_chunk_max + 2; i++) {
      m_chunks[i] = nullptr;
    }
    m_valid = true;
  }

 private:
  /// Custom destroy functor for deallocating array chunks along with a linked
  /// allocation
  template <typename Space>
  struct Destroy {
    Destroy()               = default;
    Destroy(Destroy&&)      = default;
    Destroy(const Destroy&) = default;
    Destroy& operator=(Destroy&&) = default;
    Destroy& operator=(const Destroy&) = default;

    Destroy(std::string label, value_type** arg_chunk,
            const unsigned arg_chunk_max, const unsigned arg_chunk_size,
            value_type** arg_linked)
        : m_label(label),
          m_chunks(arg_chunk),
          m_linked(arg_linked),
          m_chunk_max(arg_chunk_max),
          m_chunk_size(arg_chunk_size) {}

    void execute() {
      // Destroy the array of chunk pointers.
      // Two entries beyond the max chunks are allocation counters.
      uintptr_t const len =
          *reinterpret_cast<uintptr_t*>(m_chunks + m_chunk_max);
      for (unsigned i = 0; i < len; i++) {
        Space().deallocate(m_label.c_str(), m_chunks[i],
                           sizeof(value_type) * m_chunk_size);
      }
      // Destroy the linked allocation if we have one.
      if (m_linked != nullptr) {
        Space().deallocate(m_label.c_str(), m_linked,
                           (sizeof(value_type*) * (m_chunk_max + 2)));
      }
    }

    void destroy_shared_allocation() { execute(); }

    std::string m_label;
    value_type** m_chunks = nullptr;
    value_type** m_linked = nullptr;
    unsigned m_chunk_max;
    unsigned m_chunk_size;
  };

 public:
  template <typename Space>
  void allocate_with_destroy(const std::string& label,
                             pointer_type* linked_allocation = nullptr) {
    using destroy_type = Destroy<Space>;
    using record_type =
        Kokkos::Impl::SharedAllocationRecord<MemorySpace, destroy_type>;

    // Allocate + 2 extra slots so that *m_chunk[m_chunk_max] ==
    // num_chunks_alloc and *m_chunk[m_chunk_max+1] == extent This must match in
    // Destroy's execute(...) method
    record_type* const record = record_type::allocate(
        MemorySpace(), label, (sizeof(pointer_type) * (m_chunk_max + 2)));
    m_chunks = static_cast<pointer_type*>(record->data());
    m_track.assign_allocated_record_to_uninitialized(record);

    record->m_destroy = destroy_type(label, m_chunks, m_chunk_max, m_chunk_size,
                                     linked_allocation);
  }

  pointer_type* get_ptr() const { return m_chunks; }

  template <typename OtherMemorySpace, typename ExecutionSpace>
  void deep_copy_to(
      const ExecutionSpace& exec_space,
      ChunkedArrayManager<OtherMemorySpace, ValueType> const& other) const {
    if (other.m_chunks != m_chunks) {
      Kokkos::Impl::DeepCopy<OtherMemorySpace, MemorySpace, ExecutionSpace>(
          exec_space, other.m_chunks, m_chunks,
          sizeof(pointer_type) * (m_chunk_max + 2));
    }
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type* operator+(int i) const { return m_chunks + i; }

  KOKKOS_INLINE_FUNCTION
  pointer_type& operator[](int i) const { return m_chunks[i]; }

  track_type const& track() const { return m_track; }

  KOKKOS_INLINE_FUNCTION
  bool valid() const { return m_valid; }

 private:
  bool m_valid           = false;
  unsigned m_chunk_max   = 0;
  pointer_type* m_chunks = nullptr;
  track_type m_track;
  unsigned m_chunk_size = 0;
};

} /* end namespace Impl */

/** \brief Dynamic views are restricted to rank-one and no layout.
 *         Resize only occurs on host outside of parallel_regions.
 *         Subviews are not allowed.
 */
template <typename DataType, typename... P>
class DynamicView : public Kokkos::ViewTraits<DataType, P...> {
 public:
  using traits = Kokkos::ViewTraits<DataType, P...>;

  using value_type   = typename traits::value_type;
  using device_space = typename traits::memory_space;
  using host_space =
      typename Kokkos::Impl::HostMirror<device_space>::Space::memory_space;
  using device_accessor = Impl::ChunkedArrayManager<device_space, value_type>;
  using host_accessor   = Impl::ChunkedArrayManager<host_space, value_type>;

 private:
  template <class, class...>
  friend class DynamicView;

  using track_type = Kokkos::Impl::SharedAllocationTracker;

  static_assert(traits::rank == 1 && traits::rank_dynamic == 1,
                "DynamicView must be rank-one");

  // It is assumed that the value_type is trivially copyable;
  // when this is not the case, potential problems can occur.
  static_assert(std::is_void<typename traits::specialize>::value,
                "DynamicView only implemented for non-specialized View type");

 private:
  device_accessor m_chunks;
  host_accessor m_chunks_host;
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
    uintptr_t n =
        *reinterpret_cast<const uintptr_t*>(m_chunks_host + m_chunk_max);
    return (n << m_chunk_shift);
  }

  KOKKOS_INLINE_FUNCTION
  size_t chunk_size() const noexcept { return m_chunk_size; }

  KOKKOS_INLINE_FUNCTION
  size_t chunk_max() const noexcept { return m_chunk_max; }

  KOKKOS_INLINE_FUNCTION
  size_t size() const noexcept {
    size_t extent_0 =
        *reinterpret_cast<const size_t*>(m_chunks_host + m_chunk_max + 1);
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
  int use_count() const { return m_chunks_host.track().use_count(); }

  inline const std::string label() const {
    return m_chunks_host.track().template get_label<host_space>();
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

    Kokkos::Impl::runtime_check_memory_access_violation<
        typename traits::memory_space>(
        "Kokkos::DynamicView ERROR: attempt to access inaccessible memory "
        "space");

    // Which chunk is being indexed.
    const uintptr_t ic = uintptr_t(i0) >> m_chunk_shift;

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    const uintptr_t n = *reinterpret_cast<uintptr_t*>(m_chunks + m_chunk_max);
    if (n <= ic) Kokkos::abort("Kokkos::DynamicView array bounds error");
#endif

    typename traits::value_type** const ch = m_chunks + ic;
    return (*ch)[i0 & m_chunk_mask];
  }

  //----------------------------------------
  /** \brief  Resizing in serial can grow or shrink the array size
   *          up to the maximum number of chunks
   * */
  template <typename IntType>
  inline void resize_serial(IntType const& n) {
    using local_value_type   = typename traits::value_type;
    using value_pointer_type = local_value_type*;

    const uintptr_t NC =
        (n + m_chunk_mask) >>
        m_chunk_shift;  // New total number of chunks needed for resize

    if (m_chunk_max < NC) {
      Kokkos::abort("DynamicView::resize_serial exceeded maximum size");
    }

    // *m_chunks[m_chunk_max] stores the current number of chunks being used
    uintptr_t* const pc =
        reinterpret_cast<uintptr_t*>(m_chunks_host + m_chunk_max);
    std::string _label = m_chunks_host.track().template get_label<host_space>();

    if (*pc < NC) {
      while (*pc < NC) {
        m_chunks_host[*pc] =
            reinterpret_cast<value_pointer_type>(device_space().allocate(
                _label.c_str(), sizeof(local_value_type) << m_chunk_shift));
        ++*pc;
      }
    } else {
      while (NC + 1 <= *pc) {
        --*pc;
        device_space().deallocate(_label.c_str(), m_chunks_host[*pc],
                                  sizeof(local_value_type) << m_chunk_shift);
        m_chunks_host[*pc] = nullptr;
      }
    }
    // *m_chunks_host[m_chunk_max+1] stores the 'extent' requested by resize
    *(pc + 1) = n;

    typename device_space::execution_space exec{};
    m_chunks_host.deep_copy_to(exec, m_chunks);
    exec.fence(
        "DynamicView::resize_serial: Fence after copying chunks to the device");
  }

  KOKKOS_INLINE_FUNCTION bool is_allocated() const {
    if (m_chunks_host.valid()) {
      // *m_chunks_host[m_chunk_max] stores the current number of chunks being
      // used
      uintptr_t* const pc =
          reinterpret_cast<uintptr_t*>(m_chunks_host + m_chunk_max);
      return (*(pc + 1) > 0);
    } else {
      return false;
    }
  }

  KOKKOS_FUNCTION const device_accessor& impl_get_chunks() const {
    return m_chunks;
  }

  KOKKOS_FUNCTION device_accessor& impl_get_chunks() { return m_chunks; }

  //----------------------------------------------------------------------

  ~DynamicView()                  = default;
  DynamicView()                   = default;
  DynamicView(DynamicView&&)      = default;
  DynamicView(const DynamicView&) = default;
  DynamicView& operator=(DynamicView&&) = default;
  DynamicView& operator=(const DynamicView&) = default;

  template <class RT, class... RP>
  DynamicView(const DynamicView<RT, RP...>& rhs)
      : m_chunks(rhs.m_chunks),
        m_chunks_host(rhs.m_chunks_host),
        m_chunk_shift(rhs.m_chunk_shift),
        m_chunk_mask(rhs.m_chunk_mask),
        m_chunk_max(rhs.m_chunk_max),
        m_chunk_size(rhs.m_chunk_size) {
    using SrcTraits = typename DynamicView<RT, RP...>::traits;
    using Mapping   = Kokkos::Impl::ViewMapping<traits, SrcTraits, void>;
    static_assert(Mapping::is_assignable,
                  "Incompatible DynamicView copy construction");
  }

  /**\brief  Allocation constructor
   *
   *  Memory is allocated in chunks
   *  A maximum size is required in order to allocate a
   *  chunk-pointer array.
   */
  template <class... Prop>
  DynamicView(const Kokkos::Impl::ViewCtorProp<Prop...>& arg_prop,
              const unsigned min_chunk_size,
              const unsigned max_extent)
      :  // The chunk size is guaranteed to be a power of two
        m_chunk_shift(Kokkos::Impl::integral_power_of_two_that_contains(
            min_chunk_size))  // div ceil(log2(min_chunk_size))
        ,
        m_chunk_mask((1 << m_chunk_shift) - 1)  // mod
        ,
        m_chunk_max((max_extent + m_chunk_mask) >>
                    m_chunk_shift)  // max num pointers-to-chunks in array
        ,
        m_chunk_size(2 << (m_chunk_shift - 1)) {
    m_chunks = device_accessor(m_chunk_max, m_chunk_size);

    const std::string& label =
        static_cast<Kokkos::Impl::ViewCtorProp<void, std::string> const&>(
            arg_prop)
            .value;

    if (device_accessor::template IsAccessibleFrom<host_space>::value) {
      m_chunks.template allocate_with_destroy<device_space>(label);
      m_chunks.initialize();
      m_chunks_host =
          device_accessor::template create_mirror<host_space>(m_chunks);
    } else {
      m_chunks.allocate_device(label);
      m_chunks_host =
          device_accessor::template create_mirror<host_space>(m_chunks);
      m_chunks_host.template allocate_with_destroy<device_space>(
          label, m_chunks.get_ptr());
      m_chunks_host.initialize();

      // Add some properties if not provided to avoid need for if constexpr
      using alloc_prop_input = Kokkos::Impl::ViewCtorProp<Prop...>;
      using alloc_prop       = Kokkos::Impl::ViewCtorProp<
          Prop..., std::conditional_t<alloc_prop_input::has_execution_space,
                                      std::integral_constant<unsigned int, 15>,
                                      typename device_space::execution_space>>;
      alloc_prop arg_prop_copy(arg_prop);

      const auto& exec = static_cast<const Kokkos::Impl::ViewCtorProp<
          void, typename alloc_prop::execution_space>&>(arg_prop_copy)
                             .value;
      m_chunks_host.deep_copy_to(exec, m_chunks);
      if (!alloc_prop_input::has_execution_space)
        exec.fence(
            "DynamicView::DynamicView(): Fence after copying chunks to the "
            "device");
    }
  }

  DynamicView(const std::string& arg_label, const unsigned min_chunk_size,
              const unsigned max_extent)
      : DynamicView(Kokkos::view_alloc(arg_label), min_chunk_size, max_extent) {
  }
};

}  // namespace Experimental

template <class>
struct is_dynamic_view : public std::false_type {};

template <class D, class... P>
struct is_dynamic_view<Kokkos::Experimental::DynamicView<D, P...>>
    : public std::true_type {};

}  // namespace Kokkos

namespace Kokkos {

namespace Impl {

// Deduce Mirror Types
template <class Space, class T, class... P>
struct MirrorDynamicViewType {
  // The incoming view_type
  using src_view_type = typename Kokkos::Experimental::DynamicView<T, P...>;
  // The memory space for the mirror view
  using memory_space = typename Space::memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  using array_layout = typename src_view_type::array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.)
  using data_type = typename src_view_type::non_const_data_type;
  // The destination view type if it is not the same memory space
  using dest_view_type =
      Kokkos::Experimental::DynamicView<data_type, array_layout, Space>;
  // If it is the same memory_space return the existing view_type
  // This will also keep the unmanaged trait if necessary
  using view_type =
      std::conditional_t<is_same_memspace, src_view_type, dest_view_type>;
};
}  // namespace Impl

namespace Impl {
template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(
    const Kokkos::Experimental::DynamicView<T, P...>& src,
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    std::enable_if_t<!Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space>* =
        nullptr) {
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

  static_assert(
      !alloc_prop_input::has_label,
      "The view constructor arguments passed to Kokkos::create_mirror "
      "must not include a label!");
  static_assert(
      !alloc_prop_input::has_pointer,
      "The view constructor arguments passed to Kokkos::create_mirror must "
      "not include a pointer!");
  static_assert(
      !alloc_prop_input::allow_padding,
      "The view constructor arguments passed to Kokkos::create_mirror must "
      "not explicitly allow padding!");

  using alloc_prop = Impl::ViewCtorProp<ViewCtorArgs..., std::string>;
  alloc_prop prop_copy(arg_prop);
  static_cast<Impl::ViewCtorProp<void, std::string>&>(prop_copy).value =
      std::string(src.label()).append("_mirror");

  auto ret = typename Kokkos::Experimental::DynamicView<T, P...>::HostMirror(
      prop_copy, src.chunk_size(), src.chunk_max() * src.chunk_size());

  ret.resize_serial(src.extent(0));

  return ret;
}

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(
    const Kokkos::Experimental::DynamicView<T, P...>& src,
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    std::enable_if_t<Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space>* =
        nullptr) {
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

  static_assert(
      !alloc_prop_input::has_label,
      "The view constructor arguments passed to Kokkos::create_mirror "
      "must not include a label!");
  static_assert(
      !alloc_prop_input::has_pointer,
      "The view constructor arguments passed to Kokkos::create_mirror must "
      "not include a pointer!");
  static_assert(
      !alloc_prop_input::allow_padding,
      "The view constructor arguments passed to Kokkos::create_mirror must "
      "not explicitly allow padding!");

  using MemorySpace = typename alloc_prop_input::memory_space;
  using alloc_prop  = Impl::ViewCtorProp<ViewCtorArgs..., std::string>;
  alloc_prop prop_copy(arg_prop);
  static_cast<Impl::ViewCtorProp<void, std::string>&>(prop_copy).value =
      std::string(src.label()).append("_mirror");

  auto ret = typename Kokkos::Impl::MirrorDynamicViewType<
      MemorySpace, T, P...>::view_type(prop_copy, src.chunk_size(),
                                       src.chunk_max() * src.chunk_size());

  ret.resize_serial(src.extent(0));

  return ret;
}
}  // namespace Impl

// Create a mirror in host space
template <class T, class... P>
inline auto create_mirror(
    const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror(src, Impl::ViewCtorProp<>{});
}

template <class T, class... P>
inline auto create_mirror(
    Kokkos::Impl::WithoutInitializing_t wi,
    const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror(src, Kokkos::view_alloc(wi));
}

// Create a mirror in a new space
template <class Space, class T, class... P>
inline auto create_mirror(
    const Space&, const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror(
      src, Kokkos::view_alloc(typename Space::memory_space{}));
}

template <class Space, class T, class... P>
typename Kokkos::Impl::MirrorDynamicViewType<Space, T, P...>::view_type
create_mirror(Kokkos::Impl::WithoutInitializing_t wi, const Space&,
              const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror(
      src, Kokkos::view_alloc(wi, typename Space::memory_space{}));
}

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror(src, arg_prop);
}

namespace Impl {

template <class T, class... P, class... ViewCtorArgs>
inline std::enable_if_t<
    !Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space &&
        (std::is_same<
             typename Kokkos::Experimental::DynamicView<T, P...>::memory_space,
             typename Kokkos::Experimental::DynamicView<
                 T, P...>::HostMirror::memory_space>::value &&
         std::is_same<
             typename Kokkos::Experimental::DynamicView<T, P...>::data_type,
             typename Kokkos::Experimental::DynamicView<
                 T, P...>::HostMirror::data_type>::value),
    typename Kokkos::Experimental::DynamicView<T, P...>::HostMirror>
create_mirror_view(const Kokkos::Experimental::DynamicView<T, P...>& src,
                   const Impl::ViewCtorProp<ViewCtorArgs...>&) {
  return src;
}

template <class T, class... P, class... ViewCtorArgs>
inline std::enable_if_t<
    !Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space &&
        !(std::is_same<
              typename Kokkos::Experimental::DynamicView<T, P...>::memory_space,
              typename Kokkos::Experimental::DynamicView<
                  T, P...>::HostMirror::memory_space>::value &&
          std::is_same<
              typename Kokkos::Experimental::DynamicView<T, P...>::data_type,
              typename Kokkos::Experimental::DynamicView<
                  T, P...>::HostMirror::data_type>::value),
    typename Kokkos::Experimental::DynamicView<T, P...>::HostMirror>
create_mirror_view(const Kokkos::Experimental::DynamicView<T, P...>& src,
                   const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop) {
  return Kokkos::create_mirror(arg_prop, src);
}

template <class T, class... P, class... ViewCtorArgs,
          class = std::enable_if_t<
              Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space>>
std::enable_if_t<Impl::MirrorDynamicViewType<
                     typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space,
                     T, P...>::is_same_memspace,
                 typename Impl::MirrorDynamicViewType<
                     typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space,
                     T, P...>::view_type>
create_mirror_view(const Kokkos::Experimental::DynamicView<T, P...>& src,
                   const Impl::ViewCtorProp<ViewCtorArgs...>&) {
  return src;
}

template <class T, class... P, class... ViewCtorArgs,
          class = std::enable_if_t<
              Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space>>
std::enable_if_t<!Impl::MirrorDynamicViewType<
                     typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space,
                     T, P...>::is_same_memspace,
                 typename Impl::MirrorDynamicViewType<
                     typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space,
                     T, P...>::view_type>
create_mirror_view(const Kokkos::Experimental::DynamicView<T, P...>& src,
                   const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop) {
  return Kokkos::Impl::create_mirror(src, arg_prop);
}
}  // namespace Impl

// Create a mirror view in host space
template <class T, class... P>
inline auto create_mirror_view(
    const typename Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror_view(src, Impl::ViewCtorProp<>{});
}

template <class T, class... P>
inline auto create_mirror_view(
    Kokkos::Impl::WithoutInitializing_t wi,
    const typename Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror_view(src, Kokkos::view_alloc(wi));
}

// Create a mirror in a new space
template <class Space, class T, class... P>
inline auto create_mirror_view(
    const Space&, const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror_view(src,
                                  view_alloc(typename Space::memory_space{}));
}

template <class Space, class T, class... P>
inline auto create_mirror_view(
    Kokkos::Impl::WithoutInitializing_t wi, const Space&,
    const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror_view(
      src, Kokkos::view_alloc(wi, typename Space::memory_space{}));
}

template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror_view(
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    const Kokkos::Experimental::DynamicView<T, P...>& src) {
  return Impl::create_mirror_view(src, arg_prop);
}

template <class T, class... DP, class... SP>
inline void deep_copy(const Kokkos::Experimental::DynamicView<T, DP...>& dst,
                      const Kokkos::Experimental::DynamicView<T, SP...>& src) {
  using dst_type = Kokkos::Experimental::DynamicView<T, DP...>;
  using src_type = Kokkos::Experimental::DynamicView<T, SP...>;

  using dst_execution_space = typename ViewTraits<T, DP...>::execution_space;
  using src_execution_space = typename ViewTraits<T, SP...>::execution_space;
  using dst_memory_space    = typename ViewTraits<T, DP...>::memory_space;
  using src_memory_space    = typename ViewTraits<T, SP...>::memory_space;

  constexpr bool DstExecCanAccessSrc =
      Kokkos::SpaceAccessibility<dst_execution_space,
                                 src_memory_space>::accessible;
  constexpr bool SrcExecCanAccessDst =
      Kokkos::SpaceAccessibility<src_execution_space,
                                 dst_memory_space>::accessible;

  if (DstExecCanAccessSrc)
    Kokkos::Impl::ViewRemap<dst_type, src_type, dst_execution_space>(dst, src);
  else if (SrcExecCanAccessDst)
    Kokkos::Impl::ViewRemap<dst_type, src_type, src_execution_space>(dst, src);
  else
    src.impl_get_chunks().deep_copy_to(dst_execution_space{},
                                       dst.impl_get_chunks());
  Kokkos::fence("Kokkos::deep_copy(DynamicView)");
}

template <class ExecutionSpace, class T, class... DP, class... SP>
inline void deep_copy(const ExecutionSpace& exec,
                      const Kokkos::Experimental::DynamicView<T, DP...>& dst,
                      const Kokkos::Experimental::DynamicView<T, SP...>& src) {
  using dst_type = Kokkos::Experimental::DynamicView<T, DP...>;
  using src_type = Kokkos::Experimental::DynamicView<T, SP...>;

  using dst_execution_space = typename ViewTraits<T, DP...>::execution_space;
  using src_execution_space = typename ViewTraits<T, SP...>::execution_space;
  using dst_memory_space    = typename ViewTraits<T, DP...>::memory_space;
  using src_memory_space    = typename ViewTraits<T, SP...>::memory_space;

  constexpr bool DstExecCanAccessSrc =
      Kokkos::SpaceAccessibility<dst_execution_space,
                                 src_memory_space>::accessible;
  constexpr bool SrcExecCanAccessDst =
      Kokkos::SpaceAccessibility<src_execution_space,
                                 dst_memory_space>::accessible;

  // FIXME use execution space
  if (DstExecCanAccessSrc)
    Kokkos::Impl::ViewRemap<dst_type, src_type, dst_execution_space>(dst, src);
  else if (SrcExecCanAccessDst)
    Kokkos::Impl::ViewRemap<dst_type, src_type, src_execution_space>(dst, src);
  else
    src.impl_get_chunks().deep_copy_to(exec, dst.impl_get_chunks());
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
        Kokkos::SpaceAccessibility<dst_execution_space,
                                   src_memory_space>::accessible
  };

  if (DstExecCanAccessSrc) {
    // Copying data between views in accessible memory spaces and either
    // non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap<dst_type, src_type>(dst, src);
    Kokkos::fence("Kokkos::deep_copy(DynamicView)");
  } else {
    Kokkos::Impl::throw_runtime_exception(
        "deep_copy given views that would require a temporary allocation");
  }
}

template <class T, class... DP, class... SP>
inline void deep_copy(const Kokkos::Experimental::DynamicView<T, DP...>& dst,
                      const View<T, SP...>& src) {
  using dst_type = Kokkos::Experimental::DynamicView<T, DP...>;
  using src_type = View<T, SP...>;

  using dst_execution_space = typename ViewTraits<T, DP...>::execution_space;
  using src_memory_space    = typename ViewTraits<T, SP...>::memory_space;

  enum {
    DstExecCanAccessSrc =
        Kokkos::SpaceAccessibility<dst_execution_space,
                                   src_memory_space>::accessible
  };

  if (DstExecCanAccessSrc) {
    // Copying data between views in accessible memory spaces and either
    // non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap<dst_type, src_type>(dst, src);
    Kokkos::fence("Kokkos::deep_copy(DynamicView)");
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

template <class... ViewCtorArgs, class T, class... P>
auto create_mirror_view_and_copy(
    const Impl::ViewCtorProp<ViewCtorArgs...>&,
    const Kokkos::Experimental::DynamicView<T, P...>& src,
    std::enable_if_t<
        std::is_void<typename ViewTraits<T, P...>::specialize>::value &&
        Impl::MirrorDynamicViewType<
            typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space, T,
            P...>::is_same_memspace>* = nullptr) {
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;
  static_assert(
      alloc_prop_input::has_memory_space,
      "The view constructor arguments passed to "
      "Kokkos::create_mirror_view_and_copy must include a memory space!");
  static_assert(!alloc_prop_input::has_pointer,
                "The view constructor arguments passed to "
                "Kokkos::create_mirror_view_and_copy must "
                "not include a pointer!");
  static_assert(!alloc_prop_input::allow_padding,
                "The view constructor arguments passed to "
                "Kokkos::create_mirror_view_and_copy must "
                "not explicitly allow padding!");

  // same behavior as deep_copy(src, src)
  if (!alloc_prop_input::has_execution_space)
    fence(
        "Kokkos::create_mirror_view_and_copy: fence before returning src view");
  return src;
}

template <class... ViewCtorArgs, class T, class... P>
auto create_mirror_view_and_copy(
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    const Kokkos::Experimental::DynamicView<T, P...>& src,
    std::enable_if_t<
        std::is_void<typename ViewTraits<T, P...>::specialize>::value &&
        !Impl::MirrorDynamicViewType<
            typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space, T,
            P...>::is_same_memspace>* = nullptr) {
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;
  static_assert(
      alloc_prop_input::has_memory_space,
      "The view constructor arguments passed to "
      "Kokkos::create_mirror_view_and_copy must include a memory space!");
  static_assert(!alloc_prop_input::has_pointer,
                "The view constructor arguments passed to "
                "Kokkos::create_mirror_view_and_copy must "
                "not include a pointer!");
  static_assert(!alloc_prop_input::allow_padding,
                "The view constructor arguments passed to "
                "Kokkos::create_mirror_view_and_copy must "
                "not explicitly allow padding!");
  using Space = typename alloc_prop_input::memory_space;
  using Mirror =
      typename Impl::MirrorDynamicViewType<Space, T, P...>::view_type;

  // Add some properties if not provided to avoid need for if constexpr
  using alloc_prop = Impl::ViewCtorProp<
      ViewCtorArgs...,
      std::conditional_t<alloc_prop_input::has_label,
                         std::integral_constant<unsigned int, 12>, std::string>,
      std::conditional_t<!alloc_prop_input::initialize,
                         std::integral_constant<unsigned int, 13>,
                         Impl::WithoutInitializing_t>,
      std::conditional_t<alloc_prop_input::has_execution_space,
                         std::integral_constant<unsigned int, 14>,
                         typename Space::execution_space>>;
  alloc_prop arg_prop_copy(arg_prop);

  std::string& label =
      static_cast<Impl::ViewCtorProp<void, std::string>&>(arg_prop_copy).value;
  if (label.empty()) label = src.label();
  auto mirror = typename Mirror::non_const_type(
      arg_prop_copy, src.chunk_size(), src.chunk_max() * src.chunk_size());
  mirror.resize_serial(src.extent(0));
  if (alloc_prop_input::has_execution_space) {
    using ExecutionSpace = typename alloc_prop::execution_space;
    deep_copy(
        static_cast<Impl::ViewCtorProp<void, ExecutionSpace>&>(arg_prop_copy)
            .value,
        mirror, src);
  } else
    deep_copy(mirror, src);
  return mirror;
}

template <class Space, class T, class... P>
auto create_mirror_view_and_copy(
    const Space&, const Kokkos::Experimental::DynamicView<T, P...>& src,
    std::string const& name = "") {
  return create_mirror_view_and_copy(
      Kokkos::view_alloc(typename Space::memory_space{}, name), src);
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNAMICVIEW
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNAMICVIEW
#endif
#endif /* #ifndef KOKKOS_DYNAMIC_VIEW_HPP */
