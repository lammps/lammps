//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

/// \file Kokkos_UnorderedMap.hpp
/// \brief Declaration and definition of Kokkos::UnorderedMap.
///
/// This header file declares and defines Kokkos::UnorderedMap and its
/// related nonmember functions.

#ifndef KOKKOS_UNORDERED_MAP_HPP
#define KOKKOS_UNORDERED_MAP_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_UNORDEREDMAP
#endif

#include <Kokkos_Core.hpp>
#include <Kokkos_Functional.hpp>

#include <Kokkos_Bitset.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_UnorderedMap_impl.hpp>
#include <impl/Kokkos_ViewCtor.hpp>

#include <cstdint>

namespace Kokkos {

enum : unsigned { UnorderedMapInvalidIndex = ~0u };

/// \brief First element of the return value of UnorderedMap::insert().
///
/// Inserting an element into an UnorderedMap is not guaranteed to
/// succeed.  There are three possible conditions:
/// <ol>
/// <li> <tt>INSERT_FAILED</tt>: The insert failed.  This usually
///      means that the UnorderedMap ran out of space. </li>
/// <li> <tt>INSERT_SUCCESS</tt>: The insert succeeded, and the key
///      did <i>not</i> exist in the table before. </li>
/// <li> <tt>INSERT_EXISTING</tt>: The insert succeeded, and the key
///      <i>did</i> exist in the table before.  The new value was
///      ignored and the old value was left in place. </li>
/// </ol>

class UnorderedMapInsertResult {
 private:
  enum Status : uint32_t {
    SUCCESS          = 1u << 31,
    EXISTING         = 1u << 30,
    FREED_EXISTING   = 1u << 29,
    LIST_LENGTH_MASK = ~(SUCCESS | EXISTING | FREED_EXISTING)
  };

 public:
  /// Did the map successful insert the key/value pair
  KOKKOS_FORCEINLINE_FUNCTION
  bool success() const { return (m_status & SUCCESS); }

  /// Was the key already present in the map
  KOKKOS_FORCEINLINE_FUNCTION
  bool existing() const { return (m_status & EXISTING); }

  /// Did the map fail to insert the key due to insufficient capacity
  KOKKOS_FORCEINLINE_FUNCTION
  bool failed() const { return m_index == UnorderedMapInvalidIndex; }

  /// Did the map lose a race condition to insert a dupulicate key/value pair
  /// where an index was claimed that needed to be released
  KOKKOS_FORCEINLINE_FUNCTION
  bool freed_existing() const { return (m_status & FREED_EXISTING); }

  /// How many iterations through the insert loop did it take before the
  /// map returned
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t list_position() const { return (m_status & LIST_LENGTH_MASK); }

  /// Index where the key can be found as long as the insert did not fail
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t index() const { return m_index; }

  KOKKOS_FORCEINLINE_FUNCTION
  UnorderedMapInsertResult() : m_index(UnorderedMapInvalidIndex), m_status(0) {}

  KOKKOS_FORCEINLINE_FUNCTION
  void increment_list_position() {
    m_status += (list_position() < LIST_LENGTH_MASK) ? 1u : 0u;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_existing(uint32_t i, bool arg_freed_existing) {
    m_index = i;
    m_status =
        EXISTING | (arg_freed_existing ? FREED_EXISTING : 0u) | list_position();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_success(uint32_t i) {
    m_index  = i;
    m_status = SUCCESS | list_position();
  }

 private:
  uint32_t m_index;
  uint32_t m_status;
};

/// \class UnorderedMapInsertOpTypes
///
/// \brief Operations applied to the values array upon subsequent insertions.
///
/// The default behavior when a k,v pair already exists in the UnorderedMap is
/// to perform no operation. Alternatively, the caller may select to
/// instantiate the UnorderedMap with the AtomicAdd insert operator such that
/// duplicate keys accumulate values into the given values array entry.
/// \tparam ValueTypeView The UnorderedMap value array type.
/// \tparam ValuesIdxType The index type for lookups in the value array.
///
/// Supported operations:
///   NoOp:      the first key inserted stores the associated value.
///   AtomicAdd: duplicate key insertions sum values together.
template <class ValueTypeView, class ValuesIdxType>
struct UnorderedMapInsertOpTypes {
  using value_type = typename ValueTypeView::non_const_value_type;
  struct NoOp {
    KOKKOS_FUNCTION
    void op(ValueTypeView, ValuesIdxType, const value_type) const {}
  };
  struct AtomicAdd {
    KOKKOS_FUNCTION
    void op(ValueTypeView values, ValuesIdxType values_idx,
            const value_type v) const {
      Kokkos::atomic_add(values.data() + values_idx, v);
    }
  };
};

/// \class UnorderedMap
/// \brief Thread-safe, performance-portable lookup table.
///
/// This class provides a lookup table.  In terms of functionality,
/// this class compares to std::unordered_map (new in C++11).
/// "Unordered" means that keys are not stored in any particular
/// order, unlike (for example) std::map.  "Thread-safe" means that
/// lookups, insertion, and deletion are safe to call by multiple
/// threads in parallel.  "Performance-portable" means that parallel
/// performance of these operations is reasonable, on multiple
/// hardware platforms.  Platforms on which performance has been
/// tested include conventional Intel x86 multicore processors, Intel
/// Xeon Phi ("MIC"), and NVIDIA GPUs.
///
/// Parallel performance portability entails design decisions that
/// might differ from one's expectation for a sequential interface.
/// This particularly affects insertion of single elements.  In an
/// interface intended for sequential use, insertion might reallocate
/// memory if the original allocation did not suffice to hold the new
/// element.  In this class, insertion does <i>not</i> reallocate
/// memory.  This means that it might fail.  insert() returns an enum
/// which indicates whether the insert failed.  There are three
/// possible conditions:
/// <ol>
/// <li> <tt>INSERT_FAILED</tt>: The insert failed.  This usually
///      means that the UnorderedMap ran out of space. </li>
/// <li> <tt>INSERT_SUCCESS</tt>: The insert succeeded, and the key
///      did <i>not</i> exist in the table before. </li>
/// <li> <tt>INSERT_EXISTING</tt>: The insert succeeded, and the key
///      <i>did</i> exist in the table before.  The new value was
///      ignored and the old value was left in place. </li>
/// </ol>
///
/// \tparam Key Type of keys of the lookup table.  If \c const, users
///   are not allowed to add or remove keys, though they are allowed
///   to change values.  In that case, the implementation may make
///   optimizations specific to the <tt>Device</tt>.  For example, if
///   <tt>Device</tt> is \c Cuda, it may use texture fetches to access
///   keys.
///
/// \tparam Value Type of values stored in the lookup table.  You may use
///   \c void here, in which case the table will be a set of keys.  If
///   \c const, users are not allowed to change entries.
///   In that case, the implementation may make
///   optimizations specific to the \c Device, such as using texture
///   fetches to access values.
///
/// \tparam Device The Kokkos Device type.
///
/// \tparam Hasher Definition of the hash function for instances of
///   <tt>Key</tt>.  The default will calculate a bitwise hash.
///
/// \tparam EqualTo Definition of the equality function for instances of
///   <tt>Key</tt>.  The default will do a bitwise equality comparison.
///
template <typename Key, typename Value,
          typename Device  = Kokkos::DefaultExecutionSpace,
          typename Hasher  = pod_hash<std::remove_const_t<Key>>,
          typename EqualTo = pod_equal_to<std::remove_const_t<Key>>>
class UnorderedMap {
 private:
  using host_mirror_space =
      typename ViewTraits<Key, Device, void, void>::host_mirror_space;

 public:
  //! \name Public types and constants
  //@{
  // key_types
  using declared_key_type = Key;
  using key_type          = std::remove_const_t<declared_key_type>;
  using const_key_type    = std::add_const_t<key_type>;

  // value_types
  using declared_value_type = Value;
  using value_type          = std::remove_const_t<declared_value_type>;
  using const_value_type    = std::add_const_t<value_type>;

  using device_type     = Device;
  using execution_space = typename Device::execution_space;
  using hasher_type     = Hasher;
  using equal_to_type   = EqualTo;
  using size_type       = uint32_t;

  // map_types
  using declared_map_type =
      UnorderedMap<declared_key_type, declared_value_type, device_type,
                   hasher_type, equal_to_type>;
  using insertable_map_type = UnorderedMap<key_type, value_type, device_type,
                                           hasher_type, equal_to_type>;
  using modifiable_map_type =
      UnorderedMap<const_key_type, value_type, device_type, hasher_type,
                   equal_to_type>;
  using const_map_type = UnorderedMap<const_key_type, const_value_type,
                                      device_type, hasher_type, equal_to_type>;

  static constexpr bool is_set = std::is_void_v<value_type>;
  static constexpr bool has_const_key =
      std::is_same_v<const_key_type, declared_key_type>;
  static constexpr bool has_const_value =
      is_set || std::is_same_v<const_value_type, declared_value_type>;

  static constexpr bool is_insertable_map =
      !has_const_key && (is_set || !has_const_value);
  static constexpr bool is_modifiable_map = has_const_key && !has_const_value;
  static constexpr bool is_const_map      = has_const_key && has_const_value;

  using insert_result = UnorderedMapInsertResult;

  using HostMirror =
      UnorderedMap<Key, Value, host_mirror_space, Hasher, EqualTo>;

  using histogram_type = Impl::UnorderedMapHistogram<const_map_type>;
  //@}

 private:
  enum : size_type { invalid_index = ~static_cast<size_type>(0) };

  using impl_value_type = std::conditional_t<is_set, int, declared_value_type>;

  using key_type_view = std::conditional_t<
      is_insertable_map, View<key_type *, device_type>,
      View<const key_type *, device_type, MemoryTraits<RandomAccess>>>;

  using value_type_view = std::conditional_t<
      is_insertable_map || is_modifiable_map,
      View<impl_value_type *, device_type>,
      View<const impl_value_type *, device_type, MemoryTraits<RandomAccess>>>;

  using size_type_view = std::conditional_t<
      is_insertable_map, View<size_type *, device_type>,
      View<const size_type *, device_type, MemoryTraits<RandomAccess>>>;

  using bitset_type = std::conditional_t<is_insertable_map, Bitset<Device>,
                                         ConstBitset<Device>>;

  enum { modified_idx = 0, erasable_idx = 1, failed_insert_idx = 2 };
  enum { num_scalars = 3 };
  using scalars_view = View<int[num_scalars], LayoutLeft, device_type>;

 public:
  //! \name Public member functions
  //@{
  using default_op_type =
      typename UnorderedMapInsertOpTypes<value_type_view, uint32_t>::NoOp;

  /// \brief Constructor
  ///
  /// \param capacity_hint [in] Initial guess of how many unique keys will be
  ///                           inserted into the map.
  /// \param hash          [in] Hasher function for \c Key instances.  The
  ///                           default value usually suffices.
  /// \param equal_to      [in] The operator used for determining if two
  ///                           keys are equal.
  UnorderedMap(size_type capacity_hint = 0, hasher_type hasher = hasher_type(),
               equal_to_type equal_to = equal_to_type())
      : UnorderedMap(Kokkos::view_alloc(), capacity_hint, hasher, equal_to) {}

  template <class... P>
  UnorderedMap(const Impl::ViewCtorProp<P...> &arg_prop,
               size_type capacity_hint = 0, hasher_type hasher = hasher_type(),
               equal_to_type equal_to = equal_to_type())
      : m_bounded_insert(true), m_hasher(hasher), m_equal_to(equal_to) {
    if (!is_insertable_map) {
      Kokkos::Impl::throw_runtime_exception(
          "Cannot construct a non-insertable (i.e. const key_type) "
          "unordered_map");
    }

    //! Ensure that allocation properties are consistent.
    using alloc_prop_t = std::decay_t<decltype(arg_prop)>;
    static_assert(alloc_prop_t::initialize,
                  "Allocation property 'initialize' should be true.");
    static_assert(
        !alloc_prop_t::has_pointer,
        "Allocation properties should not contain the 'pointer' property.");

    /// Update allocation properties with 'label' and 'without initializing'
    /// properties.
    const auto prop_copy =
        Impl::with_properties_if_unset(arg_prop, std::string("UnorderedMap"));
    const auto prop_copy_noinit =
        Impl::with_properties_if_unset(prop_copy, Kokkos::WithoutInitializing);

    //! Initialize member views.
    m_size = shared_size_t(Kokkos::view_alloc(
        Kokkos::DefaultHostExecutionSpace{},
        Impl::get_property<Impl::LabelTag>(prop_copy) + " - size"));

    m_available_indexes =
        bitset_type(Kokkos::Impl::append_to_label(prop_copy, " - bitset"),
                    calculate_capacity(capacity_hint));

    m_hash_lists = size_type_view(
        Kokkos::Impl::append_to_label(prop_copy_noinit, " - hash list"),
        Impl::find_hash_size(capacity()));

    m_next_index = size_type_view(
        Kokkos::Impl::append_to_label(prop_copy_noinit, " - next index"),
        capacity() + 1);  // +1 so that the *_at functions can always return a
                          // valid reference

    m_keys = key_type_view(Kokkos::Impl::append_to_label(prop_copy, " - keys"),
                           capacity());

    m_values =
        value_type_view(Kokkos::Impl::append_to_label(prop_copy, " - values"),
                        is_set ? 0 : capacity());

    m_scalars =
        scalars_view(Kokkos::Impl::append_to_label(prop_copy, " - scalars"));

    /**
     * Deep copies should also be done using the space instance if given.
     * Instead of the if/else we could use the
     * @c get_property_or_default, but giving even the default execution space
     * instance will change the behavior of @c deep_copy.
     */
    if constexpr (alloc_prop_t::has_execution_space) {
      const auto &space = Impl::get_property<Impl::ExecutionSpaceTag>(arg_prop);
      Kokkos::deep_copy(space, m_hash_lists, invalid_index);
      Kokkos::deep_copy(space, m_next_index, invalid_index);
    } else {
      Kokkos::deep_copy(m_hash_lists, invalid_index);
      Kokkos::deep_copy(m_next_index, invalid_index);
    }
  }

  void reset_failed_insert_flag() { reset_flag(failed_insert_idx); }

  histogram_type get_histogram() { return histogram_type(*this); }

  //! Clear all entries in the table.
  void clear() {
    m_bounded_insert = true;

    if (capacity() == 0) return;

    m_available_indexes.clear();

    Kokkos::deep_copy(m_hash_lists, invalid_index);
    Kokkos::deep_copy(m_next_index, invalid_index);
    {
      const key_type tmp = key_type();
      Kokkos::deep_copy(m_keys, tmp);
    }
    Kokkos::deep_copy(m_scalars, 0);
    m_size() = 0;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool is_allocated() const {
    return (m_keys.is_allocated() && (is_set || m_values.is_allocated()) &&
            m_scalars.is_allocated());
  }

  /// \brief Change the capacity of the the map
  ///
  /// If there are no failed inserts the current size of the map will
  /// be used as a lower bound for the input capacity.
  /// If the map is not empty and does not have failed inserts
  /// and the capacity changes then the current data is copied
  /// into the resized / rehashed map.
  ///
  /// This is <i>not</i> a device function; it may <i>not</i> be
  /// called in a parallel kernel.
  bool rehash(size_type requested_capacity = 0) {
    const bool bounded_insert = (capacity() == 0) || (size() == 0u);
    return rehash(requested_capacity, bounded_insert);
  }

  bool rehash(size_type requested_capacity, bool bounded_insert) {
    if (!is_insertable_map) return false;

    const size_type curr_size = size();
    requested_capacity =
        (requested_capacity < curr_size) ? curr_size : requested_capacity;

    insertable_map_type tmp(requested_capacity, m_hasher, m_equal_to);

    if (curr_size) {
      tmp.m_bounded_insert = false;
      Impl::UnorderedMapRehash<insertable_map_type> f(tmp, *this);
      f.apply();
    }
    tmp.m_bounded_insert = bounded_insert;

    *this = tmp;

    return true;
  }

  /// \brief The number of entries in the table.
  ///
  /// This method has undefined behavior when erasable() is true.
  ///
  /// Note that this is <i>not</i> a device function; it cannot be called in
  /// a parallel kernel.  The value is not stored as a variable; it
  /// must be computed. m_size is a mutable cache of that value.
  size_type size() const {
    if (capacity() == 0u) return 0u;
    if (modified()) {
      m_size() = m_available_indexes.count();
      reset_flag(modified_idx);
    }
    return m_size();
  }

  /// \brief The current number of failed insert() calls.
  ///
  /// This is <i>not</i> a device function; it may <i>not</i> be
  /// called in a parallel kernel.  The value is not stored as a
  /// variable; it must be computed.
  bool failed_insert() const { return get_flag(failed_insert_idx); }

  bool erasable() const {
    return is_insertable_map ? get_flag(erasable_idx) : false;
  }

  bool begin_erase() {
    bool result = !erasable();
    if (is_insertable_map && result) {
      execution_space().fence(
          "Kokkos::UnorderedMap::begin_erase: fence before setting erasable "
          "flag");
      set_flag(erasable_idx);
    }
    return result;
  }

  bool end_erase() {
    bool result = erasable();
    if (is_insertable_map && result) {
      execution_space().fence(
          "Kokkos::UnorderedMap::end_erase: fence before erasing");
      Impl::UnorderedMapErase<declared_map_type> f(*this);
      f.apply();
      execution_space().fence(
          "Kokkos::UnorderedMap::end_erase: fence after erasing");
      reset_flag(erasable_idx);
    }
    return result;
  }

  /// \brief The maximum number of entries that the table can hold.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_FORCEINLINE_FUNCTION
  size_type capacity() const { return m_available_indexes.size(); }

  /// \brief The number of hash table "buckets."
  ///
  /// This is different than the number of entries that the table can
  /// hold.  Each key hashes to an index in [0, hash_capacity() - 1].
  /// That index can hold zero or more entries.  This class decides
  /// what hash_capacity() should be, given the user's upper bound on
  /// the number of entries the table must be able to hold.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  size_type hash_capacity() const { return m_hash_lists.extent(0); }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.  As discussed in the class documentation, it need not
  /// succeed.  The return value tells you if it did.
  ///
  /// \param k [in] The key to attempt to insert.
  /// \param v [in] The corresponding value to attempt to insert.  If
  ///   using this class as a set (with Value = void), then you need not
  ///   provide this value.
  /// \param insert_op [in] The operator used for combining values if a
  ///                       key already exists. See
  ///                       Kokkos::UnorderedMapInsertOpTypes for more ops.
  template <typename InsertOpType = default_op_type>
  KOKKOS_INLINE_FUNCTION insert_result
  insert(key_type const &k, impl_value_type const &v = impl_value_type(),
         [[maybe_unused]] InsertOpType arg_insert_op = InsertOpType()) const {
    if constexpr (is_set) {
      static_assert(std::is_same_v<InsertOpType, default_op_type>,
                    "Insert Operations are not supported on sets.");
    }

    insert_result result;

    if (!is_insertable_map || capacity() == 0u ||
        m_scalars((int)erasable_idx)) {
      return result;
    }

    if (!m_scalars((int)modified_idx)) {
      m_scalars((int)modified_idx) = true;
    }

    int volatile &failed_insert_ref = m_scalars((int)failed_insert_idx);

    const size_type hash_value = m_hasher(k);
    const size_type hash_list  = hash_value % m_hash_lists.extent(0);

    size_type *curr_ptr = &m_hash_lists[hash_list];
    size_type new_index = invalid_index;

    // Force integer multiply to long
    size_type index_hint = static_cast<size_type>(
        (static_cast<double>(hash_list) * capacity()) / m_hash_lists.extent(0));

    size_type find_attempts = 0;

    enum : unsigned { bounded_find_attempts = 32u };
    const size_type max_attempts =
        (m_bounded_insert &&
         (bounded_find_attempts < m_available_indexes.max_hint()))
            ? bounded_find_attempts
            : m_available_indexes.max_hint();

    bool not_done = true;

#if defined(__MIC__)
#pragma noprefetch
#endif
    while (not_done) {
      // Continue searching the unordered list for this key,
      // list will only be appended during insert phase.
      // Need volatile_load as other threads may be appending.

      // FIXME_SYCL replacement for memory_fence
#ifdef KOKKOS_ENABLE_SYCL
      size_type curr = Kokkos::atomic_load(curr_ptr);
#else
      size_type curr = volatile_load(curr_ptr);
#endif

      KOKKOS_NONTEMPORAL_PREFETCH_LOAD(
          &m_keys[curr != invalid_index ? curr : 0]);
#if defined(__MIC__)
#pragma noprefetch
#endif
      while (curr != invalid_index && !m_equal_to(
#ifdef KOKKOS_ENABLE_SYCL
                                          Kokkos::atomic_load(&m_keys[curr])
#else
                                          volatile_load(&m_keys[curr])
#endif
                                              ,
                                          k)) {
        result.increment_list_position();
        index_hint = curr;
        curr_ptr   = &m_next_index[curr];
#ifdef KOKKOS_ENABLE_SYCL
        curr = Kokkos::atomic_load(curr_ptr);
#else
        curr = volatile_load(curr_ptr);
#endif
        KOKKOS_NONTEMPORAL_PREFETCH_LOAD(
            &m_keys[curr != invalid_index ? curr : 0]);
      }

      //------------------------------------------------------------
      // If key already present then return that index.
      if (curr != invalid_index) {
        const bool free_existing = new_index != invalid_index;
        if (free_existing) {
          // Previously claimed an unused entry that was not inserted.
          // Release this unused entry immediately.
          if (!m_available_indexes.reset(new_index)) {
            Kokkos::printf("Unable to free existing\n");
          }
        }

        result.set_existing(curr, free_existing);
        if constexpr (!is_set) {
          arg_insert_op.op(m_values, curr, v);
        }
        not_done = false;
      }
      //------------------------------------------------------------
      // Key is not currently in the map.
      // If the thread has claimed an entry try to insert now.
      else {
        //------------------------------------------------------------
        // If have not already claimed an unused entry then do so now.
        if (new_index == invalid_index) {
          bool found = false;
          // use the hash_list as the flag for the search direction
          Kokkos::tie(found, index_hint) =
              m_available_indexes.find_any_unset_near(index_hint, hash_list);

          // found and index and this thread set it
          if (!found && ++find_attempts >= max_attempts) {
            failed_insert_ref = true;
            not_done          = false;
          } else if (m_available_indexes.set(index_hint)) {
            new_index = index_hint;
            // Set key and value
            KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_keys[new_index]);
// FIXME_SYCL replacement for memory_fence
#ifdef KOKKOS_ENABLE_SYCL
            Kokkos::atomic_store(&m_keys[new_index], k);
#else
            m_keys[new_index] = k;
#endif

            if (!is_set) {
              KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_values[new_index]);
#ifdef KOKKOS_ENABLE_SYCL
              Kokkos::atomic_store(&m_values[new_index], v);
#else
              m_values[new_index] = v;
#endif
            }

#ifndef KOKKOS_ENABLE_SYCL
            // Do not proceed until key and value are updated in global memory
            memory_fence();
#endif
          }
        } else if (failed_insert_ref) {
          not_done = false;
        }

        // Attempt to append claimed entry into the list.
        // Another thread may also be trying to append the same list so protect
        // with atomic.
        if (new_index != invalid_index &&
            curr == atomic_compare_exchange(
                        curr_ptr, static_cast<size_type>(invalid_index),
                        new_index)) {
          // Succeeded in appending
          result.set_success(new_index);
          not_done = false;
        }
      }
    }  // while ( not_done )

    return result;
  }

  KOKKOS_INLINE_FUNCTION
  bool erase(key_type const &k) const {
    bool result = false;

    if (is_insertable_map && 0u < capacity() && m_scalars((int)erasable_idx)) {
      if (!m_scalars((int)modified_idx)) {
        m_scalars((int)modified_idx) = true;
      }

      size_type index = find(k);
      if (valid_at(index)) {
        m_available_indexes.reset(index);
        result = true;
      }
    }

    return result;
  }

  /// \brief Find the given key \c k, if it exists in the table.
  ///
  /// \return If the key exists in the table, the index of the
  ///   value corresponding to that key; otherwise, an invalid index.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  size_type find(const key_type &k) const {
    size_type curr = 0u < capacity()
                         ? m_hash_lists(m_hasher(k) % m_hash_lists.extent(0))
                         : invalid_index;

    KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
    while (curr != invalid_index && !m_equal_to(m_keys[curr], k)) {
      KOKKOS_NONTEMPORAL_PREFETCH_LOAD(
          &m_keys[curr != invalid_index ? curr : 0]);
      curr = m_next_index[curr];
    }

    return curr;
  }

  /// \brief Does the key exist in the map
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  bool exists(const key_type &k) const { return valid_at(find(k)); }

  /// \brief Get the value with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  ///
  /// 'const value_type' via Cuda texture fetch must return by value.
  template <typename Dummy = value_type>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      !std::is_void<Dummy>::value,  // !is_set
      std::conditional_t<has_const_value, impl_value_type, impl_value_type &>>
  value_at(size_type i) const {
    KOKKOS_EXPECTS(i < capacity());
    return m_values[i];
  }

  /// \brief Get the key with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_FORCEINLINE_FUNCTION
  key_type key_at(size_type i) const {
    KOKKOS_EXPECTS(i < capacity());
    return m_keys[i];
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool valid_at(size_type i) const { return m_available_indexes.test(i); }

  template <typename SKey, typename SValue>
  UnorderedMap(
      UnorderedMap<SKey, SValue, Device, Hasher, EqualTo> const &src,
      std::enable_if_t<
          Impl::UnorderedMapCanAssign<declared_key_type, declared_value_type,
                                      SKey, SValue>::value,
          int> = 0)
      : m_bounded_insert(src.m_bounded_insert),
        m_hasher(src.m_hasher),
        m_equal_to(src.m_equal_to),
        m_size(src.m_size),
        m_available_indexes(src.m_available_indexes),
        m_hash_lists(src.m_hash_lists),
        m_next_index(src.m_next_index),
        m_keys(src.m_keys),
        m_values(src.m_values),
        m_scalars(src.m_scalars) {}

  template <typename SKey, typename SValue>
  std::enable_if_t<
      Impl::UnorderedMapCanAssign<declared_key_type, declared_value_type, SKey,
                                  SValue>::value,
      declared_map_type &>
  operator=(UnorderedMap<SKey, SValue, Device, Hasher, EqualTo> const &src) {
    m_bounded_insert    = src.m_bounded_insert;
    m_hasher            = src.m_hasher;
    m_equal_to          = src.m_equal_to;
    m_size              = src.m_size;
    m_available_indexes = src.m_available_indexes;
    m_hash_lists        = src.m_hash_lists;
    m_next_index        = src.m_next_index;
    m_keys              = src.m_keys;
    m_values            = src.m_values;
    m_scalars           = src.m_scalars;
    return *this;
  }

  // Re-allocate the views of the calling UnorderedMap according to src
  // capacity, and deep copy the src data.
  template <typename SKey, typename SValue, typename SDevice>
  std::enable_if_t<std::is_same<std::remove_const_t<SKey>, key_type>::value &&
                   std::is_same<std::remove_const_t<SValue>, value_type>::value>
  create_copy_view(
      UnorderedMap<SKey, SValue, SDevice, Hasher, EqualTo> const &src) {
    if (m_hash_lists.data() != src.m_hash_lists.data()) {
      allocate_view(src);
      deep_copy_view(src);
    }
  }

  // Allocate views of the calling UnorderedMap with the same capacity as the
  // src.
  template <typename SKey, typename SValue, typename SDevice>
  std::enable_if_t<std::is_same<std::remove_const_t<SKey>, key_type>::value &&
                   std::is_same<std::remove_const_t<SValue>, value_type>::value>
  allocate_view(
      UnorderedMap<SKey, SValue, SDevice, Hasher, EqualTo> const &src) {
    insertable_map_type tmp;

    tmp.m_bounded_insert    = src.m_bounded_insert;
    tmp.m_hasher            = src.m_hasher;
    tmp.m_equal_to          = src.m_equal_to;
    tmp.m_size()            = src.m_size();
    tmp.m_available_indexes = bitset_type(src.capacity());
    tmp.m_hash_lists        = size_type_view(
        view_alloc(WithoutInitializing, "UnorderedMap hash list"),
        src.m_hash_lists.extent(0));
    tmp.m_next_index = size_type_view(
        view_alloc(WithoutInitializing, "UnorderedMap next index"),
        src.m_next_index.extent(0));
    tmp.m_keys =
        key_type_view(view_alloc(WithoutInitializing, "UnorderedMap keys"),
                      src.m_keys.extent(0));
    tmp.m_values =
        value_type_view(view_alloc(WithoutInitializing, "UnorderedMap values"),
                        src.m_values.extent(0));
    tmp.m_scalars = scalars_view("UnorderedMap scalars");

    *this = tmp;
  }

  // Deep copy view data from src. This requires that the src capacity is
  // identical to the capacity of the calling UnorderedMap.
  template <typename SKey, typename SValue, typename SDevice>
  std::enable_if_t<std::is_same<std::remove_const_t<SKey>, key_type>::value &&
                   std::is_same<std::remove_const_t<SValue>, value_type>::value>
  deep_copy_view(
      UnorderedMap<SKey, SValue, SDevice, Hasher, EqualTo> const &src) {
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_4
    // To deep copy UnorderedMap, capacity must be identical
    KOKKOS_EXPECTS(capacity() == src.capacity());
#else
    if (capacity() != src.capacity()) {
      allocate_view(src);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
      Kokkos::Impl::log_warning(
          "Warning: deep_copy_view() allocating views is deprecated. Must call "
          "with UnorderedMaps of identical capacity, or use "
          "create_copy_view().\n");
#endif
    }
#endif

    if (m_hash_lists.data() != src.m_hash_lists.data()) {
      Kokkos::deep_copy(m_available_indexes, src.m_available_indexes);

      using raw_deep_copy =
          Kokkos::Impl::DeepCopy<typename device_type::memory_space,
                                 typename SDevice::memory_space>;

      raw_deep_copy(m_hash_lists.data(), src.m_hash_lists.data(),
                    sizeof(size_type) * src.m_hash_lists.extent(0));
      raw_deep_copy(m_next_index.data(), src.m_next_index.data(),
                    sizeof(size_type) * src.m_next_index.extent(0));
      raw_deep_copy(m_keys.data(), src.m_keys.data(),
                    sizeof(key_type) * src.m_keys.extent(0));
      if (!is_set) {
        raw_deep_copy(m_values.data(), src.m_values.data(),
                      sizeof(impl_value_type) * src.m_values.extent(0));
      }
      raw_deep_copy(m_scalars.data(), src.m_scalars.data(),
                    sizeof(int) * num_scalars);

      Kokkos::fence(
          "Kokkos::UnorderedMap::deep_copy_view: fence after copy to dst.");
    }
  }

  //@}
 private:  // private member functions
  bool modified() const { return get_flag(modified_idx); }

  void set_flag(int flag) const {
    using raw_deep_copy =
        Kokkos::Impl::DeepCopy<typename device_type::memory_space,
                               Kokkos::HostSpace>;
    const int true_ = true;
    raw_deep_copy(m_scalars.data() + flag, &true_, sizeof(int));
    Kokkos::fence(
        "Kokkos::UnorderedMap::set_flag: fence after copying flag from "
        "HostSpace");
  }

  void reset_flag(int flag) const {
    using raw_deep_copy =
        Kokkos::Impl::DeepCopy<typename device_type::memory_space,
                               Kokkos::HostSpace>;
    const int false_ = false;
    raw_deep_copy(m_scalars.data() + flag, &false_, sizeof(int));
    Kokkos::fence(
        "Kokkos::UnorderedMap::reset_flag: fence after copying flag from "
        "HostSpace");
  }

  bool get_flag(int flag) const {
    using raw_deep_copy =
        Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                               typename device_type::memory_space>;
    int result = false;
    raw_deep_copy(&result, m_scalars.data() + flag, sizeof(int));
    Kokkos::fence(
        "Kokkos::UnorderedMap::get_flag: fence after copy to return value in "
        "HostSpace");
    return result;
  }

  static uint32_t calculate_capacity(uint32_t capacity_hint) {
    // increase by 16% and round to nears multiple of 128
    return capacity_hint
               ? ((static_cast<uint32_t>(7ull * capacity_hint / 6u) + 127u) /
                  128u) *
                     128u
               : 128u;
  }

 private:  // private members
  bool m_bounded_insert;
  hasher_type m_hasher;
  equal_to_type m_equal_to;
  using shared_size_t = View<size_type, Kokkos::DefaultHostExecutionSpace>;
  shared_size_t m_size;
  bitset_type m_available_indexes;
  size_type_view m_hash_lists;
  size_type_view m_next_index;
  key_type_view m_keys;
  value_type_view m_values;
  scalars_view m_scalars;

  template <typename KKey, typename VValue, typename DDevice, typename HHash,
            typename EEqualTo>
  friend class UnorderedMap;

  template <typename UMap>
  friend struct Impl::UnorderedMapErase;

  template <typename UMap>
  friend struct Impl::UnorderedMapHistogram;

  template <typename UMap>
  friend struct Impl::UnorderedMapPrint;
};

// Specialization of deep_copy() for two UnorderedMap objects.
template <typename DKey, typename DT, typename DDevice, typename SKey,
          typename ST, typename SDevice, typename Hasher, typename EqualTo>
inline void deep_copy(
    UnorderedMap<DKey, DT, DDevice, Hasher, EqualTo> &dst,
    const UnorderedMap<SKey, ST, SDevice, Hasher, EqualTo> &src) {
  dst.deep_copy_view(src);
}

// Specialization of create_mirror() for an UnorderedMap object.
template <typename Key, typename ValueType, typename Device, typename Hasher,
          typename EqualTo>
typename UnorderedMap<Key, ValueType, Device, Hasher, EqualTo>::HostMirror
create_mirror(
    const UnorderedMap<Key, ValueType, Device, Hasher, EqualTo> &src) {
  typename UnorderedMap<Key, ValueType, Device, Hasher, EqualTo>::HostMirror
      dst;
  dst.allocate_view(src);
  return dst;
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_UNORDEREDMAP
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_UNORDEREDMAP
#endif
#endif  // KOKKOS_UNORDERED_MAP_HPP
