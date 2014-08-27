/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

/// \file Kokkos_UnorderedMap.hpp
/// \brief Declaration and definition of Kokkos::UnorderedMap.
///
/// This header file declares and defines Kokkos::UnorderedMap and its
/// related nonmember functions.

#ifndef KOKKOS_UNORDERED_MAP_HPP
#define KOKKOS_UNORDERED_MAP_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Functional.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_HostSpace.hpp>

#include <Kokkos_Bitset.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_UnorderedMap_impl.hpp>


#include <iostream>

#include <stdint.h>
#include <stdexcept>


namespace Kokkos {

enum { UnorderedMapInvalidIndex = ~0u };

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

class UnorderedMapInsertResult
{
private:
  enum Status{
     SUCCESS = 1u << 31
   , EXISTING = 1u << 30
   , FREED_EXISTING = 1u << 29
   , LIST_LENGTH_MASK = ~(SUCCESS | EXISTING | FREED_EXISTING)
  };

public:
  /// Did the map successful insert the key/value pair
  KOKKOS_FORCEINLINE_FUNCTION
  bool success() const { return (m_status & SUCCESS); }

  /// Was the key already present in the map
  KOKKOS_FORCEINLINE_FUNCTION
  bool existing() const { return (m_status & EXISTING); }

  /// Did the map fail to insert the key due to insufficent capacity
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
  UnorderedMapInsertResult()
    : m_index(UnorderedMapInvalidIndex)
    , m_status(0)
  {}

  KOKKOS_FORCEINLINE_FUNCTION
  void increment_list_position()
  {
    m_status += (list_position() < LIST_LENGTH_MASK) ? 1u : 0u;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_existing(uint32_t i, bool arg_freed_existing)
  {
    m_index = i;
    m_status = EXISTING | (arg_freed_existing ? FREED_EXISTING : 0u) | list_position();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_success(uint32_t i)
  {
    m_index = i;
    m_status = SUCCESS | list_position();
  }

private:
  uint32_t m_index;
  uint32_t m_status;
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
template <   typename Key
           , typename Value
           , typename Device = Kokkos::DefaultExecutionSpace
           , typename Hasher = pod_hash<typename Impl::remove_const<Key>::type>
           , typename EqualTo = pod_equal_to<typename Impl::remove_const<Key>::type>
        >
class UnorderedMap
{
public:
  //! \name Public types and constants
  //@{

  //key_types
  typedef Key declared_key_type;
  typedef typename Impl::remove_const<declared_key_type>::type key_type;
  typedef typename Impl::add_const<key_type>::type const_key_type;

  //value_types
  typedef Value declared_value_type;
  typedef typename Impl::remove_const<declared_value_type>::type value_type;
  typedef typename Impl::add_const<value_type>::type const_value_type;

  typedef Device device_type;
  typedef Hasher hasher_type;
  typedef EqualTo  equal_to_type;
  typedef uint32_t size_type;

  //map_types
  typedef UnorderedMap<declared_key_type,declared_value_type,device_type,hasher_type,equal_to_type> declared_map_type;
  typedef UnorderedMap<key_type,value_type,device_type,hasher_type,equal_to_type>                   insertable_map_type;
  typedef UnorderedMap<const_key_type,value_type,device_type,hasher_type,equal_to_type>             modifiable_map_type;
  typedef UnorderedMap<const_key_type,const_value_type,device_type,hasher_type,equal_to_type>       const_map_type;

  static const bool is_set = Impl::is_same<void,value_type>::value;
  static const bool has_const_key = Impl::is_same<const_key_type,declared_key_type>::value;
  static const bool has_const_value = is_set || Impl::is_same<const_value_type,declared_value_type>::value;

  static const bool is_insertable_map = !has_const_key && (is_set || !has_const_value);
  static const bool is_modifiable_map = has_const_key && !has_const_value;
  static const bool is_const_map = has_const_key && has_const_value;


  typedef UnorderedMapInsertResult insert_result;

  typedef typename Device::host_mirror_device_type host_mirror_device_type;

  typedef UnorderedMap<Key,Value,host_mirror_device_type,Hasher,EqualTo> HostMirror;

  typedef Impl::UnorderedMapHistogram<const_map_type> histogram_type;

  //@}

private:
  enum { invalid_index = ~static_cast<size_type>(0) };

  typedef typename Impl::if_c< is_set, int, declared_value_type>::type impl_value_type;

  typedef typename Impl::if_c<   is_insertable_map
                               , View< key_type *, device_type>
                               , View< const key_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type key_type_view;

  typedef typename Impl::if_c<   is_insertable_map || is_modifiable_map
                               , View< impl_value_type *, device_type>
                               , View< const impl_value_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type value_type_view;

  typedef typename Impl::if_c<   is_insertable_map
                               , View< size_type *, device_type>
                               , View< const size_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type size_type_view;

  typedef typename Impl::if_c<   is_insertable_map
                               , Bitset< device_type >
                               , ConstBitset< device_type>
                             >::type bitset_type;

  enum { modified_idx = 0, erasable_idx = 1, failed_insert_idx = 2 };
  enum { num_scalars = 3 };
  typedef View< int[num_scalars], LayoutLeft, device_type> scalars_view;

public:
  //! \name Public member functions
  //@{

  UnorderedMap()
    : m_bounded_insert()
    , m_hasher()
    , m_equal_to()
    , m_size()
    , m_available_indexes()
    , m_hash_lists()
    , m_next_index()
    , m_keys()
    , m_values()
    , m_scalars()
  {}

  /// \brief Constructor
  ///
  /// \param capacity_hint [in] Initial guess of how many unique keys will be inserted into the map
  /// \param hash [in] Hasher function for \c Key instances.  The
  ///   default value usually suffices.
  UnorderedMap(  size_type capacity_hint, hasher_type hasher = hasher_type(), equal_to_type equal_to = equal_to_type() )
    : m_bounded_insert(true)
    , m_hasher(hasher)
    , m_equal_to(equal_to)
    , m_size()
    , m_available_indexes(calculate_capacity(capacity_hint))
    , m_hash_lists(AllocateWithoutInitializing(), "UnorderedMap hash list", Impl::find_hash_size(capacity()))
    , m_next_index(AllocateWithoutInitializing(), "UnorderedMap next index", capacity()+1) // +1 so that the *_at functions can always return a valid reference
    , m_keys("UnorderedMap keys",capacity()+1)
    , m_values("UnorderedMap values",(is_set? 1 : capacity()+1))
    , m_scalars("UnorderedMap scalars")
  {
    if (!is_insertable_map) {
      throw std::runtime_error("Cannot construct a non-insertable (i.e. const key_type) unordered_map");
    }

    Kokkos::deep_copy(m_hash_lists, invalid_index);
    Kokkos::deep_copy(m_next_index, invalid_index);
  }

  void reset_failed_insert_flag()
  {
    reset_flag(failed_insert_idx);
  }

  histogram_type get_histogram()
  {
    return histogram_type(*this);
  }

  //! Clear all entries in the table.
  void clear()
  {
    m_bounded_insert = true;

    if (capacity() == 0) return;

    m_available_indexes.clear();

    Kokkos::deep_copy(m_hash_lists, invalid_index);
    Kokkos::deep_copy(m_next_index, invalid_index);
    {
      const key_type tmp = key_type();
      Kokkos::deep_copy(m_keys,tmp);
    }
    if (is_set){
      const impl_value_type tmp = impl_value_type();
      Kokkos::deep_copy(m_values,tmp);
    }
    {
      Kokkos::deep_copy(m_scalars, 0);
    }
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
  bool rehash(size_type requested_capacity = 0)
  {
    const bool bounded_insert = (capacity() == 0) || (size() == 0u);
    return rehash(requested_capacity, bounded_insert );
  }

  bool rehash(size_type requested_capacity, bool bounded_insert)
  {
    if(!is_insertable_map) return false;

    const size_type curr_size = size();
    requested_capacity = (requested_capacity < curr_size) ? curr_size : requested_capacity;

    insertable_map_type tmp(requested_capacity, m_hasher, m_equal_to);

    if (curr_size) {
      tmp.m_bounded_insert = false;
      Impl::UnorderedMapRehash<insertable_map_type> f(tmp,*this);
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
  /// Note that this is not a device function; it cannot be called in
  /// a parallel kernel.  The value is not stored as a variable; it
  /// must be computed.
  size_type size() const
  {
    if( capacity() == 0u ) return 0u;
    if (modified()) {
      m_size = m_available_indexes.count();
      reset_flag(modified_idx);
    }
    return m_size;
  }

  /// \brief The current number of failed insert() calls.
  ///
  /// This is <i>not</i> a device function; it may <i>not</i> be
  /// called in a parallel kernel.  The value is not stored as a
  /// variable; it must be computed.
  bool failed_insert() const
  {
    return get_flag(failed_insert_idx);
  }

  bool erasable() const
  {
    return is_insertable_map ? get_flag(erasable_idx) : false;
  }

  bool begin_erase()
  {
    bool result = !erasable();
    if (is_insertable_map && result) {
      device_type::fence();
      set_flag(erasable_idx);
      device_type::fence();
    }
    return result;
  }

  bool end_erase()
  {
    bool result = erasable();
    if (is_insertable_map && result) {
      device_type::fence();
      Impl::UnorderedMapErase<declared_map_type> f(*this);
      f.apply();
      device_type::fence();
      reset_flag(erasable_idx);
    }
    return result;
  }

  /// \brief The maximum number of entries that the table can hold.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_FORCEINLINE_FUNCTION
  size_type capacity() const
  { return m_available_indexes.size(); }

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
  size_type hash_capacity() const
  { return m_hash_lists.size(); }

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
  KOKKOS_INLINE_FUNCTION
  insert_result insert(key_type const& k, impl_value_type const&v = impl_value_type()) const
  {
    insert_result result;

    if ( !is_insertable_map || capacity() == 0u || m_scalars((int)erasable_idx) ) {
      return result;
    }

    if ( !m_scalars((int)modified_idx) ) {
      m_scalars((int)modified_idx) = true;
    }

    int volatile & failed_insert_ref = m_scalars((int)failed_insert_idx) ;

    const size_type hash_value = m_hasher(k);
    const size_type hash_list = hash_value % m_hash_lists.size();

    size_type * curr_ptr   = & m_hash_lists[ hash_list ];
    size_type new_index    = invalid_index ;

    // Force integer multiply to long
    size_type index_hint = static_cast<size_type>( (static_cast<double>(hash_list) * capacity()) / m_hash_lists.size());

    size_type find_attempts = 0;

    enum { bounded_find_attempts = 32u };
    const size_type max_attempts = (m_bounded_insert && (bounded_find_attempts < m_available_indexes.max_hint()) ) ?
                                    bounded_find_attempts :
                                    m_available_indexes.max_hint();

    bool not_done = true ;

#if defined( __MIC__ )
      #pragma noprefetch
#endif
    while ( not_done ) {

      // Continue searching the unordered list for this key,
      // list will only be appended during insert phase.
      // Need volatile_load as other threads may be appending.
      size_type curr = volatile_load(curr_ptr);

      KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
#if defined( __MIC__ )
      #pragma noprefetch
#endif
      while ( curr != invalid_index && ! m_equal_to( volatile_load(&m_keys[curr]), k) ) {
        result.increment_list_position();
        index_hint = curr;
        curr_ptr = &m_next_index[curr];
        curr = volatile_load(curr_ptr);
        KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
      }

      //------------------------------------------------------------
      // If key already present then return that index.
      if ( curr != invalid_index ) {

        const bool free_existing = new_index != invalid_index;
        if ( free_existing ) {
          // Previously claimed an unused entry that was not inserted.
          // Release this unused entry immediately.
          if (!m_available_indexes.reset(new_index) ) {
            printf("Unable to free existing\n");
          }

        }

        result.set_existing(curr, free_existing);
        not_done = false ;
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
          Kokkos::tie(found, index_hint) = m_available_indexes.find_any_unset_near( index_hint, hash_list );

          // found and index and this thread set it
          if ( !found && ++find_attempts >= max_attempts ) {
            failed_insert_ref = true;
            not_done = false ;
          }
          else if (m_available_indexes.set(index_hint) ) {
            new_index = index_hint;
            // Set key and value
            KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_keys[new_index]);
            m_keys[new_index] = k ;

            if (!is_set) {
              KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_values[new_index]);
              m_values[new_index] = v ;
            }

            // Do not proceed until key and value are updated in global memory
            memory_fence();
          }
        }
        else if (failed_insert_ref) {
          not_done = false;
        }

        // Attempt to append claimed entry into the list.
        // Another thread may also be trying to append the same list so protect with atomic.
        if ( new_index != invalid_index &&
             curr ==  atomic_compare_exchange(curr_ptr, static_cast<size_type>(invalid_index), new_index) ) {
          // Succeeded in appending
          result.set_success(new_index);
          not_done = false ;
        }
      }
    } // while ( not_done )

    return result ;
  }

  KOKKOS_INLINE_FUNCTION
  bool erase(key_type const& k) const
  {
    bool result = false;

    if(is_insertable_map && 0u < capacity() && m_scalars((int)erasable_idx)) {

      if ( ! m_scalars((int)modified_idx) ) {
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
  size_type find( const key_type & k) const
  {
    size_type curr = 0u < capacity() ? m_hash_lists( m_hasher(k) % m_hash_lists.size() ) : invalid_index ;

    KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
    while (curr != invalid_index && !m_equal_to( m_keys[curr], k) ) {
      KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
      curr = m_next_index[curr];
    }

    return curr;
  }

  /// \brief Does the key exist in the map
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  bool exists( const key_type & k) const
  {
    return valid_at(find(k));
  }


  /// \brief Get the value with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  ///
  /// 'const value_type' via Cuda texture fetch must return by value.
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::if_c< (is_set || has_const_value), impl_value_type, impl_value_type &>::type
  value_at(size_type i) const
  {
    return m_values[ is_set ? 0 : (i < capacity() ? i : capacity()) ];
  }

  /// \brief Get the key with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_FORCEINLINE_FUNCTION
  key_type key_at(size_type i) const
  {
    return m_keys[ i < capacity() ? i : capacity() ];
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool valid_at(size_type i) const
  {
    return m_available_indexes.test(i);
  }

  template <typename SKey, typename SValue>
  UnorderedMap( UnorderedMap<SKey,SValue,Device,Hasher,EqualTo> const& src,
                typename Impl::enable_if< Impl::UnorderedMapCanAssign<declared_key_type,declared_value_type,SKey,SValue>::value,int>::type = 0
              )
    : m_bounded_insert(src.m_bounded_insert)
    , m_hasher(src.m_hasher)
    , m_equal_to(src.m_equal_to)
    , m_size(src.m_size)
    , m_available_indexes(src.m_available_indexes)
    , m_hash_lists(src.m_hash_lists)
    , m_next_index(src.m_next_index)
    , m_keys(src.m_keys)
    , m_values(src.m_values)
    , m_scalars(src.m_scalars)
  {}


  template <typename SKey, typename SValue>
  typename Impl::enable_if< Impl::UnorderedMapCanAssign<declared_key_type,declared_value_type,SKey,SValue>::value
                           ,declared_map_type & >::type
  operator=( UnorderedMap<SKey,SValue,Device,Hasher,EqualTo> const& src)
  {
    m_bounded_insert = src.m_bounded_insert;
    m_hasher = src.m_hasher;
    m_equal_to = src.m_equal_to;
    m_size = src.m_size;
    m_available_indexes = src.m_available_indexes;
    m_hash_lists = src.m_hash_lists;
    m_next_index = src.m_next_index;
    m_keys = src.m_keys;
    m_values = src.m_values;
    m_scalars = src.m_scalars;
    return *this;
  }

  template <typename SKey, typename SValue, typename SDevice>
  typename Impl::enable_if< Impl::is_same< typename Impl::remove_const<SKey>::type, key_type>::value &&
                            Impl::is_same< typename Impl::remove_const<SValue>::type, value_type>::value
                          >::type
  create_copy_view( UnorderedMap<SKey, SValue, SDevice, Hasher,EqualTo> const& src)
  {
    if (m_hash_lists.ptr_on_device() != src.m_hash_lists.ptr_on_device()) {

      insertable_map_type tmp;

      tmp.m_bounded_insert = src.m_bounded_insert;
      tmp.m_hasher = src.m_hasher;
      tmp.m_equal_to = src.m_equal_to;
      tmp.m_size = src.size();
      tmp.m_available_indexes = bitset_type( src.capacity() );
      tmp.m_hash_lists        = size_type_view( AllocateWithoutInitializing(), "UnorderedMap hash list", src.m_hash_lists.size() );
      tmp.m_next_index        = size_type_view( AllocateWithoutInitializing(), "UnorderedMap next index", src.m_next_index.size() );
      tmp.m_keys              = key_type_view( AllocateWithoutInitializing(), "UnorderedMap keys", src.m_keys.size() );
      tmp.m_values            = value_type_view( AllocateWithoutInitializing(), "UnorderedMap values", src.m_values.size() );
      tmp.m_scalars           = scalars_view("UnorderedMap scalars");

      Kokkos::deep_copy(tmp.m_available_indexes, src.m_available_indexes);

      typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, typename SDevice::memory_space > raw_deep_copy;

      raw_deep_copy(tmp.m_hash_lists.ptr_on_device(), src.m_hash_lists.ptr_on_device(), sizeof(size_type)*src.m_hash_lists.size());
      raw_deep_copy(tmp.m_next_index.ptr_on_device(), src.m_next_index.ptr_on_device(), sizeof(size_type)*src.m_next_index.size());
      raw_deep_copy(tmp.m_keys.ptr_on_device(), src.m_keys.ptr_on_device(), sizeof(key_type)*src.m_keys.size());
      if (!is_set) {
        raw_deep_copy(tmp.m_values.ptr_on_device(), src.m_values.ptr_on_device(), sizeof(impl_value_type)*src.m_values.size());
      }
      raw_deep_copy(tmp.m_scalars.ptr_on_device(), src.m_scalars.ptr_on_device(), sizeof(int)*num_scalars );

      *this = tmp;
    }
  }

  //@}
private: // private member functions

  bool modified() const
  {
    return get_flag(modified_idx);
  }

  void set_flag(int flag) const
  {
    typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > raw_deep_copy;
    const int true_ = true;
    raw_deep_copy(m_scalars.ptr_on_device() + flag, &true_, sizeof(int));
  }

  void reset_flag(int flag) const
  {
    typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > raw_deep_copy;
    const int false_ = false;
    raw_deep_copy(m_scalars.ptr_on_device() + flag, &false_, sizeof(int));
  }

  bool get_flag(int flag) const
  {
    typedef Kokkos::Impl::DeepCopy< Kokkos::HostSpace, typename device_type::memory_space > raw_deep_copy;
    int result = false;
    raw_deep_copy(&result, m_scalars.ptr_on_device() + flag, sizeof(int));
    return result;
  }

  static uint32_t calculate_capacity(uint32_t capacity_hint)
  {
    // increase by 16% and round to nears multiple of 128
    return capacity_hint ? ((static_cast<uint32_t>(7ull*capacity_hint/6u) + 127u)/128u)*128u : 128u;
  }

private: // private members
  bool              m_bounded_insert;
  hasher_type       m_hasher;
  equal_to_type     m_equal_to;
  mutable size_type m_size;
  bitset_type       m_available_indexes;
  size_type_view    m_hash_lists;
  size_type_view    m_next_index;
  key_type_view     m_keys;
  value_type_view   m_values;
  scalars_view      m_scalars;

  template <typename KKey, typename VValue, typename DDevice, typename HHash, typename EEqualTo>
  friend class UnorderedMap;

  template <typename UMap>
  friend struct Impl::UnorderedMapErase;

  template <typename UMap>
  friend struct Impl::UnorderedMapHistogram;

  template <typename UMap>
  friend struct Impl::UnorderedMapPrint;
};

// Specialization of deep_copy for two UnorderedMap objects.
template <  typename DKey, typename DT, typename DDevice
          , typename SKey, typename ST, typename SDevice
          , typename Hasher, typename EqualTo >
inline void deep_copy(         UnorderedMap<DKey, DT, DDevice, Hasher, EqualTo> & dst
                       , const UnorderedMap<SKey, ST, SDevice, Hasher, EqualTo> & src )
{
  dst.create_copy_view(src);
}


} // namespace Kokkos

#endif //KOKKOS_UNORDERED_MAP_HPP
