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

#ifndef KOKKOS_TEST_UNORDERED_MAP_HPP
#define KOKKOS_TEST_UNORDERED_MAP_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <Kokkos_UnorderedMap.hpp>

namespace Test {

namespace Impl {

template <typename MapType, bool Near = false>
struct TestInsert
{
  typedef MapType map_type;
  typedef typename map_type::execution_space execution_space;
  typedef uint32_t value_type;

  map_type map;
  uint32_t inserts;
  uint32_t collisions;

  TestInsert( map_type arg_map, uint32_t arg_inserts, uint32_t arg_collisions)
    : map(arg_map)
    , inserts(arg_inserts)
    , collisions(arg_collisions)
  {}

  void testit( bool rehash_on_fail = true )
  {
    execution_space::fence();

    uint32_t failed_count = 0;
    do {
      failed_count = 0;
      Kokkos::parallel_reduce(inserts, *this, failed_count);

      if (rehash_on_fail && failed_count > 0u) {
        const uint32_t new_capacity = map.capacity() + ((map.capacity()*3ull)/20u) + failed_count/collisions ;
        map.rehash( new_capacity );
      }
    } while (rehash_on_fail && failed_count > 0u);

    execution_space::fence();
  }


  KOKKOS_INLINE_FUNCTION
  void init( value_type & failed_count ) const { failed_count = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & failed_count, const volatile value_type & count ) const
  { failed_count += count; }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type & failed_count) const
  {
    const uint32_t key = Near ? i/collisions : i%(inserts/collisions);
    if (map.insert(key,i).failed()) ++failed_count;
  }

};

  template <typename MapType, bool Near>
  struct TestErase
  {
    typedef TestErase<MapType, Near> self_type;

    typedef MapType map_type;
    typedef typename MapType::execution_space execution_space;

    map_type m_map;
    uint32_t m_num_erase;
    uint32_t m_num_duplicates;

    TestErase(map_type map, uint32_t num_erases, uint32_t num_duplicates)
      : m_map(map)
      , m_num_erase(num_erases)
      , m_num_duplicates(num_duplicates)
    {}

    void testit()
    {
      execution_space::fence();
      Kokkos::parallel_for(m_num_erase, *this);
      execution_space::fence();
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename execution_space::size_type i) const
    {
      if (Near) {
        m_map.erase(i/m_num_duplicates);
      }
      else {
        m_map.erase(i%(m_num_erase/m_num_duplicates));
      }

    }
  };

  template <typename MapType>
  struct TestFind
  {
    typedef MapType map_type;
    typedef typename MapType::execution_space::execution_space execution_space;
    typedef uint32_t value_type;

    map_type m_map;
    uint32_t m_num_insert;
    uint32_t m_num_duplicates;
    uint32_t m_max_key;

    TestFind(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map)
      , m_num_insert(num_inserts)
      , m_num_duplicates(num_duplicates)
      , m_max_key( ((num_inserts + num_duplicates) - 1)/num_duplicates )
    {}

    void testit(value_type &errors)
    {
      execution_space::execution_space::fence();
      Kokkos::parallel_reduce(m_map.capacity(), *this, errors);
      execution_space::execution_space::fence();
    }

    KOKKOS_INLINE_FUNCTION
    static void init( value_type & dst)
    {
      dst = 0;
    }

    KOKKOS_INLINE_FUNCTION
    static void join( volatile value_type & dst, const volatile value_type & src)
    { dst += src; }

    KOKKOS_INLINE_FUNCTION
    void operator()(typename execution_space::size_type i, value_type & errors) const
    {
      const bool expect_to_find_i = (i < m_max_key);

      const bool exists = m_map.exists(i);

      if (expect_to_find_i && !exists)  ++errors;
      if (!expect_to_find_i && exists)  ++errors;
    }
  };

} // namespace Impl



template <typename Device>
void test_insert( uint32_t num_nodes , uint32_t num_inserts , uint32_t num_duplicates , bool near )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t,const uint32_t, Device> const_map_type;

  const uint32_t expected_inserts = (num_inserts + num_duplicates -1u) / num_duplicates;

  map_type map;
  map.rehash(num_nodes,false);

  if (near) {
    Impl::TestInsert<map_type,true> test_insert(map, num_inserts, num_duplicates);
    test_insert.testit();
  } else
  {
    Impl::TestInsert<map_type,false> test_insert(map, num_inserts, num_duplicates);
    test_insert.testit();
  }

  const bool print_list = false;
  if (print_list) {
    Kokkos::Impl::UnorderedMapPrint<map_type> f(map);
    f.apply();
  }

  const uint32_t map_size = map.size();

  ASSERT_FALSE( map.failed_insert());
  {
    EXPECT_EQ(expected_inserts, map_size);

    {
      uint32_t find_errors = 0;
      Impl::TestFind<const_map_type> test_find(map, num_inserts, num_duplicates);
      test_find.testit(find_errors);
      EXPECT_EQ( 0u, find_errors);
    }

    map.begin_erase();
    Impl::TestErase<map_type,false> test_erase(map, num_inserts, num_duplicates);
    test_erase.testit();
    map.end_erase();
    EXPECT_EQ(0u, map.size());
  }
}

template <typename Device>
void test_failed_insert( uint32_t num_nodes)
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;

  map_type map(num_nodes);
  Impl::TestInsert<map_type> test_insert(map, 2u*num_nodes, 1u);
  test_insert.testit(false /*don't rehash on fail*/);
  Device::execution_space::fence();

  EXPECT_TRUE( map.failed_insert() );
}



template <typename Device>
void test_deep_copy( uint32_t num_nodes )
{
  typedef Kokkos::UnorderedMap<uint32_t,uint32_t, Device> map_type;
  typedef Kokkos::UnorderedMap<const uint32_t, const uint32_t, Device> const_map_type;

  typedef typename map_type::HostMirror host_map_type ;
  // typedef Kokkos::UnorderedMap<uint32_t, uint32_t, typename Device::host_mirror_execution_space > host_map_type;

  map_type map;
  map.rehash(num_nodes,false);

  {
    Impl::TestInsert<map_type> test_insert(map, num_nodes, 1);
    test_insert.testit();
    ASSERT_EQ( map.size(), num_nodes);
    ASSERT_FALSE( map.failed_insert() );
    {
      uint32_t find_errors = 0;
      Impl::TestFind<map_type> test_find(map, num_nodes, 1);
      test_find.testit(find_errors);
      EXPECT_EQ( find_errors, 0u);
    }

  }

  host_map_type hmap;
  Kokkos::deep_copy(hmap, map);

  ASSERT_EQ( map.size(), hmap.size());
  ASSERT_EQ( map.capacity(), hmap.capacity());
  {
    uint32_t find_errors = 0;
    Impl::TestFind<host_map_type> test_find(hmap, num_nodes, 1);
    test_find.testit(find_errors);
    EXPECT_EQ( find_errors, 0u);
  }

  map_type mmap;
  Kokkos::deep_copy(mmap, hmap);

  const_map_type cmap = mmap;

  EXPECT_EQ( cmap.size(), num_nodes);

  {
    uint32_t find_errors = 0;
    Impl::TestFind<const_map_type> test_find(cmap, num_nodes, 1);
    test_find.testit(find_errors);
    EXPECT_EQ( find_errors, 0u);
  }

}

TEST_F( TEST_CATEGORY, UnorderedMap_insert) {
  for (int i=0; i<500; ++i) {
    test_insert<TEST_EXECSPACE>(100000, 90000, 100, true);
    test_insert<TEST_EXECSPACE>(100000, 90000, 100, false);
  }
}

TEST_F( TEST_CATEGORY, UnorderedMap_failed_insert) {
  for (int i=0; i<1000; ++i)
    test_failed_insert<TEST_EXECSPACE>(10000);
}

TEST_F( TEST_CATEGORY, UnorderedMap_deep_copy) {
  for (int i=0; i<2; ++i)
    test_deep_copy<TEST_EXECSPACE>(10000);
}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP

