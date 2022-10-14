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

#ifndef KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP
#define KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP

#include <Kokkos_Timer.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

namespace Perf {

template <typename Device, bool Near>
struct UnorderedMapTest {
  using execution_space = Device;
  using map_type = Kokkos::UnorderedMap<uint32_t, uint32_t, execution_space>;
  using histogram_type = typename map_type::histogram_type;

  struct value_type {
    uint32_t failed_count;
    uint32_t max_list;
  };

  uint32_t capacity;
  uint32_t inserts;
  uint32_t collisions;
  double seconds;
  map_type map;
  histogram_type histogram;

  UnorderedMapTest(uint32_t arg_capacity, uint32_t arg_inserts,
                   uint32_t arg_collisions)
      : capacity(arg_capacity),
        inserts(arg_inserts),
        collisions(arg_collisions),
        seconds(0),
        map(capacity),
        histogram(map.get_histogram()) {
    Kokkos::Timer wall_clock;
    wall_clock.reset();

    value_type v   = {};
    int loop_count = 0;
    do {
      ++loop_count;

      v = value_type();
      Kokkos::parallel_reduce(inserts, *this, v);

      if (v.failed_count > 0u) {
        const uint32_t new_capacity = map.capacity() +
                                      ((map.capacity() * 3ull) / 20u) +
                                      v.failed_count / collisions;
        map.rehash(new_capacity);
      }
    } while (v.failed_count > 0u);

    seconds = wall_clock.seconds();

    switch (loop_count) {
      case 1u: std::cout << " \033[0;32m" << loop_count << "\033[0m "; break;
      case 2u: std::cout << " \033[1;31m" << loop_count << "\033[0m "; break;
      default: std::cout << " \033[0;31m" << loop_count << "\033[0m "; break;
    }
    std::cout << std::setprecision(2) << std::fixed << std::setw(5)
              << (1e9 * (seconds / (inserts))) << "; " << std::flush;

    histogram.calculate();
    Device().fence();
  }

  void print(std::ostream& metrics_out, std::ostream& length_out,
             std::ostream& distance_out, std::ostream& block_distance_out) {
    metrics_out << map.capacity() << " , ";
    metrics_out << inserts / collisions << " , ";
    metrics_out << (100.0 * inserts / collisions) / map.capacity() << " , ";
    metrics_out << inserts << " , ";
    metrics_out << (map.failed_insert() ? "true" : "false") << " , ";
    metrics_out << collisions << " , ";
    metrics_out << 1e9 * (seconds / inserts) << " , ";
    metrics_out << seconds << std::endl;

    length_out << map.capacity() << " , ";
    length_out << ((100.0 * inserts / collisions) / map.capacity()) << " , ";
    length_out << collisions << " , ";
    histogram.print_length(length_out);

    distance_out << map.capacity() << " , ";
    distance_out << ((100.0 * inserts / collisions) / map.capacity()) << " , ";
    distance_out << collisions << " , ";
    histogram.print_distance(distance_out);

    block_distance_out << map.capacity() << " , ";
    block_distance_out << ((100.0 * inserts / collisions) / map.capacity())
                       << " , ";
    block_distance_out << collisions << " , ";
    histogram.print_block_distance(block_distance_out);
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& v) const {
    v.failed_count = 0;
    v.max_list     = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const {
    dst.failed_count += src.failed_count;
    dst.max_list = src.max_list < dst.max_list ? dst.max_list : src.max_list;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type& v) const {
    const uint32_t key = Near ? i / collisions : i % (inserts / collisions);
    typename map_type::insert_result result = map.insert(key, i);
    v.failed_count += !result.failed() ? 0 : 1;
    v.max_list = result.list_position() < v.max_list ? v.max_list
                                                     : result.list_position();
  }
};

template <typename Device, bool Near>
void run_performance_tests(std::string const& base_file_name) {
#if 0
  std::string metrics_file_name = base_file_name + std::string("-metrics.csv");
  std::string length_file_name = base_file_name  + std::string("-length.csv");
  std::string distance_file_name = base_file_name + std::string("-distance.csv");
  std::string block_distance_file_name = base_file_name + std::string("-block_distance.csv");

  std::ofstream metrics_out( metrics_file_name.c_str(), std::ofstream::out );
  std::ofstream length_out( length_file_name.c_str(), std::ofstream::out );
  std::ofstream distance_out( distance_file_name.c_str(), std::ofstream::out );
  std::ofstream block_distance_out( block_distance_file_name.c_str(), std::ofstream::out );


  /*
  const double test_ratios[] = {
     0.50
   , 0.75
   , 0.80
   , 0.85
   , 0.90
   , 0.95
   , 1.00
   , 1.25
   , 2.00
  };
  */

  const double test_ratios[] = { 1.00 };

  const int num_ratios = sizeof(test_ratios) / sizeof(double);

  /*
  const uint32_t collisions[] {
      1
    , 4
    , 16
    , 64
  };
  */

  const uint32_t collisions[] = { 16 };

  const int num_collisions = sizeof(collisions) / sizeof(uint32_t);

  // set up file headers
  metrics_out << "Capacity , Unique , Percent Full , Attempted Inserts , Failed Inserts , Collision Ratio , Nanoseconds/Inserts, Seconds" << std::endl;
  length_out << "Capacity , Percent Full , ";
  distance_out << "Capacity , Percent Full , ";
  block_distance_out << "Capacity , Percent Full , ";

  for (int i=0; i<100; ++i) {
    length_out << i << " , ";
    distance_out << i << " , ";
    block_distance_out << i << " , ";
  }

  length_out << "\b\b\b   " << std::endl;
  distance_out << "\b\b\b   " << std::endl;
  block_distance_out << "\b\b\b   " << std::endl;

  Kokkos::Timer wall_clock ;
  for (int i=0;  i < num_collisions ; ++i) {
    wall_clock.reset();
    std::cout << "Collisions: " << collisions[i] << std::endl;
    for (int j = 0; j < num_ratios; ++j) {
      std::cout << std::setprecision(1) << std::fixed << std::setw(5) << (100.0*test_ratios[j]) << "%  " << std::flush;
      for (uint32_t capacity = 1<<14; capacity < 1<<25; capacity = capacity << 1) {
        uint32_t inserts = static_cast<uint32_t>(test_ratios[j]*(capacity));
        std::cout << capacity << std::flush;
        UnorderedMapTest<Device, Near> test(capacity, inserts*collisions[i], collisions[i]);
        Device().fence();
        test.print(metrics_out, length_out, distance_out, block_distance_out);
      }
      std::cout << "\b\b  " <<  std::endl;

    }
    std::cout << "  " << wall_clock.seconds() << " secs" << std::endl;
  }
  metrics_out.close();
  length_out.close();
  distance_out.close();
  block_distance_out.close();
#else
  (void)base_file_name;
  std::cout << "skipping test" << std::endl;
#endif
}

}  // namespace Perf

#endif  // KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP
