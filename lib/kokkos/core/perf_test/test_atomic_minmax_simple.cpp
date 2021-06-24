// export OMP_PROC_BIND=spread ; export OMP_PLACES=threads
// c++  -O2 -g -DNDEBUG  -fopenmp
// ../core/perf_test/test_atomic_minmax_simple.cpp -I../core/src/ -I. -o
// test_atomic_minmax_simple.x  containers/src/libkokkoscontainers.a
// core/src/libkokkoscore.a -ldl && OMP_NUM_THREADS=1
// ./test_atomic_minmax_simple.x 10000000

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <typeinfo>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

using exec_space = Kokkos::DefaultExecutionSpace;

template <typename T>
void test(const int length) {
  Kokkos::Impl::Timer timer;

  using vector = Kokkos::View<T*, exec_space>;

  vector inp("input", length);
  T max = std::numeric_limits<T>::max();
  T min = std::numeric_limits<T>::lowest();

  // input is max values - all min atomics will replace
  {
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) { inp(i) = max; });
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) {
          (void)Kokkos::atomic_fetch_min(&(inp(i)), (T)i);
        });
    Kokkos::fence();
    double time = timer.seconds();

    int errors(0);
    Kokkos::parallel_reduce(
        length,
        KOKKOS_LAMBDA(const int i, int& inner) { inner += (inp(i) != (T)i); },
        errors);
    Kokkos::fence();

    if (errors) {
      std::cerr << "Error in 100% min replacements: " << errors << std::endl;
      std::cerr << "inp(0)=" << inp(0) << std::endl;
    }
    std::cout << "Time for 100% min replacements: " << time << std::endl;
  }

  // input is min values - all max atomics will replace
  {
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) { inp(i) = min; });
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) {
          (void)Kokkos::atomic_max_fetch(&(inp(i)), (T)i);
        });
    Kokkos::fence();
    double time = timer.seconds();

    int errors(0);
    Kokkos::parallel_reduce(
        length,
        KOKKOS_LAMBDA(const int i, int& inner) { inner += (inp(i) != (T)i); },
        errors);
    Kokkos::fence();

    if (errors) {
      std::cerr << "Error in 100% max replacements: " << errors << std::endl;
      std::cerr << "inp(0)=" << inp(0) << std::endl;
    }
    std::cout << "Time for 100% max replacements: " << time << std::endl;
  }

  // input is max values - all max atomics will early exit
  {
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) { inp(i) = max; });
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) {
          (void)Kokkos::atomic_max_fetch(&(inp(i)), (T)i);
        });
    Kokkos::fence();
    double time = timer.seconds();

    int errors(0);
    Kokkos::parallel_reduce(
        length,
        KOKKOS_LAMBDA(const int i, int& inner) {
          T ref = max;
          inner += (inp(i) != ref);
        },
        errors);
    Kokkos::fence();

    if (errors) {
      std::cerr << "Error in 100% max early exits: " << errors << std::endl;
      std::cerr << "inp(0)=" << inp(0) << std::endl;
    }
    std::cout << "Time for 100% max early exits: " << time << std::endl;
  }

  // input is min values - all min atomics will early exit
  {
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) { inp(i) = min; });
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(
        length, KOKKOS_LAMBDA(const int i) {
          (void)Kokkos::atomic_min_fetch(&(inp(i)), (T)i);
        });
    Kokkos::fence();
    double time = timer.seconds();

    int errors(0);
    Kokkos::parallel_reduce(
        length,
        KOKKOS_LAMBDA(const int i, int& inner) {
          T ref = min;
          inner += (inp(i) != ref);
        },
        errors);
    Kokkos::fence();

    if (errors) {
      std::cerr << "Error in 100% min early exits: " << errors << std::endl;
      std::cerr << "inp(0)=" << inp(0) << std::endl;
      if (length > 9) std::cout << "inp(9)=" << inp(9) << std::endl;
    }
    std::cout << "Time for 100% min early exits: " << time << std::endl;
  }

  // limit iterations for contentious test, takes ~50x longer for same length
  auto con_length = length / 5;
  // input is min values - some max atomics will replace
  {
    Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) { inp(i) = min; });
    Kokkos::fence();

    T current(0);
    timer.reset();
    Kokkos::parallel_reduce(
        con_length,
        KOKKOS_LAMBDA(const int i, T& inner) {
          inner = Kokkos::atomic_max_fetch(&(inp(0)), inner + 1);
          if (i == con_length - 1) {
            Kokkos::atomic_max_fetch(&(inp(0)), max);
            inner = max;
          }
        },
        Kokkos::Max<T>(current));
    Kokkos::fence();
    double time = timer.seconds();

    if (current < max) {
      std::cerr << "Error in contentious max replacements: " << std::endl;
      std::cerr << "final=" << current << " inp(0)=" << inp(0) << " max=" << max
                << std::endl;
    }
    std::cout << "Time for contentious max " << con_length
              << " replacements: " << time << std::endl;
  }

  // input is max values - some min atomics will replace
  {
    Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) { inp(i) = max; });
    Kokkos::fence();

    timer.reset();
    T current(100000000);
    Kokkos::parallel_reduce(
        con_length,
        KOKKOS_LAMBDA(const int i, T& inner) {
          inner = Kokkos::atomic_min_fetch(&(inp(0)), inner - 1);
          if (i == con_length - 1) {
            Kokkos::atomic_min_fetch(&(inp(0)), min);
            inner = min;
          }
        },
        Kokkos::Min<T>(current));
    Kokkos::fence();
    double time = timer.seconds();

    if (current > min) {
      std::cerr << "Error in contentious min replacements: " << std::endl;
      std::cerr << "final=" << current << " inp(0)=" << inp(0) << " min=" << min
                << std::endl;
    }
    std::cout << "Time for contentious min " << con_length
              << " replacements: " << time << std::endl;
  }
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    int length = 1000000;
    if (argc == 2) {
      length = std::stoi(argv[1]);
    }

    if (length < 1) {
      throw std::invalid_argument("");
    }

    std::cout << "================ int" << std::endl;
    test<int>(length);
    std::cout << "================ long" << std::endl;
    test<long>(length);
    std::cout << "================ long long" << std::endl;
    test<long long>(length);

    std::cout << "================ unsigned int" << std::endl;
    test<unsigned int>(length);
    std::cout << "================ unsigned long" << std::endl;
    test<unsigned long>(length);
    std::cout << "================ unsigned long long" << std::endl;
    test<unsigned long long>(length);

    std::cout << "================ float" << std::endl;
    test<float>(length);
    std::cout << "================ double" << std::endl;
    test<double>(length);
  }
  Kokkos::finalize();
  return 0;
}
