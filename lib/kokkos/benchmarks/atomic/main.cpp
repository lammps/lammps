#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Random.hpp>

template <class Scalar>
double test_atomic(int L, int N, int M, int K, int R,
                   Kokkos::View<const int*> offsets) {
  Kokkos::View<Scalar*> output("Output", N);
  Kokkos::Impl::Timer timer;

  for (int r = 0; r < R; r++)
    Kokkos::parallel_for(
        L, KOKKOS_LAMBDA(const int& i) {
          Scalar s = 2;
          for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) s = s * s + s;
            const int idx = (i + offsets(i, m)) % N;
            Kokkos::atomic_add(&output(idx), s);
          }
        });
  Kokkos::fence();
  double time = timer.seconds();

  return time;
}

template <class Scalar>
double test_no_atomic(int L, int N, int M, int K, int R,
                      Kokkos::View<const int*> offsets) {
  Kokkos::View<Scalar*> output("Output", N);
  Kokkos::Impl::Timer timer;
  for (int r = 0; r < R; r++)
    Kokkos::parallel_for(
        L, KOKKOS_LAMBDA(const int& i) {
          Scalar s = 2;
          for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) s = s * s + s;
            const int idx = (i + offsets(i, m)) % N;
            output(idx) += s;
          }
        });
  Kokkos::fence();
  double time = timer.seconds();
  return time;
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    if (argc < 8) {
      printf("Arguments: L N M D K R T\n");
      printf("  L:   Number of iterations to run\n");
      printf("  N:   Length of array to do atomics into\n");
      printf("  M:   Number of atomics per iteration to do\n");
      printf("  D:   Distance from index i to do atomics into (randomly)\n");
      printf("  K:   Number of FMAD per atomic\n");
      printf("  R:   Number of repeats of the experiments\n");
      printf("  T:   Type of atomic\n");
      printf("       1 - int\n");
      printf("       2 - long\n");
      printf("       3 - float\n");
      printf("       4 - double\n");
      printf("       5 - complex<double>\n");
      printf("Example Input GPU:\n");
      printf("  Histogram : 1000000 1000 1 1000 1 10 1\n");
      printf("  MD Force : 100000 100000 100 1000 20 10 4\n");
      printf("  Matrix Assembly : 100000 1000000 50 1000 20 10 4\n");
      Kokkos::finalize();
      return 0;
    }

    int L    = atoi(argv[1]);
    int N    = atoi(argv[2]);
    int M    = atoi(argv[3]);
    int D    = atoi(argv[4]);
    int K    = atoi(argv[5]);
    int R    = atoi(argv[6]);
    int type = atoi(argv[7]);

    Kokkos::View<int*> offsets("Offsets", L, M);
    Kokkos::Random_XorShift64_Pool<> pool(12371);
    Kokkos::fill_random(offsets, pool, D);
    double time = 0;
    if (type == 1) time = test_atomic<int>(L, N, M, K, R, offsets);
    if (type == 2) time = test_atomic<long>(L, N, M, K, R, offsets);
    if (type == 3) time = test_atomic<float>(L, N, M, K, R, offsets);
    if (type == 4) time = test_atomic<double>(L, N, M, K, R, offsets);
    if (type == 5)
      time = test_atomic<Kokkos::complex<double> >(L, N, M, K, R, offsets);

    double time2 = 1;
    if (type == 1) time2 = test_no_atomic<int>(L, N, M, K, R, offsets);
    if (type == 2) time2 = test_no_atomic<long>(L, N, M, K, R, offsets);
    if (type == 3) time2 = test_no_atomic<float>(L, N, M, K, R, offsets);
    if (type == 4) time2 = test_no_atomic<double>(L, N, M, K, R, offsets);
    if (type == 5)
      time2 = test_no_atomic<Kokkos::complex<double> >(L, N, M, K, R, offsets);

    int size = 0;
    if (type == 1) size = sizeof(int);
    if (type == 2) size = sizeof(long);
    if (type == 3) size = sizeof(float);
    if (type == 4) size = sizeof(double);
    if (type == 5) size = sizeof(Kokkos::complex<double>);

    printf("%i\n", size);
    printf(
        "Time: %s %i %i %i %i %i %i (t_atomic: %e t_nonatomic: %e ratio: %lf "
        ")( GUpdates/s: %lf GB/s: %lf )\n",
        (type == 1)
            ? "int"
            : ((type == 2)
                   ? "long"
                   : ((type == 3) ? "float"
                                  : ((type == 4) ? "double" : "complex"))),
        L, N, M, D, K, R, time, time2, time / time2, 1.e-9 * L * R * M / time,
        1.0 * L * R * M * 2 * size / time / 1024 / 1024 / 1024);
  }
  Kokkos::finalize();
}
