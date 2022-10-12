#include <Kokkos_Core.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    int N = (argc > 1) ? std::stoi(argv[1]) : 10000;
    int M = (argc > 2) ? std::stoi(argv[2]) : 10000;
    int R = (argc > 3) ? std::stoi(argv[3]) : 10;

    printf("Called with: %i %i %i\n", N, M, R);
  }
  Kokkos::finalize();
}
