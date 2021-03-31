#include <Kokkos_Core.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  Kokkos::finalize();
}
