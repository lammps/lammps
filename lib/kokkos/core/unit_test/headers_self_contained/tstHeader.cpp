#define KOKKOS_HEADER_TEST_STRINGIZE_IMPL(x) #x
#define KOKKOS_HEADER_TEST_STRINGIZE(x) KOKKOS_HEADER_TEST_STRINGIZE_IMPL(x)

#define KOKKOS_HEADER_TO_TEST \
  KOKKOS_HEADER_TEST_STRINGIZE(KOKKOS_HEADER_TEST_NAME)

#define KOKKOS_IMPL_PUBLIC_INCLUDE

// include header twice to see if the include guards are set correctly
#include KOKKOS_HEADER_TO_TEST
#include KOKKOS_HEADER_TO_TEST

#if !defined(KOKKOS_MACROS_HPP)
#error "This header does not include Kokkos_Macros.hpp"
#endif

int main() { return 0; }
