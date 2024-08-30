# gcc bundles libquadmath and doesn't need any extra link or include directories
# (which would not be contained in CMake's search paths anyway).
# Hence, try if the compiler supports libquadmath natively first before doing
# the standard package search.
SET(CMAKE_REQUIRED_LIBRARIES "quadmath")
INCLUDE(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
  #include <quadmath.h>
  int main(void){
      __float128 foo = ::sqrtq(123.456);
      return foo;
  }"
  KOKKOS_QUADMATH_COMPILER_SUPPORT)
IF (KOKKOS_QUADMATH_COMPILER_SUPPORT)
  KOKKOS_CREATE_IMPORTED_TPL(LIBQUADMATH INTERFACE LINK_LIBRARIES quadmath)
ELSE()
  KOKKOS_FIND_IMPORTED(LIBQUADMATH HEADER quadmath.h LIBRARY quadmath)
ENDIF()
