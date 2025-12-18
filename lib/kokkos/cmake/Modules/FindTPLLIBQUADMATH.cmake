# gcc bundles libquadmath and doesn't need any extra link or include directories
# (which would not be contained in CMake's search paths anyway).
# Hence, try if the compiler supports libquadmath natively first before doing
# the standard package search.
set(CMAKE_REQUIRED_LIBRARIES "quadmath")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "
  #include <quadmath.h>
  int main(void){
      __float128 foo = ::sqrtq(123.456);
      return foo;
  }"
  KOKKOS_QUADMATH_COMPILER_SUPPORT
)
unset(CMAKE_REQUIRED_LIBRARIES)

if(KOKKOS_QUADMATH_COMPILER_SUPPORT)
  kokkos_create_imported_tpl(LIBQUADMATH INTERFACE LINK_LIBRARIES quadmath)
else()
  kokkos_find_imported(LIBQUADMATH HEADER quadmath.h LIBRARY quadmath)
endif()
