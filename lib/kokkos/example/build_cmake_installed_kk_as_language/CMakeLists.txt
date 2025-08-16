# Kokkos minimally requires 3.16 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.19)

# Projects can safely mix languages - must have C++ support
# Kokkos flags will only apply to C++ files
project(Example CXX Fortran)

# Look for an installed Kokkos
find_package(Kokkos)
enable_language(${Kokkos_COMPILE_LANGUAGE})

#the language has to be set on the files
set_source_files_properties(cmake_example.cpp PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE})
add_executable(example cmake_example.cpp bar.cpp foo.f)

# The architectures and std have to be set on the target
set_property(
  TARGET example PROPERTY ${Kokkos_COMPILE_LANGUAGE}_ARCHITECTURES ${Kokkos_${Kokkos_COMPILE_LANGUAGE}_ARCHITECTURES}
)
set_property(TARGET example PROPERTY ${Kokkos_COMPILE_LANGUAGE}_STANDARD ${Kokkos_${Kokkos_COMPILE_LANGUAGE}_STANDARD})
target_link_libraries(example Kokkos::kokkos)

enable_testing()
add_test(NAME KokkosInTree_Verify COMMAND example 10)
