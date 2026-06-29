# Kokkos minimally requires 3.25.2 for cuda language
cmake_minimum_required(VERSION 3.22)

# Projects can safely mix languages - must have C++ support
# Kokkos flags will only apply to C++ files
project(Example CXX)

# Look for an installed Kokkos
find_package(Kokkos)
if(NOT ${Kokkos_ENABLE_MULTIPLE_CMAKE_LANGUAGES})
  message(FATAL_ERROR "This example needs Kokkos_ENABLE_MULTIPLE_CMAKE_LANGUAGES=ON")
endif()

#use Kokkos in a CUDA/HIP language project
if(Kokkos_ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_ARCHITECTURES ${Kokkos_CUDA_ARCHITECTURES})

  add_library(language_lib library.cu)
  target_link_libraries(language_lib PUBLIC Kokkos::kokkos)

elseif(Kokkos_ENABLE_HIP)
  enable_language(HIP)
  set(CMAKE_HIP_ARCHITECTURES ${Kokkos_HIP_ARCHITECTURES})

  add_library(language_lib library.hip)
  target_link_libraries(language_lib PUBLIC Kokkos::kokkos)
endif()

#use Kokkos in a CXX project
add_library(cxx_lib library.cpp)
target_link_libraries(cxx_lib PUBLIC Kokkos::kokkos)

#use both in an executable
add_executable(example cmake_example.cpp)

target_link_libraries(example PUBLIC language_lib cxx_lib)

enable_testing()
add_test(NAME KokkosInTree_Verify COMMAND example 10)
