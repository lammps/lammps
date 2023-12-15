# Create tests that contain each header separately. We do not  run these tests
# but we just try to compile them.

# Globbing all the header filenames to test for self-containment and presence of header guards
SET(BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../")
file(GLOB KOKKOS_CORE_HEADERS RELATIVE ${BASE_DIR}/core/src
     ${BASE_DIR}/core/src/*.hpp ${BASE_DIR}/core/src/*.h)
file(GLOB KOKKOS_CONTAINERS_HEADERS RELATIVE ${BASE_DIR}/containers/src
     ${BASE_DIR}/containers/src/*.hpp)
file(GLOB KOKKOS_ALGORITHMS_HEADERS RELATIVE  ${BASE_DIR}/algorithms/src
     ${BASE_DIR}/algorithms/src/*.hpp)

if(NOT Kokkos_ENABLE_DEPRECATED_CODE_4)
  list(REMOVE_ITEM KOKKOS_CONTAINERS_HEADERS "Kokkos_Vector.hpp")
endif()

foreach (_header ${KOKKOS_CORE_HEADERS} ${KOKKOS_CONTAINERS_HEADERS} ${KOKKOS_ALGORITHMS_HEADERS})
  string(REGEX REPLACE "[\./]" "_" header_test_name ${_header})
  set(header_test_name Kokkos_HeaderSelfContained_${header_test_name})
  set_source_files_properties(tstHeader.cpp PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})
  add_executable(${header_test_name} tstHeader.cpp)
  target_link_libraries(${header_test_name} PRIVATE Kokkos::kokkos)
  target_compile_definitions(${header_test_name} PRIVATE KOKKOS_HEADER_TEST_NAME=${_header})
endforeach()
