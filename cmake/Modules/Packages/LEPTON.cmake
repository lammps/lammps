# avoid including this file twice
if(LEPTON_SOURCE_DIR)
   return()
endif()
set(LEPTON_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/lepton)

file(GLOB LEPTON_SOURCES CONFIGURE_DEPENDS ${LEPTON_SOURCE_DIR}/src/[^.]*.cpp)

if((CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "amd64") OR
   (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "AMD64") OR
   (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
   option(LEPTON_ENABLE_JIT "Enable Just-In-Time compiler for Lepton" ON)
else()
   option(LEPTON_ENABLE_JIT "Enable Just-In-Time compiler for Lepton" OFF)
endif()

if(LEPTON_ENABLE_JIT)
  file(GLOB ASMJIT_SOURCES CONFIGURE_DEPENDS ${LEPTON_SOURCE_DIR}/asmjit/*/[^.]*.cpp)
endif()

add_library(lepton STATIC ${LEPTON_SOURCES} ${ASMJIT_SOURCES})
set_target_properties(lepton PROPERTIES OUTPUT_NAME lammps_lepton${LAMMPS_MACHINE})
target_compile_definitions(lepton PUBLIC LEPTON_BUILDING_STATIC_LIBRARY=1)
target_include_directories(lepton PUBLIC ${LEPTON_SOURCE_DIR}/include)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  find_library(LIB_RT rt QUIET)
  target_link_libraries(lepton PUBLIC ${LIB_RT})
endif()

if(LEPTON_ENABLE_JIT)
  target_compile_definitions(lepton PUBLIC "LEPTON_USE_JIT=1;ASMJIT_BUILD_X86=1;ASMJIT_STATIC=1;ASMJIT_BUILD_RELEASE=1")
  target_include_directories(lepton PUBLIC ${LEPTON_SOURCE_DIR})
endif()

target_link_libraries(lammps PRIVATE lepton)
