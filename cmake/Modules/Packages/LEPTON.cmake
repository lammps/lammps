set(LEPTON_SOURCE_DIR ${LAMMPS_LIB_SOURCE_DIR}/lepton)

file(GLOB LEPTON_SOURCES ${LEPTON_SOURCE_DIR}/src/[^.]*.cpp)

if((CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64") OR
   (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "amd64") OR
   (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "AMD64"))
   option(LEPTON_ENABLE_JIT "Enable Just-In-Time compiler for Lepton" ON)
else()
   option(LEPTON_ENABLE_JIT "Enable Just-In-Time compiler for Lepton" OFF)
endif()

if(LEPTON_ENABLE_JIT)
  file(GLOB ASMJIT_SOURCES ${LEPTON_SOURCE_DIR}/asmjit/*/[^.]*.cpp)
endif()

add_library(lmplepton STATIC ${LEPTON_SOURCES} ${ASMJIT_SOURCES})
set_target_properties(lmplepton PROPERTIES OUTPUT_NAME lammps_lmplepton${LAMMPS_MACHINE})
target_compile_definitions(lmplepton PUBLIC -DLEPTON_BUILDING_STATIC_LIBRARY=1)
target_include_directories(lmplepton PUBLIC ${LEPTON_SOURCE_DIR}/include)

if(LEPTON_ENABLE_JIT)
  target_compile_definitions(lmplepton PUBLIC "LEPTON_USE_JIT=1;ASMJIT_BUILD_X86=1;ASMJIT_EMBED=1;ASMJIT_BUILD_RELEASE=1")
  target_include_directories(lmplepton PUBLIC ${LEPTON_SOURCE_DIR})
endif()

target_link_libraries(lammps PRIVATE lmplepton)
