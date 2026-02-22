#These are tribits wrappers only ever called by Kokkos itself

include(CMakeParseArguments)
include(CTest)
include(GNUInstallDirs)

message(STATUS "The project name is: ${PROJECT_NAME}")

if(GTest_FOUND)
  set(KOKKOS_GTEST_LIB GTest::gtest)
  message(STATUS "Using gtest found in ${GTest_DIR}")
else() # fallback to internal gtest
  set(KOKKOS_GTEST_LIB kokkos_gtest)
  message(STATUS "Using internal gtest for testing")
endif()

function(VERIFY_EMPTY CONTEXT)
  if(${ARGN})
    message(FATAL_ERROR "Kokkos does not support all of Tribits. Unhandled arguments in ${CONTEXT}:\n${ARGN}")
  endif()
endfunction()

macro(KOKKOS_PROCESS_SUBPACKAGES)
  add_subdirectory(core)
  add_subdirectory(containers)
  add_subdirectory(algorithms)
  add_subdirectory(simd)
  add_subdirectory(example)
  add_subdirectory(benchmarks)
endmacro()

macro(KOKKOS_INTERNAL_ADD_LIBRARY_INSTALL LIBRARY_NAME)
  kokkos_lib_type(${LIBRARY_NAME} INCTYPE)
  target_include_directories(${LIBRARY_NAME} ${INCTYPE} $<INSTALL_INTERFACE:${KOKKOS_HEADER_DIR}>)

  if(Kokkos_ENABLE_EXPERIMENTAL_CXX20_MODULES)
    install(
      TARGETS ${LIBRARY_NAME}
      EXPORT KokkosTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} FILE_SET ${LIBRARY_NAME}_file_set
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/src
    )
  else()
    install(
      TARGETS ${LIBRARY_NAME}
      EXPORT KokkosTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  endif()

  verify_empty(KOKKOS_ADD_LIBRARY ${PARSE_UNPARSED_ARGUMENTS})
endmacro()

function(KOKKOS_ADD_EXECUTABLE ROOT_NAME)
  cmake_parse_arguments(PARSE "TESTONLY" "" "SOURCES;TESTONLYLIBS" ${ARGN})

  set_source_files_properties(${PARSE_SOURCES} PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})

  set(EXE_NAME ${PACKAGE_NAME}_${ROOT_NAME})
  add_executable(${EXE_NAME} ${PARSE_SOURCES})
  if(PARSE_TESTONLYLIBS)
    target_link_libraries(${EXE_NAME} PRIVATE ${PARSE_TESTONLYLIBS})
  endif()
  verify_empty(KOKKOS_ADD_EXECUTABLE ${PARSE_UNPARSED_ARGUMENTS})
  #All executables must link to all the kokkos targets
  #This is just private linkage because exe is final
  target_link_libraries(${EXE_NAME} PRIVATE Kokkos::kokkos)
endfunction()

function(KOKKOS_ADD_EXECUTABLE_AND_TEST ROOT_NAME)
  cmake_parse_arguments(PARSE "" "" "SOURCES;CATEGORIES;ARGS" ${ARGN})
  verify_empty(KOKKOS_ADD_EXECUTABLE_AND_TEST ${PARSE_UNPARSED_ARGUMENTS})

  kokkos_add_test_executable(${ROOT_NAME} SOURCES ${PARSE_SOURCES})
  if(PARSE_ARGS)
    set(TEST_NUMBER 0)
    foreach(ARG_STR ${PARSE_ARGS})
      # This is passed as a single string blob to match TriBITS behavior
      # We need this to be turned into a list
      string(REPLACE " " ";" ARG_STR_LIST ${ARG_STR})
      list(APPEND TEST_NAME "${ROOT_NAME}${TEST_NUMBER}")
      math(EXPR TEST_NUMBER "${TEST_NUMBER} + 1")
      kokkos_add_test(
        NAME
        ${TEST_NAME}
        EXE
        ${ROOT_NAME}
        FAIL_REGULAR_EXPRESSION
        "  FAILED  "
        ARGS
        ${ARG_STR_LIST}
      )
    endforeach()
  else()
    kokkos_add_test(NAME ${ROOT_NAME} EXE ${ROOT_NAME} FAIL_REGULAR_EXPRESSION "  FAILED  ")
  endif()
  # We noticed problems with -fvisibility=hidden for inline static variables
  # if Kokkos was built as shared library.
  if(BUILD_SHARED_LIBS AND NOT ${TEST_NAME}_DISABLE)
    set_property(TARGET ${EXE_NAME} PROPERTY VISIBILITY_INLINES_HIDDEN ON)
    set_property(TARGET ${EXE_NAME} PROPERTY CXX_VISIBILITY_PRESET hidden)
  endif()
  if(NOT
     (Kokkos_INSTALL_TESTING
      OR Kokkos_ENABLE_EXPERIMENTAL_CXX20_MODULES
      OR Kokkos_ENABLE_SYCL
      OR Kokkos_ENABLE_HPX
      OR (KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA" AND KOKKOS_CXX_HOST_COMPILER_ID STREQUAL "MSVC"))
  )
    if(MSVC)
      target_compile_options(${PACKAGE_NAME}_${ROOT_NAME} PRIVATE "/GR-")
    else()
      target_compile_options(${PACKAGE_NAME}_${ROOT_NAME} PRIVATE "-fno-rtti")
    endif()
  endif()
endfunction()

function(KOKKOS_SET_EXE_PROPERTY ROOT_NAME)
  set(TARGET_NAME ${PACKAGE_NAME}_${ROOT_NAME})
  if(NOT TARGET ${TARGET_NAME})
    message(SEND_ERROR "No target ${TARGET_NAME} exists - cannot set target properties")
  endif()
  set_property(TARGET ${TARGET_NAME} PROPERTY ${ARGN})
endfunction()

macro(KOKKOS_SETUP_BUILD_ENVIRONMENT)
  # This is needed for both regular build and install tests
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_compiler_id.cmake)
  #set an internal option, if not already set
  set(Kokkos_INSTALL_TESTING OFF CACHE INTERNAL "Whether to build tests and examples against installation")
  if(Kokkos_INSTALL_TESTING)
    set(KOKKOS_ENABLE_TESTS ON)
    set(KOKKOS_ENABLE_BENCHMARKS ON)
    set(KOKKOS_ENABLE_EXAMPLES ON)
    # This looks a little weird, but what we are doing
    # is to NOT build Kokkos but instead look for an
    # installed Kokkos - then build examples and tests
    # against that installed Kokkos
    find_package(Kokkos REQUIRED)
    # Just grab the configuration from the installation
    foreach(DEV ${Kokkos_DEVICES})
      set(KOKKOS_ENABLE_${DEV} ON)
    endforeach()
    foreach(OPT ${Kokkos_OPTIONS})
      set(KOKKOS_ENABLE_${OPT} ON)
    endforeach()
    foreach(TPL ${Kokkos_TPLS})
      set(KOKKOS_ENABLE_${TPL} ON)
    endforeach()
    foreach(ARCH ${Kokkos_ARCH})
      set(KOKKOS_ARCH_${ARCH} ON)
    endforeach()
  else()
    include(${KOKKOS_SRC_PATH}/cmake/kokkos_enable_devices.cmake)
    include(${KOKKOS_SRC_PATH}/cmake/kokkos_enable_options.cmake)
    include(${KOKKOS_SRC_PATH}/cmake/kokkos_test_cxx_std.cmake)
    include(${KOKKOS_SRC_PATH}/cmake/kokkos_arch.cmake)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Kokkos_SOURCE_DIR}/cmake/Modules/")
    include(${KOKKOS_SRC_PATH}/cmake/kokkos_tpls.cmake)
  endif()
endmacro()

macro(KOKKOS_ADD_TEST_EXECUTABLE ROOT_NAME)
  cmake_parse_arguments(PARSE "" "" "SOURCES" ${ARGN})
  # Don't do anything if the user disabled the test
  if(NOT ${PACKAGE_NAME}_${ROOT_NAME}_DISABLE)
    kokkos_add_executable(
      ${ROOT_NAME} SOURCES ${PARSE_SOURCES} ${PARSE_UNPARSED_ARGUMENTS} TESTONLYLIBS ${KOKKOS_GTEST_LIB}
    )
    set(EXE_NAME ${PACKAGE_NAME}_${ROOT_NAME})
  endif()
endmacro()

## KOKKOS_CONFIGURE_CORE  Configure/Generate header files for core content based
##                        on enabled backends.
##                        KOKKOS_FWD is the forward declare set
##                        KOKKOS_SETUP  is included in Kokkos_Macros.hpp and include prefix includes/defines
##                        KOKKOS_DECLARE is the declaration set
##                        KOKKOS_POST_INCLUDE is included at the end of Kokkos_Core.hpp
macro(KOKKOS_CONFIGURE_CORE)
  message(STATUS "Kokkos Backends: ${KOKKOS_ENABLED_DEVICES}")
  kokkos_config_header(
    KokkosCore_Config_HeaderSet.in KokkosCore_Config_FwdBackend.hpp "KOKKOS_FWD" "fwd/Kokkos_Fwd"
    "${KOKKOS_ENABLED_DEVICES}"
  )
  kokkos_config_header(
    KokkosCore_Config_HeaderSet.in KokkosCore_Config_SetupBackend.hpp "KOKKOS_SETUP" "setup/Kokkos_Setup"
    "${DEVICE_SETUP_LIST}"
  )
  kokkos_config_header(
    KokkosCore_Config_HeaderSet.in KokkosCore_Config_DeclareBackend.hpp "KOKKOS_DECLARE" "decl/Kokkos_Declare"
    "${KOKKOS_ENABLED_DEVICES}"
  )
  if(BUILD_SHARED_LIBS)
    set(KOKKOS_IMPL_BUILD_SHARED_LIBS ON)
  endif()
  configure_file(cmake/KokkosCore_config.h.in KokkosCore_config.h @ONLY)
endmacro()

## KOKKOS_INSTALL_ADDITIONAL_FILES - instruct cmake to install files in target destination.
##                        Includes generated header files, scripts such as nvcc_wrapper and hpcbind,
##                        as well as other files provided through plugins.
macro(KOKKOS_INSTALL_ADDITIONAL_FILES)

  # kokkos_launch_compiler is used by Kokkos to prefix compiler commands so that they forward to original kokkos compiler
  # if nvcc_wrapper was not used as CMAKE_CXX_COMPILER, configure the original compiler into kokkos_launch_compiler
  if(NOT "${CMAKE_CXX_COMPILER}" MATCHES "nvcc_wrapper")
    set(NVCC_WRAPPER_DEFAULT_COMPILER "${CMAKE_CXX_COMPILER}")
  else()
    if(NOT "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}" STREQUAL "")
      set(NVCC_WRAPPER_DEFAULT_COMPILER "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}")
    endif()
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/bin/kokkos_launch_compiler ${PROJECT_BINARY_DIR}/temp/kokkos_launch_compiler @ONLY
  )

  install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/bin/nvcc_wrapper" "${CMAKE_CURRENT_SOURCE_DIR}/bin/hpcbind"
                   "${PROJECT_BINARY_DIR}/temp/kokkos_launch_compiler" DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_config.h"
          "${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_Config_FwdBackend.hpp"
          "${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_Config_SetupBackend.hpp"
          "${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_Config_DeclareBackend.hpp"
    DESTINATION ${KOKKOS_HEADER_DIR}
  )
endmacro()

function(KOKKOS_SET_LIBRARY_PROPERTIES LIBRARY_NAME)
  cmake_parse_arguments(PARSE "PLAIN_STYLE" "" "" ${ARGN})

  #allow multiple languages to be used downstream
  set(Kokkos_LANGUAGES "${KOKKOS_COMPILE_LANGUAGE}")
  if(Kokkos_ENABLE_MULTIPLE_CMAKE_LANGUAGES)
    if(Kokkos_ENABLE_HIP)
      set(Kokkos_LANGUAGES "HIP,CXX")
    elseif(Kokkos_ENABLE_CUDA)
      set(Kokkos_LANGUAGES "CUDA,CXX")
    endif()
  endif()

  if(NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
    #I can use link options
    #check for CXX linkage using the simple 3.18 way
    target_link_options(${LIBRARY_NAME} PUBLIC $<$<LINK_LANGUAGE:${Kokkos_LANGUAGES}>:${KOKKOS_LINK_OPTIONS}>)
  else()
    #I can use link options
    #just assume CXX linkage
    target_link_options(${LIBRARY_NAME} PUBLIC ${KOKKOS_LINK_OPTIONS})
  endif()

  #exclude case of compiler_launcher. The launcher forwards to nvcc_wrapper and shadow the CXX compiler that CMake sees (compiler_launcher changes the compiler).
  #The CXX compiler CMake will invoke for the check is not able to consume the cuda flags if it is not nvcc_wrapper or clang+cuda.
  #FIXME_NVHPC nvc++ is failing the check spuriously with various version numbers.
  if(NOT (KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
     AND (NOT (KOKKOS_ENABLE_CUDA) OR ("${CMAKE_CXX_COMPILER}" MATCHES "nvcc_wrapper") OR (${KOKKOS_CXX_COMPILER_ID}
                                                                                           STREQUAL Clang))
  )
    kokkos_check_flags(LINKER LANGUAGE ${KOKKOS_COMPILE_LANGUAGE} FLAGS ${KOKKOS_LINK_OPTIONS})
  endif()

  list(APPEND ALL_KOKKOS_COMPILER_FLAGS ${KOKKOS_COMPILE_OPTIONS})
  target_compile_options(${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${Kokkos_LANGUAGES}>:${KOKKOS_COMPILE_OPTIONS}>)

  target_compile_definitions(
    ${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${Kokkos_LANGUAGES}>:${KOKKOS_COMPILE_DEFINITIONS}>
  )

  target_link_libraries(${LIBRARY_NAME} PUBLIC ${KOKKOS_LINK_LIBRARIES})

  if(KOKKOS_ENABLE_CUDA)
    target_compile_options(${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${Kokkos_LANGUAGES}>:${KOKKOS_CUDA_OPTIONS}>)
    set(NODEDUP_CUDAFE_OPTIONS)
    foreach(OPT ${KOKKOS_CUDAFE_OPTIONS})
      list(APPEND NODEDUP_CUDAFE_OPTIONS -Xcudafe ${OPT})
    endforeach()
    target_compile_options(
      ${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${KOKKOS_COMPILE_LANGUAGE}>:${NODEDUP_CUDAFE_OPTIONS}>
    )

    list(APPEND ALL_KOKKOS_COMPILER_FLAGS ${KOKKOS_CUDA_OPTIONS})
    list(APPEND ALL_KOKKOS_COMPILER_FLAGS ${NODEDUP_CUDAFE_OPTIONS})
  endif()

  if(KOKKOS_ENABLE_HIP)
    target_compile_options(${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${Kokkos_LANGUAGES}>:${KOKKOS_AMDGPU_OPTIONS}>)
    list(APPEND ALL_KOKKOS_COMPILER_FLAGS ${KOKKOS_AMDGPU_OPTIONS})
  endif()

  list(LENGTH KOKKOS_XCOMPILER_OPTIONS XOPT_LENGTH)
  if(XOPT_LENGTH GREATER 1)
    message(
      FATAL_ERROR
        "CMake deduplication does not allow multiple -Xcompiler flags (${KOKKOS_XCOMPILER_OPTIONS}): will require Kokkos to upgrade to minimum 3.12"
    )
  endif()
  if(KOKKOS_XCOMPILER_OPTIONS)
    set(NODEDUP_XCOMPILER_OPTIONS)
    foreach(OPT ${KOKKOS_XCOMPILER_OPTIONS})
      #I have to do this for now because we can't guarantee 3.12 support
      #I really should do this with the shell option
      list(APPEND NODEDUP_XCOMPILER_OPTIONS -Xcompiler)
      list(APPEND NODEDUP_XCOMPILER_OPTIONS ${OPT})
    endforeach()
    target_compile_options(
      ${LIBRARY_NAME} PUBLIC $<$<COMPILE_LANGUAGE:${KOKKOS_COMPILE_LANGUAGE}>:${NODEDUP_XCOMPILER_OPTIONS}>
    )
    list(APPEND ALL_KOKKOS_COMPILER_FLAGS ${NODEDUP_XCOMPILER_OPTIONS})
  endif()

  #exclude case of compiler_launcher. The launcher forwards to nvcc_wrapper and shadow the CXX compiler that CMake sees (compiler_launcher changes the compiler).
  #The CXX compiler CMake will invoke for the check is not able to consume the cuda flags if it is not nvcc_wrapper or clang+cuda.
  #FIXME_NVHPC nvc++ is failing the check spuriously with various version numbers.
  #FIXME CLANG+RDC with std=c++20 and CMake 3.22+ the runtime is not linked in the CheckCompilerFlag leading to false positives
  if(NOT (KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
     AND (NOT (KOKKOS_ENABLE_CUDA) OR ("${CMAKE_CXX_COMPILER}" MATCHES "nvcc_wrapper")
          OR (${KOKKOS_CXX_COMPILER_ID} STREQUAL Clang AND NOT KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE))
  )
    kokkos_check_flags(
      COMPILER
      LANGUAGE
      ${KOKKOS_COMPILE_LANGUAGE}
      FLAGS
      ${ALL_KOKKOS_COMPILER_FLAGS}
      LINKER_FLAGS
      ${KOKKOS_LINK_OPTIONS}
    )
  endif()

  if(KOKKOS_CXX_STANDARD_FEATURE)
    #GREAT! I can do this the right way
    target_compile_features(${LIBRARY_NAME} PUBLIC ${KOKKOS_CXX_STANDARD_FEATURE})
    if(NOT KOKKOS_USE_CXX_EXTENSIONS)
      set_target_properties(${LIBRARY_NAME} PROPERTIES CXX_EXTENSIONS OFF)
    endif()
  else()
    #OH, well, no choice but the wrong way
    target_compile_options(${LIBRARY_NAME} PUBLIC ${KOKKOS_CXX_STANDARD_FLAG})
  endif()
endfunction()

function(KOKKOS_ADD_LIBRARY LIBRARY_NAME)
  cmake_parse_arguments(PARSE "ADD_BUILD_OPTIONS;STATIC;SHARED" "" "HEADERS;SOURCES;MODULE_INTERFACE" ${ARGN})

  if(PARSE_HEADERS)
    list(REMOVE_DUPLICATES PARSE_HEADERS)
  endif()
  if(PARSE_SOURCES)
    list(REMOVE_DUPLICATES PARSE_SOURCES)
  endif()
  if(PARSE_MODULE_INTERFACE)
    list(REMOVE_DUPLICATES PARSE_MODULE_INTERFACE)
  endif()
  foreach(source ${PARSE_SOURCES})
    set_source_files_properties(${source} PROPERTIES LANGUAGE ${KOKKOS_COMPILE_LANGUAGE})
  endforeach()

  if(PARSE_STATIC)
    set(LINK_TYPE STATIC)
  endif()

  if(PARSE_SHARED)
    set(LINK_TYPE SHARED)
  endif()

  # MSVC and other platforms want to have the headers included as source files for better dependency detection
  add_library(${LIBRARY_NAME} ${LINK_TYPE} ${PARSE_HEADERS} ${PARSE_SOURCES})
  if(PARSE_STATIC)
    # We want position independent code to be able to link static libraries into shared ones.
    set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  endif()

  if(Kokkos_ENABLE_EXPERIMENTAL_CXX20_MODULES)
    if(PARSE_MODULE_INTERFACE)
      target_sources(
        ${LIBRARY_NAME} PUBLIC FILE_SET ${LIBRARY_NAME}_file_set TYPE CXX_MODULES FILES ${PARSE_MODULE_INTERFACE}
      )
    endif()
  endif()

  if(PARSE_SHARED OR BUILD_SHARED_LIBS)
    set_target_properties(
      ${LIBRARY_NAME} PROPERTIES VERSION ${Kokkos_VERSION} SOVERSION ${Kokkos_VERSION_MAJOR}.${Kokkos_VERSION_MINOR}
                                 WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
    target_compile_definitions(${LIBRARY_NAME} PRIVATE KOKKOS_IMPL_EXPORT_SYMBOLS)
  endif()

  kokkos_internal_add_library_install(${LIBRARY_NAME})

  #In case we are building in-tree, add an alias name
  #that matches the install Kokkos:: name
  add_library(Kokkos::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

  if(PARSE_ADD_BUILD_OPTIONS)
    kokkos_set_library_properties(${LIBRARY_NAME})
  endif()
endfunction()

function(KOKKOS_LIB_INCLUDE_DIRECTORIES TARGET)
  kokkos_lib_type(${TARGET} INCTYPE)
  foreach(DIR ${ARGN})
    target_include_directories(${TARGET} ${INCTYPE} $<BUILD_INTERFACE:${DIR}>)
  endforeach()
endfunction()

function(KOKKOS_LIB_COMPILE_OPTIONS TARGET)
  kokkos_lib_type(${TARGET} INCTYPE)
  target_compile_options(${${PROJECT_NAME}_LIBRARY_NAME_PREFIX}${TARGET} ${INCTYPE} ${ARGN})
endfunction()

macro(KOKKOS_ADD_TEST_DIRECTORIES)
  if(KOKKOS_ENABLE_TESTS)
    foreach(TEST_DIR ${ARGN})
      add_subdirectory(${TEST_DIR})
    endforeach()
  endif()
endmacro()

macro(KOKKOS_ADD_EXAMPLE_DIRECTORIES)
  if(KOKKOS_ENABLE_EXAMPLES)
    foreach(EXAMPLE_DIR ${ARGN})
      add_subdirectory(${EXAMPLE_DIR})
    endforeach()
  endif()
endmacro()

macro(KOKKOS_ADD_BENCHMARK_DIRECTORIES)
  if(KOKKOS_ENABLE_BENCHMARKS)
    foreach(BENCHMARK_DIR ${ARGN})
      add_subdirectory(${BENCHMARK_DIR})
    endforeach()
  endif()
endmacro()
