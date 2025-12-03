# PACE library support for ML-PACE package

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2024.9.11.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "8409ed1a7dceb15de0ff2637a06369f6" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

# LOCAL_ML-PACE points to top-level dir with local lammps-user-pace repo,
# to make it easier to check local build without going through the public github releases
if(LOCAL_ML-PACE)
  message(STATUS "Using LOCAL ML-PACE ${LOCAL_ML-PACE}")
  set(lib-pace "${LOCAL_ML-PACE}")
else()
  # download library sources to build folder
  if(EXISTS ${CMAKE_BINARY_DIR}/libpace.tar.gz)
    file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
  endif()
  if(NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}")
    message(STATUS "Downloading ${PACELIB_URL}")
    file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
    file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
    if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}"))
      message(WARNING "Download from primary URL ${PACELIB_URL} failed\nTrying fallback URL ${PACELIB_FALLBACK}")
      file(DOWNLOAD ${PACELIB_FALLBACK} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH MD5=${PACELIB_MD5} SHOW_PROGRESS)
    endif()
  else()
    message(STATUS "Using already downloaded archive ${CMAKE_BINARY_DIR}/libpace.tar.gz")
  endif()


  # uncompress downloaded sources
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
    COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)
endif()

add_subdirectory(${lib-pace} build-pace)
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

# GRACE/TensorFlow is compiled by default
if(NOT DEFINED NO_GRACE_TF)
  # We will compile with TF support

  # Check, if TF_LIB_FILE is provided  
  # TODO: add to doc TF_LIB_FILE and   
  if(TF_LIB_FILE) 
    message("User-defined TF_LIB_FILE is provided: ${TF_LIB_FILE}")
  else()
    # 1) try to find TensorFlow library  from Python installation (for older versions of TF)
    
    # get default python
    if(NOT PACE_PYTHON_EXEC)
      find_package(Python COMPONENTS Interpreter QUIET)
      set(PACE_PYTHON_EXEC ${Python_EXECUTABLE})
    endif()
    message("Python interpreter found: ${PACE_PYTHON_EXEC}")
    execute_process(
      COMMAND ${PACE_PYTHON_EXEC} -c "import os;import pkgutil;package = pkgutil.get_loader('tensorflow');print(os.path.dirname(package.get_filename()))"
      OUTPUT_VARIABLE TF_DISCOVER
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # message("TF_DISCOVER=${TF_DISCOVER}")
    string(STRIP "${TF_DISCOVER}" TF_DISCOVER)
    set(TF_PATH ${TF_DISCOVER})

    if(APPLE)
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.2.dylib")
    elseif(WIN32)
      set(TF_LIB_FILE "${TF_PATH}/tensorflow.dll")
    else()
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.so.2")
    endif()
    # setup include path
    set(TF_INCLUDE_PATH "${TF_PATH}/include")
  

    # 2) If not found, download it 
    if(NOT EXISTS ${TF_LIB_FILE})
      # Define URLs for TensorFlow C++ library for different platforms
      set(TF_URL_WINDOWS "https://storage.googleapis.com/tensorflow/versions/2.18.1/libtensorflow-cpu-windows-x86_64.zip")
      set(TF_URL_LINUX   "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-gpu-linux-x86_64.tar.gz")
      set(TF_URL_MACOS   "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-cpu-darwin-arm64.tar.gz")
      set(TF_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/tensorflow-library-download")
      
      message(STATUS "TensorFlow library not found via Python discovery. Attempting to download.")

      if(WIN32)
        set(TF_URL ${TF_URL_WINDOWS})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.zip")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xf)
      elseif(APPLE)
        set(TF_URL ${TF_URL_MACOS})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xzf)
      else() # linux
        set(TF_URL ${TF_URL_LINUX})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xzf)
      endif()

      if(NOT EXISTS ${TF_ARCHIVE})
        message(STATUS "Downloading TensorFlow C library from ${TF_URL}")
        file(DOWNLOAD ${TF_URL} ${TF_ARCHIVE} SHOW_PROGRESS STATUS DL_STATUS)
        if(NOT DL_STATUS EQUAL 0)
          message(FATAL_ERROR "Failed to download TensorFlow from ${TF_URL}")
        endif()
      else()
        message(STATUS "Using already downloaded archive ${TF_ARCHIVE}")
      endif()

      message(STATUS "Clean folder for archive ${TF_DOWNLOAD_DIR}...")
      # file(REMOVE_RECURSE ${TF_DOWNLOAD_DIR})
      message(STATUS "Extracting TensorFlow library archive ${TF_ARCHIVE} to current working dir ${CMAKE_BINARY_DIR}...")
      execute_process(
        COMMAND ${EXTRACT_COMMAND} ${TF_ARCHIVE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
      # message("Rename ${TF_DOWNLOAD_DIR} to ${CMAKE_BINARY_DIR}/lib")
      # file(RENAME ${TF_DOWNLOAD_DIR} ${CMAKE_BINARY_DIR}/lib)
      set(TF_PATH ${CMAKE_BINARY_DIR})

      # setup library path
      if(WIN32)
        set(TF_LIB_FILE "${TF_PATH}/lib/tensorflow.dll") # Path inside downloaded archive
        string(REPLACE ".dll" ".lib" TF_IMPORTS_LIB_FILE "${TF_LIB_FILE}")
      elseif(APPLE)
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.2.dylib")
      else() # linux
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.so.2")      
      endif()

      # setup include path
      set(TF_INCLUDE_PATH "${TF_PATH}/include")
    endif()
  endif()


  # 3) Finally, import library or fail
  if(EXISTS ${TF_LIB_FILE})
    message("-- TensorFlow library is FOUND at ${TF_LIB_FILE}")
    add_library(tensorflow SHARED IMPORTED)
    if(WIN32)
      set_target_properties(tensorflow PROPERTIES
              IMPORTED_LOCATION "${TF_LIB_FILE}"
              IMPORTED_IMPLIB "${TF_IMPORTS_LIB_FILE}"
              INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    elseif(APPLE)
      set_target_properties(tensorflow PROPERTIES
              IMPORTED_IMPLIB "${TF_PATH}/lib/tensorflow.lib"
              INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    else()
      set_target_properties(tensorflow PROPERTIES
              IMPORTED_LOCATION ${TF_LIB_FILE}
              INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    endif()
    
    ###############################
    # download cppflow
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/cppflow-2.0.0)
        if(NOT EXISTS ${CMAKE_BINARY_DIR}/libcppflow.tar.gz)
            set(CPPFLOW_URL "https://github.com/serizba/cppflow/archive/refs/tags/v2.0.0.tar.gz" CACHE STRING "URL for cppflow")
            message(STATUS "Downloading ${CPPFLOW_URL}")
            file(DOWNLOAD ${CPPFLOW_URL} ${CMAKE_BINARY_DIR}/libcppflow.tar.gz STATUS DL_CPPFLOW_STATUS)
            # uncompress downloaded sources
            execute_process(
                    COMMAND ${CMAKE_COMMAND} -E remove_directory cppflow-*
                    COMMAND ${CMAKE_COMMAND} -E tar xzf libcppflow.tar.gz
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
        else()
          message(STATUS "Using already downloaded archive  ${CMAKE_BINARY_DIR}/libcppflow.tar.gz")
        endif()
    else()
      message(STATUS "Using already existing CPPFLOW ${CMAKE_BINARY_DIR}/cppflow-2.0.0")
    endif()

    set(cppflow_path "${CMAKE_BINARY_DIR}/cppflow-2.0.0")
    add_library(cppflow INTERFACE)
    target_include_directories(cppflow
            INTERFACE
            ${tensorflow_INCLUDE_DIRS}
            $<BUILD_INTERFACE:${cppflow_path}/include>
    )

    target_compile_features(cppflow INTERFACE cxx_std_17)
    target_link_libraries(cppflow INTERFACE
            ${tensorflow_LIBRARIES}
    )

    set(PACE_TP ON)
    find_package(OpenMP)
  else()
    message("-- TensorFlow library is NEITHER found at ${TF_LIB_FILE} NOR downloaded/extracted")
  endif()
else()
  message("-- NO GRACE/TensorFlow will be compiled (because flag NO_GRACE_TF is set)")
  add_definitions(-DNO_GRACE_TF=1)
endif() # if(NOT DEFINED NO_GRACE_TF)

if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)

  if(DEFINED PACE_TP)
    add_definitions(-DPACE_TP)
    target_link_libraries(lammps PRIVATE tensorflow)
    target_link_libraries(lammps PRIVATE cppflow)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(lammps PUBLIC OpenMP::OpenMP_CXX)
    endif()
  endif()

  if(WIN32)
    if(BUILD_SHARED_LIBS)
      # Building lammps AND pace/yaml-cpp as shared libs (.dll)
      # Tell lammps (pair_grace.cpp) to IMPORT symbols from the DLL.
      message(STATUS "ML-PACE: Configuring 'lammps' for shared library import (YAML_CPP_DLL)")
      target_compile_definitions(lammps PRIVATE YAML_CPP_DLL)
    else()
      # Building lammps AND pace/yaml-cpp as static libs (.lib)
      # Tell lammps (pair_grace.cpp) it's a STATIC lib.
      message(STATUS "ML-PACE: Configuring 'lammps' for static library link (YAML_CPP_STATIC_DEFINE)")
      target_compile_definitions(lammps PRIVATE YAML_CPP_STATIC_DEFINE)
    endif()
  endif()

endif()
