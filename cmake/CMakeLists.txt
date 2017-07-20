########################################
# CMake build system
# This file is part of LAMMPS
# Created by Christoph Junghans and Richard Berger
cmake_minimum_required(VERSION 3.1)

project(lammps)
set(SOVERSION 0)
set(LAMMPS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../src)
set(LAMMPS_LIB_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../lib)
set(LAMMPS_LIB_BINARY_DIR ${CMAKE_BINARY_DIR}/lib)

#To not conflict with old Makefile build system, we build everything here
file(GLOB LIB_SOURCES ${LAMMPS_SOURCE_DIR}/*.cpp)
file(GLOB LMP_SOURCES ${LAMMPS_SOURCE_DIR}/main.cpp)
list(REMOVE_ITEM LIB_SOURCES ${LMP_SOURCES})

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/Modules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_C_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_C_FLAGS)

foreach(STYLE_FILE style_angle.h style_atom.h style_body.h style_bond.h style_command.h style_compute.h style_dihedral.h style_dump.h
             style_fix.h style_improper.h style_integrate.h style_kspace.h style_minimize.h style_nbin.h style_npair.h style_nstencil.h
             style_ntopo.h style_pair.h style_reader.h style_region.h)
  if(EXISTS ${LAMMPS_SOURCE_DIR}/${STYLE_FILE})
    message(FATAL_ERROR "There is a ${STYLE_FILE} in ${LAMMPS_SOURCE_DIR}, please clean up the source directory first")
  endif()
endforeach()

enable_language(CXX)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################
include(CheckCCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" OFF)
option(INSTALL_LIB "Install lammps library and header" ON)
include(GNUInstallDirs)

set(LAMMPS_LINK_LIBS)
option(ENABLE_MPI "Build MPI version" OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_C_INCLUDE_PATH})
  list(APPEND LAMMPS_LINK_LIBS ${MPI_CXX_LIBRARIES})
  option(LAMMPS_LONGLONG_TO_LONG "Workaround if your system or MPI version does not recognize 'long long' data types" OFF)
  if(LAMMPS_LONGLONG_TO_LONG)
    add_definitions(-DLAMMPS_LONGLONG_TO_LONG)
  endif()
else()
  file(GLOB MPI_SOURCES ${LAMMPS_SOURCE_DIR}/STUBS/mpi.c)
  list(APPEND LIB_SOURCES ${MPI_SOURCES})
  include_directories(${LAMMPS_SOURCE_DIR}/STUBS)
endif()

set(LAMMPS_SIZE_LIMIT "LAMMPS_SMALLBIG" CACHE STRING "Lammps size limit")
set_property(CACHE LAMMPS_SIZE_LIMIT PROPERTY STRINGS LAMMPS_SMALLBIG LAMMPS_BIGBIG LAMMPS_SMALLSMALL)
add_definitions(-D${LAMMPS_SIZE_LIMIT})

set(LAMMPS_MEMALIGN "64" CACHE STRING "enables the use of the posix_memalign() call instead of malloc() when large chunks or memory are allocated by LAMMPS")
add_definitions(-DLAMMPS_MEMALIGN=${LAMMPS_MEMALIGN})

option(CMAKE_VERBOSE_MAKEFILE "Verbose makefile" OFF)

option(ENABLE_ALL "Build all packages" OFF)
set(PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE GRANULAR
  KIM KSPACE MANYBODY MC MEAM MISC MOLECULE MSCG MPIIO PERI POEMS PYTHON QEQ
  REAX REPLICA RIGID SHOCK SNAP SRD VORONOI)
set(USER-PACKAGES USER-ATC USER-AWPMD USER-CGDNA
  USER-CGSDK USER-COLVARS USER-DIFFRACTION USER-DPD USER-DRUDE USER-EFF
  USER-FEP USER-H5MD USER-LB USER-MANIFOLD USER-MEAMC USER-MGPT USER-MISC
  USER-MOLFILE USER-NETCDF USER-PHONON USER-QTB USER-REAXC USER-SMD 
  USER-SMTBQ USER-SPH USER-TALLY USER-VTK USER-QUIP USER-QMMM)
set(ACCEL_PACKAGES USER-OMP KOKKOS OPT USER-INTEL GPU)
foreach(PKG ${PACKAGES})
  option(ENABLE_${PKG} "Build ${PKG} Package" ${ENABLE_ALL})
endforeach()
foreach(PKG ${ACCEL_PACKAGES} ${USER-PACKAGES})
  option(ENABLE_${PKG} "Build ${PKG} Package" OFF)
endforeach()

macro(pkg_depends PKG1 PKG2)
  if(ENABLE_${PKG1} AND NOT ENABLE_${PKG2})
    message(FATAL_ERROR "${PKG1} package needs LAMMPS to be build with ${PKG2}")
  endif()
endmacro()

pkg_depends(MPIIO MPI)
pkg_depends(QEQ MANYBODY)
pkg_depends(USER-ATC MANYBODY)
pkg_depends(USER-H5MD MPI)
pkg_depends(USER-LB MPI)
pkg_depends(USER-MISC MANYBODY)
pkg_depends(USER-PHONON KSPACE)

######################################################
# packages with special compiler needs or external libs
######################################################
if(ENABLE_REAX OR ENABLE_MEAM OR ENABLE_USER-QUIP OR ENABLE_USER-QMMM)
  enable_language(Fortran)
endif()

if(ENABLE_KOKKOS OR ENABLE_MSCG)
  # starting with CMake 3.1 this is all you have to do to enforce C++11
  set(CMAKE_CXX_STANDARD 11) # C++11...
  set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
  set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
endif()

if(ENABLE_USER-OMP OR ENABLE_KOKKOS OR ENABLE_USER-INTEL)
  find_package(OpenMP REQUIRED)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(ENABLE_KSPACE)
  set(FFT "KISSFFT" CACHE STRING "FFT library for KSPACE package")
  set_property(CACHE FFT PROPERTY STRINGS KISSFFT FFTW3 MKL FFTW2)
  if(NOT FFT STREQUAL "KISSFFT")
    find_package(${FFT} REQUIRED)
    add_definitions(-DFFT_${FFT})
    include_directories(${${FFT}_INCLUDE_DIRS})
    list(APPEND LAMMPS_LINK_LIBS ${${FFT}_LIBRARIES})
  endif()
  set(PACK_OPTIMIZATION "PACK_ARRAY" CACHE STRING "Optimization for FFT")
  set_property(CACHE LAMMPS_SIZE_LIMIT PROPERTY STRINGS PACK_ARRAY PACK_POINTER PACK_MEMCPY)
  if(NOT PACK_OPTIMIZATION STREQUAL "PACK_ARRAY")
    add_definitions(-D${PACK_OPTIMIZATION})
  endif()
endif()

if(ENABLE_MISC)
  option(LAMMPS_XDR "include XDR compatibility files for doing particle dumps in XTC format" OFF)
  if(LAMMPS_XDR)
    add_definitions(-DLAMMPS_XDR)
  endif()
endif()

if(ENABLE_MSCG OR ENABLE_USER-ATC OR ENABLE_USER-AWPMD OR ENABLE_USER-QUIP)
  find_package(LAPACK)
  if(LAPACK_FOUND)
    list(APPEND LAMMPS_LINK_LIBS ${LAPACK_LIBRARIES})
  else()
    enable_language(Fortran)
    file(GLOB LAPACK_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/linalg/*.f)
    list(APPEND LIB_SOURCES ${LAPACK_SOURCES})
  endif()
endif()

if(ENABLE_PYTHON)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  add_definitions(-DLMP_PYTHON)
  include_directories(${PYTHON_INCLUDE_DIR})
  list(APPEND LAMMPS_LINK_LIBS ${PYTHON_LIBRARY})
  execute_process(COMMAND ${PYTHON_EXECUTABLE}
	  -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  install(FILES ${CMAKE_SOURCE_DIR}/../python/lammps.py DESTINATION ${PYTHON_INSTDIR})
  if(NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "Python package need lammps to be build shared, -DBUILD_SHARED_LIBS=ON")
  endif()
endif()

find_package(JPEG)
if(JPEG_FOUND)
  add_definitions(-DLAMMPS_JPEG)
  include_directories(${JPEG_INCLUDE_DIR})
  list(APPEND LAMMPS_LINK_LIBS ${JPEG_LIBRARIES})
endif()

find_package(PNG)
find_package(ZLIB)
if(PNG_FOUND AND ZLIB_FOUND)
  include_directories(${PNG_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${PNG_LIBRARIES} ${ZLIB_LIBRARIES})
  add_definitions(-DLAMMPS_PNG)
endif()

find_program(GZIP gzip)
if(GZIP)
  add_definitions(-DLAMMPS_GZIP)
endif()

find_program(FFMPEG ffmpeg)
if(FFMPEG)
  add_definitions(-DLAMMPS_FFMPEG)
endif()

if(ENABLE_VORONOI)
  find_package(VORO REQUIRED) #some distros
  include_directories(${VORO_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${VORO_LIBRARIES})
endif()

if(ENABLE_USER-MOLFILE)
  list(APPEND LAMMPS_LINK_LIBS ${CMAKE_DL_LIBS})
endif()

if(ENABLE_USER-NETCDF)
  find_package(NetCDF REQUIRED)
  include_directories(NETCDF_INCLUDE_DIR)
  list(APPEND LAMMPS_LINK_LIBS ${NETCDF_LIBRARY})
  add_definitions(-DLMP_HAS_NETCDF -DNC_64BIT_DATA=0x0020)
endif()

if(ENABLE_USER-SMD)
  find_package(Eigen3 REQUIRED)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

if(ENABLE_USER-QUIP)
  find_package(QUIP REQUIRED)
  list(APPEND LAMMPS_LINK_LIBS ${QUIP_LIBRARIES} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()

if(ENABLE_USER-QMMM)
  find_package(QE REQUIRED)
  include_directories(${QE_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${QE_LIBRARIES} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()

if(ENABLE_USER-AWPMD)
  include_directories(${LAMMPS_LIB_SOURCE_DIR}/awpmd/systems/interact
    ${LAMMPS_LIB_SOURCE_DIR}/awpmd/ivutils/include)
endif()

if(ENABLE_USER-H5MD)
  find_package(HDF5 REQUIRED)
  list(APPEND LAMMPS_LINK_LIBS ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIRS} ${LAMMPS_LIB_SOURCE_DIR}/h5md/include)
endif()

if(ENABLE_USER-VTK)
  find_package(VTK REQUIRED NO_MODULE)
  include(${VTK_USE_FILE})
  add_definitions(-DLAMMPS_VTK)
  list(APPEND LAMMPS_LINK_LIBS ${VTK_LIBRARIES})
endif()

if(ENABLE_KIM)
  find_package(KIM REQUIRED)
  list(APPEND LAMMPS_LINK_LIBS ${KIM_LIBRARIES})
  include_directories(${KIM_INCLUDE_DIRS})
endif()

if(ENABLE_MSCG)
  find_package(GSL REQUIRED)
  set(LAMMPS_LIB_MSCG_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/mscg)
  set(MSCG_TARBALL ${LAMMPS_LIB_MSCG_BIN_DIR}/MS-CG-master.zip)
  set(LAMMPS_LIB_MSCG_BIN_DIR ${LAMMPS_LIB_MSCG_BIN_DIR}/MSCG-release-master/src)
  if(NOT EXISTS ${LAMMPS_LIB_MSCG_BIN_DIR})
    if(NOT EXISTS ${MSCG_TARBALL})
      message(STATUS "Downloading ${MSCG_TARBALL}")
      file(DOWNLOAD
        https://github.com/uchicago-voth/MSCG-release/archive/master.zip
        ${MSCG_TARBALL} SHOW_PROGRESS) #EXPECTED_MD5 cannot be due due to master
    endif()
    message(STATUS "Unpacking ${MSCG_TARBALL}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ${MSCG_TARBALL}
      WORKING_DIRECTORY ${LAMMPS_LIB_BINARY_DIR}/mscg)
  endif()
  file(GLOB MSCG_SOURCES ${LAMMPS_LIB_MSCG_BIN_DIR}/*.cpp)
  list(APPEND LIB_SOURCES ${MSCG_SOURCES})
  foreach(MSCG_SOURCE ${MSCG_SOURCES})
    set_property(SOURCE ${MSCG_SOURCE} APPEND PROPERTY COMPILE_DEFINITIONS
      DIMENSION=3 _exclude_gromacs=1)
  endforeach()
  include_directories(${LAMMPS_LIB_MSCG_BIN_DIR} ${GSL_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${GSL_LIBRARIES})
endif()

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER math.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sin cos)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)
list(APPEND LAMMPS_LINK_LIBS ${MATH_LIBRARIES})

######################################
# Generate Basic Style files 
######################################
include(StyleHeaderUtils)
RegisterStyles(${LAMMPS_SOURCE_DIR})

##############################################
# add sources of enabled packages
############################################
foreach(PKG ${PACKAGES} ${USER-PACKAGES})
  if(ENABLE_${PKG})
    set(${PKG}_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/${PKG})

    # detects styles in package and adds them to global list
    RegisterStyles(${${PKG}_SOURCES_DIR})

    file(GLOB ${PKG}_SOURCES ${${PKG}_SOURCES_DIR}/*.cpp)
    list(APPEND LIB_SOURCES ${${PKG}_SOURCES})
    include_directories(${${PKG}_SOURCES_DIR})
  endif()
endforeach()

##############################################
# add lib sources of (simple) enabled packages
############################################
foreach(SIMPLE_LIB REAX MEAM POEMS USER-ATC USER-AWPMD USER-COLVARS USER-H5MD
  USER-MOLFILE USER-QMMM)
  if(ENABLE_${SIMPLE_LIB})
    string(REGEX REPLACE "^USER-" "" SIMPLE_LIB "${SIMPLE_LIB}")
    string(TOLOWER "${SIMPLE_LIB}" INC_DIR)
    file(GLOB_RECURSE ${SIMPLE_LIB}_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/${INC_DIR}/*.F
      ${LAMMPS_LIB_SOURCE_DIR}/${INC_DIR}/*.c ${LAMMPS_LIB_SOURCE_DIR}/${INC_DIR}/*.cpp)
    list(APPEND LIB_SOURCES ${${SIMPLE_LIB}_SOURCES})
    include_directories(${LAMMPS_LIB_SOURCE_DIR}/${INC_DIR})
  endif()
endforeach()

######################################################################
# packages which selectively include variants based on enabled styles
# e.g. accelerator packages
######################################################################
if(ENABLE_USER-OMP)
    set(USER-OMP_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-OMP)
    set(USER-OMP_SOURCES ${USER-OMP_SOURCES_DIR}/thr_data.cpp
                         ${USER-OMP_SOURCES_DIR}/thr_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/fix_nh_omp.cpp
                         ${USER-OMP_SOURCES_DIR}/fix_nh_sphere_omp.cpp)
    set_property(GLOBAL PROPERTY "OMP_SOURCES" "${USER-OMP_SOURCES}")

    # detects styles which have USER-OMP version
    RegisterStylesExt(${USER-OMP_SOURCES_DIR} omp OMP_SOURCES)

    get_property(USER-OMP_SOURCES GLOBAL PROPERTY OMP_SOURCES)

    list(APPEND LIB_SOURCES ${USER-OMP_SOURCES})
    include_directories(${USER-OMP_SOURCES_DIR})
endif()

if(ENABLE_KOKKOS)
  set(LAMMPS_LIB_KOKKOS_SRC_DIR ${LAMMPS_LIB_SOURCE_DIR}/kokkos)
  set(LAMMPS_LIB_KOKKOS_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/kokkos)
  add_definitions(-DLMP_KOKKOS)
  add_subdirectory(${LAMMPS_LIB_KOKKOS_SRC_DIR} ${LAMMPS_LIB_KOKKOS_BIN_DIR})

  set(Kokkos_INCLUDE_DIRS ${LAMMPS_LIB_KOKKOS_SRC_DIR}/core/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/containers/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/algorithms/src
                          ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  include_directories(${Kokkos_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS kokkos)

  set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
  set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/atom_vec_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/comm_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/comm_tiled_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neighbor_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neigh_list_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/neigh_bond_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/fix_nh_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/domain_kokkos.cpp
                         ${KOKKOS_PKG_SOURCES_DIR}/modify_kokkos.cpp)
  set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")

  # detects styles which have KOKKOS version
  RegisterStylesExt(${KOKKOS_PKG_SOURCES_DIR} kokkos KOKKOS_PKG_SOURCES)

  get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)

  list(APPEND LIB_SOURCES ${KOKKOS_PKG_SOURCES})
  include_directories(${KOKKOS_PKG_SOURCES_DIR})
endif()

if(ENABLE_OPT)
    set(OPT_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/OPT)
    set(OPT_SOURCES)
    set_property(GLOBAL PROPERTY "OPT_SOURCES" "${OPT_SOURCES}")

    # detects styles which have OPT version
    RegisterStylesExt(${OPT_SOURCES_DIR} opt OPT_SOURCES)

    get_property(OPT_SOURCES GLOBAL PROPERTY OPT_SOURCES)

    list(APPEND LIB_SOURCES ${OPT_SOURCES})
    include_directories(${OPT_SOURCES_DIR})
endif()

if(ENABLE_USER-INTEL)
    set(USER-INTEL_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-INTEL)
    set(USER-INTEL_SOURCES ${USER-INTEL_SOURCES_DIR}/intel_preprocess.h
                           ${USER-INTEL_SOURCES_DIR}/intel_buffers.h
                           ${USER-INTEL_SOURCES_DIR}/intel_buffers.cpp
                           ${USER-INTEL_SOURCES_DIR}/math_extra_intel.h
                           ${USER-INTEL_SOURCES_DIR}/nbin_intel.h
                           ${USER-INTEL_SOURCES_DIR}/nbin_intel.cpp
                           ${USER-INTEL_SOURCES_DIR}/npair_intel.h
                           ${USER-INTEL_SOURCES_DIR}/npair_intel.cpp
                           ${USER-INTEL_SOURCES_DIR}/intel_simd.h
                           ${USER-INTEL_SOURCES_DIR}/intel_intrinsics.h)

    set_property(GLOBAL PROPERTY "USER-INTEL_SOURCES" "${USER-INTEL_SOURCES}")

    # detects styles which have USER-INTEL version
    RegisterStylesExt(${USER-INTEL_SOURCES_DIR} opt USER-INTEL_SOURCES)

    get_property(USER-INTEL_SOURCES GLOBAL PROPERTY USER-INTEL_SOURCES)

    list(APPEND LIB_SOURCES ${USER-INTEL_SOURCES})
    include_directories(${USER-INTEL_SOURCES_DIR})
endif()

if(ENABLE_GPU)
    find_package(CUDA REQUIRED)
    find_program(BIN2C bin2c)
    if(NOT BIN2C)
      message(FATAL_ERROR "Couldn't find bin2c")
    endif()
    include_directories(${CUDA_TOOLKIT_INCLUDE})
    set(GPU_PREC "SINGLE_DOUBLE" CACHE STRING "Lammps gpu precision size")
    set_property(CACHE GPU_PREC PROPERTY STRINGS SINGLE_DOUBLE SINGLE_SINGLE DOUBLE_DOUBLE)
    add_definitions(-D_${GPU_PREC})
    add_definitions(-DNV_KERNEL -DUSE_CUDPP -DUCL_CUDADR)

    set(GPU_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/GPU)
    set(GPU_SOURCES ${GPU_SOURCES_DIR}/gpu_extra.h)

    set_property(GLOBAL PROPERTY "GPU_SOURCES" "${GPU_SOURCES}")

    # detects styles which have GPU version
    RegisterStylesExt(${GPU_SOURCES_DIR} opt GPU_SOURCES)

    get_property(GPU_SOURCES GLOBAL PROPERTY GPU_SOURCES)
    
    file(GLOB_RECURSE GPU_LIB_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/gpu/*.cpp)
    file(GLOB_RECURSE GPU_LIB_CU ${LAMMPS_LIB_SOURCE_DIR}/gpu/*.cu ${CMAKE_SOURCE_DIR}/gpu/*.cu)
    file(GLOB_RECURSE GPU_NOT_LIB_CU ${LAMMPS_LIB_SOURCE_DIR}/gpu/lal_pppm.cu)
    list(REMOVE_ITEM GPU_LIB_CU ${GPU_NOT_LIB_CU})
    add_custom_target(gpu_objs)
    cuda_wrap_srcs(gpu_objs OBJ GPU_OBJS ${GPU_LIB_CU})
    file(MAKE_DIRECTORY ${LAMMPS_LIB_BINARY_DIR}/gpu)
    foreach(CU_OBJ ${GPU_OBJS})
      get_filename_component(CU_NAME ${CU_OBJ} NAME_WE)
      string(REGEX REPLACE "^.*_lal_" "" CU_NAME "${CU_NAME}")
      add_custom_command(OUTPUT ${LAMMPS_LIB_BINARY_DIR}/gpu/${CU_NAME}_cubin.h
        COMMAND ${BIN2C} -c -n ${CU_NAME} ${CU_OBJ} > ${LAMMPS_LIB_BINARY_DIR}/gpu/${CU_NAME}_cubin.h
        DEPENDS ${CU_OBJ} 
        COMMENT "Generating ${CU_NAME}_cubin.h")
      list(APPEND LIB_SOURCES ${LAMMPS_LIB_BINARY_DIR}/gpu/${CU_NAME}_cubin.h)
    endforeach()
    list(APPEND LIB_SOURCES ${GPU_SOURCES} ${GPU_LIB_SOURCES} ${GPU_OBJS})
    include_directories(${GPU_SOURCES_DIR} ${LAMMPS_LIB_SOURCE_DIR}/gpu
      ${LAMMPS_LIB_SOURCE_DIR}/gpu/cudpp_mini ${LAMMPS_LIB_BINARY_DIR}/gpu)
    set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${AMMPS_LIB_BINARY_DIR}/gpu/*_cubin.h")
endif()

######################################################
# Generate style headers based on global list of
# styles registered during package selection           
######################################################
set(LAMMPS_STYLE_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/styles)

GenerateStyleHeaders(${LAMMPS_STYLE_HEADERS_DIR})

include_directories(${LAMMPS_SOURCE_DIR})
include_directories(${LAMMPS_STYLE_HEADERS_DIR})

###########################################
# Actually add executable and lib to build
############################################
add_library(lammps ${LIB_SOURCES})
target_link_libraries(lammps ${LAMMPS_LINK_LIBS})
set_target_properties(lammps PROPERTIES SOVERSION ${SOVERSION})
if(INSTALL_LIB)
  install(TARGETS lammps LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES ${LAMMPS_SOURCE_DIR}/lammps.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
elseif(BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Shared library has to be installed, use -DBUILD_SHARED_LIBS=ON to install lammps with a library")
endif()

add_executable(lmp ${LMP_SOURCES})
target_link_libraries(lmp lammps)
install(TARGETS lmp DESTINATION ${CMAKE_INSTALL_BINDIR})

##################################
# Print package summary
##################################
foreach(PKG ${PACKAGES} ${USER-PACKAGES} ${ACCEL_PACKAGES})
  if(ENABLE_${PKG})
    message(STATUS "Building package: ${PKG}")
  endif()
endforeach()
