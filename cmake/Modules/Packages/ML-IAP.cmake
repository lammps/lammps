# if PYTHON package is included we may also include Python support in ML-IAP
set(MLIAP_ENABLE_PYTHON_DEFAULT OFF)
if(PKG_PYTHON)
  find_package(Cythonize QUIET)
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
    find_package(Python COMPONENTS NumPy QUIET)
  else()
    # assume we have NumPy
    set(Python_NumPy_FOUND ON)
  endif()
  if(Cythonize_FOUND AND Python_NumPy_FOUND)
    set(MLIAP_ENABLE_PYTHON_DEFAULT ON)
  endif()
endif()

option(MLIAP_ENABLE_PYTHON "Build ML-IAP package with Python support" ${MLIAP_ENABLE_PYTHON_DEFAULT})

if(MLIAP_ENABLE_PYTHON)
  find_package(Cythonize REQUIRED)
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
    find_package(Python COMPONENTS NumPy REQUIRED)
  endif()
  if(NOT PKG_PYTHON)
    message(FATAL_ERROR "Must enable PYTHON package for including Python support in ML-IAP")
  endif()
  if(CMAKE_VERSION VERSION_LESS 3.12)
    if(PYTHONLIBS_VERSION_STRING VERSION_LESS 3.6)
      message(FATAL_ERROR "Python support in ML-IAP requires Python 3.6 or later")
    endif()
  else()
    if(Python_VERSION VERSION_LESS 3.6)
      message(FATAL_ERROR "Python support in ML-IAP requires Python 3.6 or later")
    endif()
  endif()

  set(MLIAP_BINARY_DIR ${CMAKE_BINARY_DIR}/cython)
  file(GLOB MLIAP_CYTHON_SRC ${CONFIGURE_DEPENDS} ${LAMMPS_SOURCE_DIR}/ML-IAP/*.pyx)
  file(MAKE_DIRECTORY ${MLIAP_BINARY_DIR})
  foreach(MLIAP_CYTHON_FILE ${MLIAP_CYTHON_SRC})
    get_filename_component(MLIAP_CYTHON_BASE ${MLIAP_CYTHON_FILE} NAME_WE)
    add_custom_command(OUTPUT  ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.h
            COMMAND            ${CMAKE_COMMAND} -E copy_if_different ${MLIAP_CYTHON_FILE} ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
            COMMAND            ${Cythonize_EXECUTABLE} -3 ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.pyx
            WORKING_DIRECTORY  ${MLIAP_BINARY_DIR}
            MAIN_DEPENDENCY    ${MLIAP_CYTHON_FILE}
            COMMENT "Generating C++ sources with cythonize...")
    target_sources(lammps PRIVATE ${MLIAP_BINARY_DIR}/${MLIAP_CYTHON_BASE}.cpp)
  endforeach()
  target_compile_definitions(lammps PRIVATE -DMLIAP_PYTHON)
  target_include_directories(lammps PRIVATE ${MLIAP_BINARY_DIR})
endif()
