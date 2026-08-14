# GFN1-xTB and GFN2-xTB QM/MM coupling with a self-consistent periodic shift.
# libxtb must be built with INSTALL_MODULES=ON using the same Fortran compiler
# as LAMMPS, because the public C API cannot inject an atom-dependent shift at
# every SCC iteration.

if(NOT PKG_KSPACE)
  message(FATAL_ERROR "The QMMM-XTB package requires the KSPACE package")
endif()

enable_language(Fortran)
find_package(PkgConfig REQUIRED)
pkg_check_modules(XTB REQUIRED IMPORTED_TARGET "xtb>=6.7")
pkg_check_modules(MCTC REQUIRED IMPORTED_TARGET mctc-lib)
find_package(BLAS REQUIRED)

set(XTB_FORTRAN_MODULE_DIR "" CACHE PATH
  "Directory containing the libxtb Fortran .mod files (built with INSTALL_MODULES=ON)")

if(NOT XTB_FORTRAN_MODULE_DIR)
  foreach(_xtb_inc ${XTB_INCLUDE_DIRS})
    if(EXISTS "${_xtb_inc}/xtb_type_solvation.mod")
      set(XTB_FORTRAN_MODULE_DIR "${_xtb_inc}")
      break()
    endif()
    file(GLOB _xtb_mod_candidates LIST_DIRECTORIES TRUE
      "${_xtb_inc}/*" "${_xtb_inc}/*/*")
    foreach(_xtb_candidate ${_xtb_mod_candidates})
      if(EXISTS "${_xtb_candidate}/xtb_type_solvation.mod")
        set(XTB_FORTRAN_MODULE_DIR "${_xtb_candidate}")
        break()
      endif()
    endforeach()
  endforeach()
endif()

if(NOT EXISTS "${XTB_FORTRAN_MODULE_DIR}/xtb_type_solvation.mod")
  message(FATAL_ERROR
    "QMMM-XTB requires compatible libxtb 6.7 or newer Fortran modules. Rebuild xTB with "
    "-Dinstall_modules=true and set XTB_FORTRAN_MODULE_DIR to their include directory.")
endif()

add_library(qmmm_xtb_adapter STATIC
  ${LAMMPS_SOURCE_DIR}/QMMM-XTB/qmmm_xtb_adapter.f90)
set_target_properties(qmmm_xtb_adapter PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/qmmm_xtb_modules)
target_include_directories(qmmm_xtb_adapter PRIVATE ${XTB_FORTRAN_MODULE_DIR})
target_link_libraries(qmmm_xtb_adapter PUBLIC PkgConfig::XTB PkgConfig::MCTC BLAS::BLAS)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # xTB builds its GNU Fortran interface with these ABI-affecting options.
  # Compile the adapter identically so derived types imported from .mod files
  # have the same default-kind layout as the linked library.
  target_compile_options(qmmm_xtb_adapter PRIVATE
    -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none)
endif()

target_link_libraries(lammps PRIVATE qmmm_xtb_adapter)
