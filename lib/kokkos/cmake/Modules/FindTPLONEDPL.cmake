include(CheckIncludeFileCXX)
check_include_file_cxx(oneapi/dpl/execution KOKKOS_COMPILER_HAS_ONEDPL_EXECUTION_HEADER)
check_include_file_cxx(oneapi/dpl/algorithm KOKKOS_COMPILER_HAS_ONEDPL_ALGORITHM_HEADER)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "
  #include <iostream>

  int main()
  {
    #if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE == 9 || _GLIBCXX_RELEASE == 10)
      static_assert(false);
    #endif
    return 0;
  }"
  KOKKOS_NO_TBB_CONFLICT
)

if(KOKKOS_COMPILER_HAS_ONEDPL_EXECUTION_HEADER AND KOKKOS_COMPILER_HAS_ONEDPL_ALGORITHM_HEADER)
  if(KOKKOS_NO_TBB_CONFLICT)
    kokkos_create_imported_tpl(ONEDPL INTERFACE)
  else()
    kokkos_create_imported_tpl(
      ONEDPL
      INTERFACE
      # https://stackoverflow.com/questions/67923287/how-to-resolve-no-member-named-task-in-namespace-tbb-error-when-using-oned/
      COMPILE_DEFINITIONS
      PSTL_USE_PARALLEL_POLICIES=0
      _GLIBCXX_USE_TBB_PAR_BACKEND=0
    )
  endif()
else()
  find_package(oneDPL REQUIRED)

  if(KOKKOS_NO_TBB_CONFLICT)
    kokkos_create_imported_tpl(ONEDPL INTERFACE LINK_LIBRARIES oneDPL)
  else()
    kokkos_create_imported_tpl(
      ONEDPL
      INTERFACE
      LINK_LIBRARIES
      oneDPL
      # https://stackoverflow.com/questions/67923287/how-to-resolve-no-member-named-task-in-namespace-tbb-error-when-using-oned/
      COMPILE_DEFINITIONS
      PSTL_USE_PARALLEL_POLICIES=0
      _GLIBCXX_USE_TBB_PAR_BACKEND=0
    )
  endif()

  # Export oneDPL as a Kokkos dependency
  kokkos_export_cmake_tpl(oneDPL)
endif()
