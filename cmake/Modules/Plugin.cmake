function(lammps_add_plugin target)
  add_library(${target} MODULE ${ARGN})
  target_link_libraries(${target} PRIVATE LAMMPS::lammps)

  set_target_properties(${target} PROPERTIES PREFIX "" SUFFIX ".so")

  # MacOS seems to need this
  if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
    target_link_options(${target} PRIVATE "-Wl,-undefined,dynamic_lookup")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    # tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
    set_target_properties(${target} PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    if(CMAKE_CROSSCOMPILING)
      target_link_options(${target} PRIVATE "-Wl,--export-all-symbols")
    endif()
  else()
    target_link_options(${target} PRIVATE "-rdynamic")
  endif()
endfunction()
