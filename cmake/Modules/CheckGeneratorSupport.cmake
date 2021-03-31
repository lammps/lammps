# ninja-build<1.10 does not support fortran.
if(CMAKE_GENERATOR STREQUAL "Ninja")
  set(CMAKE_GENERATOR_SUPPORT_FORTRAN FALSE)
  execute_process(COMMAND "${CMAKE_MAKE_PROGRAM}" --version
    OUTPUT_VARIABLE NINJA_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _Ninja_version_result
  )
  if(_Ninja_version_result)
    message(WARNING "Unable to determine ninja version: ${_Ninja_version_result}, assuming fortran isn't supported")
  elseif(NINJA_VERSION VERSION_LESS "1.10")
    message(WARNING "Ninja build tool too old, to compile Fortran code, please install ninja-1.10 or newer")
  else()
    set(CMAKE_GENERATOR_SUPPORT_FORTRAN TRUE)
  endif()
else()
  set(CMAKE_GENERATOR_SUPPORT_FORTRAN TRUE)
  if(NOT CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    message(WARNING "Assuming fortran is supported for ${CMAKE_GENERATOR}")
  endif()
endif()
