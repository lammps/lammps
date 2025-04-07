# the geturl command needs libcurl

find_package(CURL QUIET)
option(WITH_CURL "Enable libcurl support" ${CURL_FOUND})
if(WITH_CURL)
  target_compile_definitions(lammps PRIVATE -DLAMMPS_CURL)

  # need to use pkgconfig for fully static bins to find custom static libs
  if (CMAKE_SYSTEM_NAME STREQUAL "LinuxMUSL")
    include(FindPkgConfig)
    pkg_check_modules(CURL IMPORTED_TARGET libcurl libssl libcrypto)
    target_link_libraries(lammps PUBLIC PkgConfig::CURL)
  else()
    find_package(CURL REQUIRED)
    target_link_libraries(lammps PRIVATE CURL::libcurl)
  endif()
endif()

