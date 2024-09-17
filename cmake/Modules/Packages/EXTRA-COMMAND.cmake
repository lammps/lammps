# the geturl command needs libcurl

find_package(CURL QUIET COMPONENTS HTTP HTTPS)
option(WITH_CURL "Enable libcurl support" ${CURL_FOUND})
if(WITH_CURL)
  find_package(CURL REQUIRED COMPONENTS HTTP HTTPS)
  target_compile_definitions(lammps PRIVATE -DLAMMPS_CURL)
  target_link_libraries(lammps PRIVATE CURL::libcurl)
endif()

