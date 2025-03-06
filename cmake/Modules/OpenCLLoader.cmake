message(STATUS "Downloading and building OpenCL loader library")
set(OPENCL_LOADER_URL "${LAMMPS_THIRDPARTY_URL}/opencl-loader-2024.05.09.tar.gz" CACHE STRING "URL for OpenCL loader tarball")
set(OPENCL_LOADER_MD5 "e7796826b21c059224fabe997e0f2075" CACHE STRING "MD5 checksum of OpenCL loader tarball")
mark_as_advanced(OPENCL_LOADER_URL)
mark_as_advanced(OPENCL_LOADER_MD5)

set(INSTALL_LIBOPENCL OFF CACHE BOOL "" FORCE)
include(ExternalCMakeProject)
ExternalCMakeProject(opencl_loader ${OPENCL_LOADER_URL} ${OPENCL_LOADER_MD5} opencl-loader . "")

