message(STATUS "Downloading and building OpenCL loader library")
set(OPENCL_LOADER_URL "${LAMMPS_THIRDPARTY_URL}/opencl-loader-2024.05.09.tar.gz" CACHE STRING "URL for OpenCL loader tarball")
set(OPENCL_LOADER_SHA256 "e33f78a92bbacc2c8639cb2f00b347b4715c8b57aaac2c14dc2cd86e836cfd34" CACHE STRING "SHA256 checksum of OpenCL loader tarball")
mark_as_advanced(OPENCL_LOADER_URL)
mark_as_advanced(OPENCL_LOADER_SHA256)

set(INSTALL_LIBOPENCL OFF CACHE BOOL "" FORCE)
include(ExternalCMakeProject)
ExternalCMakeProject(opencl_loader ${OPENCL_LOADER_URL} ${OPENCL_LOADER_SHA256} opencl-loader . "")

