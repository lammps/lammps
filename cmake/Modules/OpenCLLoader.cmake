message(STATUS "Downloading and building OpenCL loader library")
SetDownloadSettings(OPENCL_LOADER "OpenCL loader"
  "${LAMMPS_THIRDPARTY_URL}/opencl-loader-2024.05.09.tar.gz"
  "e33f78a92bbacc2c8639cb2f00b347b4715c8b57aaac2c14dc2cd86e836cfd34")

set(INSTALL_LIBOPENCL OFF CACHE BOOL "" FORCE)
include(ExternalCMakeProject)
ExternalCMakeProject(opencl_loader ${OPENCL_LOADER_URL} ${OPENCL_LOADER_SHA256} opencl-loader . "")

