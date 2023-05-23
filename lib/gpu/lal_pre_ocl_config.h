// **************************************************************************
//                               pre_ocl_config.h
//                             -------------------
//                           W. Michael Brown (ORNL)
//                           Nitin Dhamankar (Intel)
//
//  Device-side preprocessor definitions
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

//*************************************************************************
//                       Device Configuration Definitions
//                    See lal_preprocessor.h for definitions
//                           Configuration order:
//
//  {CONFIG_NAME, CONFIG_ID, SIMD_SIZE, MEM_THREADS, SHUFFLE_AVAIL, FAST_MATH,
//   THREADS_PER_ATOM, THREADS_PER_CHARGE, THREADS_PER_THREE, BLOCK_PAIR,
//   BLOCK_BIO_PAIR, BLOCK_ELLIPSE, PPPM_BLOCK_1D, BLOCK_NBOR_BUILD,
//   BLOCK_CELL_2D, BLOCK_CELL_ID, MAX_SHARED_TYPES, MAX_BIO_SHARED_TYPES,
//   PPPM_MAX_SPLINE, NBOR_PREFETCH}
//
//*************************************************************************/

const int nconfigs=6;
const char * ocl_config_names[] =
  {
   "generic",
   "nvidiagpu",
   "amdgpu",
   "intelgpu",
   "applegpu",
   "intelcpu"
  };
const char * ocl_config_strings[] =
  {
   "GENERIC,1,1,16,0,1,1,1,1,64,64,64,64,64,8,128,8,128,8,0",
   "NVIDIA_GPU,203,32,32,1,1,4,8,2,256,256,128,64,128,8,128,11,128,8,0",
   "AMD_GPU,403,64,64,0,1,4,8,2,256,256,128,64,128,8,128,11,128,8,0",
#ifdef _SINGLE_SINGLE
   "INTEL_GPU,500,8,32,1,1,4,8,2,128,128,128,128,64,8,128,8,128,8,2",
   "APPLE_GPU,600,16,16,0,1,4,8,1,64,64,64,64,64,8,128,8,128,8,0",
#else
   "INTEL_GPU,500,8,32,1,1,2,8,2,128,128,128,128,64,8,128,8,128,8,2",
   "APPLE_GPU,600,16,16,0,1,2,8,1,64,64,64,64,64,8,128,8,128,8,0",
#endif
   "INTEL_CPU,1500,8,8,1,1,1,1,1,64,64,64,64,64,8,64,8,128,8,0"
  };
