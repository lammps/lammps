/***************************************************************************
                                nvd_device.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with cuda devices

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Thu Jan 21 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef NVD_DEVICE
#define NVD_DEVICE

#include <string>
#include <vector>
#include <iostream>
#include "nvd_macros.h"
#include "ucl_types.h"

namespace ucl_cudadr {

// --------------------------------------------------------------------------
// - COMMAND QUEUE STUFF
// --------------------------------------------------------------------------
typedef CUstream command_queue;

inline void ucl_flush(command_queue &cq) {}

inline void ucl_sync(CUstream &stream) {
  CU_SAFE_CALL(cuStreamSynchronize(stream));
}

struct NVDProperties {
  int device_id;
  std::string name;
  int major;
  int minor;
  CUDA_INT_TYPE totalGlobalMem;
  int multiProcessorCount;

  int maxThreadsPerBlock;
  int maxThreadsDim[3];
  int maxGridSize[3];
  int sharedMemPerBlock;
  int totalConstantMemory;
  int SIMDWidth;
  int memPitch;
  int regsPerBlock;
  int clockRate;
  int textureAlign;

  int kernelExecTimeoutEnabled;
  int integrated;
  int canMapHostMemory;
  int concurrentKernels;
  int ECCEnabled;
  int computeMode;
};

/// Class for looking at device properties
/** \note Calls to change the device outside of the class results in incorrect
  *       behavior
  * \note There is no error checking for indexing past the number of devices **/
class UCL_Device {
 public:
  /// Collect properties for every GPU on the node
  /** \note You must set the active GPU with set() before using the device **/
  inline UCL_Device();

  inline ~UCL_Device();

  /// Returns 1 (For compatibility with OpenCL)
  inline int num_platforms() { return 1; }

  /// Return a string with name and info of the current platform
  inline std::string platform_name()
    { return "NVIDIA Corporation NVIDIA CUDA Driver"; }

  /// Delete any contexts/data and set the platform number to be used
  inline int set_platform(const int pid);

  /// Return the number of devices that support CUDA
  inline int num_devices() { return _properties.size(); }

  /// Set the CUDA device to the specified device number
  /** A context and default command queue will be created for the device
    * Returns UCL_SUCCESS if successful or UCL_ERROR if the device could not
    * be allocated for use. clear() is called to delete any contexts and
    * associated data from previous calls to set(). **/
  inline int set(int num);

  /// Delete any context and associated data stored from a call to set()
  inline void clear();

  /// Get the current device number
  inline int device_num() { return _device; }

  /// Returns the default stream for the current device
  inline command_queue & cq() { return cq(0); }

  /// Returns the stream indexed by i
  inline command_queue & cq(const int i) { return _cq[i]; }

  /// Block until all commands in the default stream have completed
  inline void sync() { sync(0); }

  /// Block until all commands in the specified stream have completed
  inline void sync(const int i) { ucl_sync(cq(i)); }

  /// Get the number of command queues currently available on device
  inline int num_queues()
    { return _cq.size(); }

  /// Add a stream for device computations
  inline void push_command_queue() {
    _cq.push_back(CUstream());
    CU_SAFE_CALL(cuStreamCreate(&_cq.back(),0));
  }

  /// Remove a stream for device computations
  /** \note You cannot delete the default stream **/
  inline void pop_command_queue() {
    if (_cq.size()<2) return;
    CU_SAFE_CALL_NS(cuStreamDestroy(_cq.back()));
    _cq.pop_back();
  }

  /// Set the default command queue (by default this is the null stream)
  /** \param i index of the command queue (as added by push_command_queue())
      If i is 0, the default command queue is set to the null stream **/
  inline void set_command_queue(const int i) {
    if (i==0) _cq[0]=0;
    else _cq[0]=_cq[i];
  }

  /// Get the current CUDA device name
  inline std::string name() { return name(_device); }
  /// Get the CUDA device name
  inline std::string name(const int i)
    { return std::string(_properties[i].name); }

  /// Get a string telling the type of the current device
  inline std::string device_type_name() { return device_type_name(_device); }
  /// Get a string telling the type of the device
  inline std::string device_type_name(const int i) { return "GPU"; }

  /// Get current device type (UCL_CPU, UCL_GPU, UCL_ACCELERATOR, UCL_DEFAULT)
  inline enum UCL_DEVICE_TYPE device_type() { return device_type(_device); }
  /// Get device type (UCL_CPU, UCL_GPU, UCL_ACCELERATOR, UCL_DEFAULT)
  inline enum UCL_DEVICE_TYPE device_type(const int i) { return UCL_GPU; }

  /// Returns true if host memory is efficiently addressable from device
  inline bool shared_memory() { return shared_memory(_device); }
  /// Returns true if host memory is efficiently addressable from device
  inline bool shared_memory(const int i) { return device_type(i)==UCL_CPU; }

  /// Returns preferred vector width
  inline int preferred_fp32_width() { return preferred_fp32_width(_device); }
  /// Returns preferred vector width
  inline int preferred_fp32_width(const int i)
    {return _properties[i].SIMDWidth;}
  /// Returns preferred vector width
  inline int preferred_fp64_width() { return preferred_fp64_width(_device); }
  /// Returns preferred vector width
  inline int preferred_fp64_width(const int i)
    {return _properties[i].SIMDWidth;}

  /// Returns true if double precision is support for the current device
  inline bool double_precision() { return double_precision(_device); }
  /// Returns true if double precision is support for the device
  inline bool double_precision(const int i) {return arch(i)>=1.3;}

  /// Get the number of compute units on the current device
  inline unsigned cus() { return cus(_device); }
  /// Get the number of compute units
  inline unsigned cus(const int i)
    { return _properties[i].multiProcessorCount; }

  /// Get the number of cores in the current device
  inline unsigned cores() { return cores(_device); }
  /// Get the number of cores
  inline unsigned cores(const int i)
    { if (arch(i)<2.0) return _properties[i].multiProcessorCount*8;
      else if (arch(i)<2.1) return _properties[i].multiProcessorCount*32;
      else if (arch(i)<3.0) return _properties[i].multiProcessorCount*48;
      else return _properties[i].multiProcessorCount*192; }

  /// Get the gigabytes of global memory in the current device
  inline double gigabytes() { return gigabytes(_device); }
  /// Get the gigabytes of global memory
  inline double gigabytes(const int i)
    { return static_cast<double>(_properties[i].totalGlobalMem)/1073741824; }

  /// Get the bytes of global memory in the current device
  inline size_t bytes() { return bytes(_device); }
  /// Get the bytes of global memory
  inline size_t bytes(const int i) { return _properties[i].totalGlobalMem; }

  // Get the gigabytes of free memory in the current device
  inline double free_gigabytes() { return free_gigabytes(_device); }
  // Get the gigabytes of free memory
  inline double free_gigabytes(const int i)
    { return static_cast<double>(free_bytes(i))/1073741824; }

  // Get the bytes of free memory in the current device
  inline size_t free_bytes() { return free_bytes(_device); }
  // Get the bytes of free memory
  inline size_t free_bytes(const int i) {
    CUDA_INT_TYPE dfree, dtotal;
    CU_SAFE_CALL_NS(cuMemGetInfo(&dfree, &dtotal));
    return static_cast<size_t>(dfree);
  }

  /// Return the GPGPU compute capability for current device
  inline double arch() { return arch(_device); }
  /// Return the GPGPU compute capability
  inline double arch(const int i)
    { return static_cast<double>(_properties[i].minor)/10+_properties[i].major;}

  /// Clock rate in GHz for current device
  inline double clock_rate() { return clock_rate(_device); }
  /// Clock rate in GHz
  inline double clock_rate(const int i)
    { return _properties[i].clockRate*1e-6;}

  /// Get the maximum number of threads per block
  inline size_t group_size() { return group_size(_device); }
  /// Get the maximum number of threads per block
  inline size_t group_size(const int i)
    { return _properties[i].maxThreadsPerBlock; }
  /// Get the maximum number of threads per block in dimension 'dim'
  inline size_t group_size_dim(const int dim)
    { return group_size_dim(_device, dim); }
  /// Get the maximum number of threads per block in dimension 'dim'
  inline size_t group_size_dim(const int i, const int dim)
    { return _properties[i].maxThreadsDim[dim]; }

  /// Get the shared local memory size in bytes
  inline size_t slm_size() { return slm_size(_device); }
  /// Get the shared local memory size in bytes
  inline size_t slm_size(const int i)
    { return _properties[i].sharedMemPerBlock; }

  /// Return the maximum memory pitch in bytes for current device
  inline size_t max_pitch() { return max_pitch(_device); }
  /// Return the maximum memory pitch in bytes
  inline size_t max_pitch(const int i) { return _properties[i].memPitch; }

  /// Returns false if accelerator cannot be shared by multiple processes
  /** If it cannot be determined, true is returned **/
  inline bool sharing_supported() { return sharing_supported(_device); }
  /// Returns false if accelerator cannot be shared by multiple processes
  /** If it cannot be determined, true is returned **/
  inline bool sharing_supported(const int i)
    { return (_properties[i].computeMode == CU_COMPUTEMODE_DEFAULT); }

  /// True if splitting device into equal subdevices supported
  inline bool fission_equal()
    { return fission_equal(_device); }
  /// True if splitting device into equal subdevices supported
  inline bool fission_equal(const int i)
    { return false; }
  /// True if splitting device into subdevices by specified counts supported
  inline bool fission_by_counts()
    { return fission_by_counts(_device); }
  /// True if splitting device into subdevices by specified counts supported
  inline bool fission_by_counts(const int i)
    { return false; }
  /// True if splitting device into subdevices by affinity domains supported
  inline bool fission_by_affinity()
    { return fission_by_affinity(_device); }
  /// True if splitting device into subdevices by affinity domains supported
  inline bool fission_by_affinity(const int i)
    { return false; }

  /// Maximum number of subdevices allowed from device fission
  inline int max_sub_devices()
    { return max_sub_devices(_device); }
  /// Maximum number of subdevices allowed from device fission
  inline int max_sub_devices(const int i)
    { return 0; }

  /// True if the device supports shuffle intrinsics
  inline bool has_shuffle_support()
    { return has_shuffle_support(_device); }
  /// True if the device supports shuffle intrinsics
  inline bool has_shuffle_support(const int i)
    { return arch(i)>=3.0; }

  /// List all devices along with all properties
  inline void print_all(std::ostream &out);

  /// For compatability with OCL API
  inline int auto_set_platform(const enum UCL_DEVICE_TYPE type=UCL_GPU,
			       const std::string vendor="",
			       const int ndevices=-1,
			       const int first_device=-1)
    { return set_platform(0); }

 private:
  int _device, _num_devices;
  std::vector<NVDProperties> _properties;
  std::vector<CUstream> _cq;
  CUdevice _cu_device;
  CUcontext _context;
};

// Grabs the properties for all devices
UCL_Device::UCL_Device() {
  CU_SAFE_CALL_NS(cuInit(0));
  CU_SAFE_CALL_NS(cuDeviceGetCount(&_num_devices));
  for (int i=0; i<_num_devices; ++i) {
    CUdevice dev;
    CU_SAFE_CALL_NS(cuDeviceGet(&dev,i));
    int major, minor;
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, dev));
    if (major==9999)
      continue;

    NVDProperties prop;
    prop.device_id = i;
    prop.major=major;
    prop.minor=minor;

    char namecstr[1024];
    CU_SAFE_CALL_NS(cuDeviceGetName(namecstr,1024,dev));
    prop.name=namecstr;

    CU_SAFE_CALL_NS(cuDeviceTotalMem(&prop.totalGlobalMem,dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.multiProcessorCount, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, dev));

    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxThreadsPerBlock, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxThreadsDim[0], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxThreadsDim[1], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxThreadsDim[2], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxGridSize[0], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxGridSize[1], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.maxGridSize[2], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.sharedMemPerBlock, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.totalConstantMemory, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.SIMDWidth, CU_DEVICE_ATTRIBUTE_WARP_SIZE, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.memPitch, CU_DEVICE_ATTRIBUTE_MAX_PITCH, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.regsPerBlock, CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.clockRate, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.textureAlign, CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT, dev));

    #if CUDA_VERSION >= 2020
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.kernelExecTimeoutEnabled, CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT,dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.integrated, CU_DEVICE_ATTRIBUTE_INTEGRATED, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.canMapHostMemory, CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.computeMode, CU_DEVICE_ATTRIBUTE_COMPUTE_MODE,dev));
    #endif
    #if CUDA_VERSION >= 3010
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.concurrentKernels, CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS, dev));
    CU_SAFE_CALL_NS(cuDeviceGetAttribute(&prop.ECCEnabled, CU_DEVICE_ATTRIBUTE_ECC_ENABLED, dev));
    #endif

    _properties.push_back(prop);
  }
  _device=-1;
  _cq.push_back(CUstream());
  _cq.back()=0;
}

UCL_Device::~UCL_Device() {
  clear();
}

int UCL_Device::set_platform(const int pid) {
  clear();
  #ifdef UCL_DEBUG
  assert(pid<num_platforms());
  #endif
  return UCL_SUCCESS;
}

// Set the CUDA device to the specified device number
int UCL_Device::set(int num) {
  clear();
  _device=_properties[num].device_id;
  CU_SAFE_CALL_NS(cuDeviceGet(&_cu_device,_device));
  CUresult err=cuCtxCreate(&_context,0,_cu_device);
  if (err!=CUDA_SUCCESS) {
    #ifndef UCL_NO_EXIT
    std::cerr << "UCL Error: Could not access accelerator number " << num
              << " for use.\n";
    UCL_GERYON_EXIT;
    #endif
    return UCL_ERROR;
  }
  return UCL_SUCCESS;
}

void UCL_Device::clear() {
  if (_device>-1) {
    for (int i=1; i<num_queues(); i++) pop_command_queue();
    cuCtxDestroy(_context);
  }
  _device=-1;
}

// List all devices along with all properties
void UCL_Device::print_all(std::ostream &out) {
  #if CUDA_VERSION >= 2020
  int driver_version;
  cuDriverGetVersion(&driver_version);
  out << "CUDA Driver Version:                           "
      << driver_version/1000 << "." << driver_version%100
                  << std::endl;
  #endif

  if (num_devices() == 0)
    out << "There is no device supporting CUDA\n";
  for (int i=0; i<num_devices(); ++i) {
    out << "\nDevice " << i << ": \"" << name(i) << "\"\n";
    out << "  Type of device:                                "
        << device_type_name(i).c_str() << std::endl;
    out << "  Compute capability:                            "
        << arch(i) << std::endl;
    out << "  Double precision support:                      ";
    if (double_precision(i))
      out << "Yes\n";
    else
      out << "No\n";
    out << "  Total amount of global memory:                 "
        << gigabytes(i) << " GB\n";
    #if CUDA_VERSION >= 2000
    out << "  Number of compute units/multiprocessors:       "
        << _properties[i].multiProcessorCount << std::endl;
    out << "  Number of cores:                               "
        << cores(i) << std::endl;
    #endif
    out << "  Total amount of constant memory:               "
        << _properties[i].totalConstantMemory << " bytes\n";
    out << "  Total amount of local/shared memory per block: "
        << _properties[i].sharedMemPerBlock << " bytes\n";
    out << "  Total number of registers available per block: "
        << _properties[i].regsPerBlock << std::endl;
    out << "  Warp size:                                     "
        << _properties[i].SIMDWidth << std::endl;
    out << "  Maximum number of threads per block:           "
        << _properties[i].maxThreadsPerBlock << std::endl;
    out << "  Maximum group size (# of threads per block)    "
        << _properties[i].maxThreadsDim[0] << " x "
        << _properties[i].maxThreadsDim[1] << " x "
        << _properties[i].maxThreadsDim[2] << std::endl;
    out << "  Maximum item sizes (# threads for each dim)    "
        << _properties[i].maxGridSize[0] << " x "
        << _properties[i].maxGridSize[1] << " x "
        << _properties[i].maxGridSize[2] << std::endl;
    out << "  Maximum memory pitch:                          "
        << max_pitch(i) << " bytes\n";
    out << "  Texture alignment:                             "
        << _properties[i].textureAlign << " bytes\n";
    out << "  Clock rate:                                    "
        << clock_rate(i) << " GHz\n";
    #if CUDA_VERSION >= 2020
    out << "  Run time limit on kernels:                     ";
    if (_properties[i].kernelExecTimeoutEnabled)
      out << "Yes\n";
    else
      out << "No\n";
    out << "  Integrated:                                    ";
    if (_properties[i].integrated)
      out << "Yes\n";
    else
      out << "No\n";
    out << "  Support host page-locked memory mapping:       ";
    if (_properties[i].canMapHostMemory)
      out << "Yes\n";
    else
      out << "No\n";
    out << "  Compute mode:                                  ";
    if (_properties[i].computeMode == CU_COMPUTEMODE_DEFAULT)
      out << "Default\n"; // multiple threads can use device
#if CUDA_VERSION >= 8000
    else if (_properties[i].computeMode == CU_COMPUTEMODE_EXCLUSIVE_PROCESS)
#else
    else if (_properties[i].computeMode == CU_COMPUTEMODE_EXCLUSIVE)
#endif
      out << "Exclusive\n"; // only thread can use device
    else if (_properties[i].computeMode == CU_COMPUTEMODE_PROHIBITED)
      out << "Prohibited\n"; // no thread can use device
    #if CUDART_VERSION >= 4000
    else if (_properties[i].computeMode == CU_COMPUTEMODE_EXCLUSIVE_PROCESS)
      out << "Exclusive Process\n"; // multiple threads 1 process
    #endif
    else
      out << "Unknown\n";
    #endif
    #if CUDA_VERSION >= 3010
    out << "  Concurrent kernel execution:                   ";
    if (_properties[i].concurrentKernels)
      out << "Yes\n";
    else
      out << "No\n";
    out << "  Device has ECC support enabled:                ";
    if (_properties[i].ECCEnabled)
      out << "Yes\n";
    else
      out << "No\n";
    #endif
  }
}

}

#endif
