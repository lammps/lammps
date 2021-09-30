/***************************************************************************
                                ocl_device.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with OpenCL devices

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon Dec 23 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef OCL_DEVICE
#define OCL_DEVICE

#include <string>
#include <vector>
#include <iostream>

#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 300
#endif

#ifdef __APPLE__
#include <OpenCL/cl.h>
#include <OpenCL/cl_platform.h>
#else
#include <CL/cl.h>
#include <CL/cl_platform.h>
#endif

#include "ocl_macros.h"
#include "ucl_types.h"

namespace ucl_opencl {

// --------------------------------------------------------------------------
// - COMMAND QUEUE STUFF
// --------------------------------------------------------------------------
typedef cl_command_queue command_queue;
typedef cl_context context_type;

inline void ucl_flush(command_queue &cq) { CL_SAFE_CALL(clFlush(cq)); }

inline void ucl_sync(cl_command_queue &cq) {
  CL_SAFE_CALL(clFinish(cq));
}

#if defined(GERYON_FORCE_SHARED_MAIN_MEM_ON)
inline bool _shared_mem_device(cl_device_id &device) { return true; }
#elif defined(GERYON_FORCE_SHARED_MAIN_MEM_OFF)
inline bool _shared_mem_device(cl_device_id &device) { return false; }
#else
inline bool _shared_mem_device(cl_device_id &device) {
  #ifdef CL_VERSION_1_2
  cl_bool br;
  CL_SAFE_CALL(clGetDeviceInfo(device, CL_DEVICE_HOST_UNIFIED_MEMORY,
                               sizeof(cl_bool), &br,NULL));
  return (br == CL_TRUE);
  #else
  cl_device_type device_type;
  CL_SAFE_CALL(clGetDeviceInfo(device,CL_DEVICE_TYPE,
                               sizeof(device_type),&device_type,NULL));
  return (device_type==CL_DEVICE_TYPE_CPU);
  #endif
}
#endif

struct OCLProperties {
  std::string name;
  cl_device_type device_type;
  bool is_subdevice;
  cl_ulong global_mem;
  cl_ulong shared_mem;
  cl_ulong const_mem;
  cl_uint compute_units;
  cl_uint clock;
  size_t work_group_size;
  size_t work_item_size[3];
  bool double_precision;
  int preferred_vector_width32, preferred_vector_width64;
  int alignment;
  size_t timer_resolution;
  bool ecc_support;
  std::string c_version;
  bool partition_equal, partition_counts, partition_affinity;
  cl_uint max_sub_devices;
  int cl_device_version;
  bool has_subgroup_support;
  bool has_shuffle_support;
};

/// Class for looking at data parallel device properties
/** \note Calls to change the device outside of the class results in incorrect
  *       behavior
  * \note There is no error checking for indexing past the number of devices **/
class UCL_Device {
 public:
  /// Collect properties for every device on the node
   /** \note You must set the active GPU with set() before using the device **/
  inline UCL_Device();

  inline ~UCL_Device();

  /// Return the number of platforms (0 if error or no platforms)
  inline int num_platforms() { return _num_platforms; }

  /// Return a string with name and info of the current platform
  inline std::string platform_name();

  /// Delete any contexts/data and set the platform number to be used
  inline int set_platform(const int pid);

  /// Return the number of devices that support OpenCL
  inline int num_devices() { return _num_devices; }

  /// Set the OpenCL device to the specified device number
  /** A context and default command queue will be created for the device *
    * Returns UCL_SUCCESS if successful or UCL_ERROR if the device could not
    * be allocated for use. clear() is called to delete any contexts and
    * associated data from previous calls to set(). **/
  inline int set(int num);

  /// Delete any context and associated data stored from a call to set()
  inline void clear();

  /// Get the current device number
  inline int device_num() { return _device; }

  /// Returns the context for the current device
  inline cl_context & context() { return _context; }

  /// Returns the default stream for the current device
  inline command_queue & cq() { return cq(_default_cq); }

  /// Returns the stream indexed by i
  inline command_queue & cq(const int i) { return _cq[i]; }

  /// Set the default command queue
  /** \param i index of the command queue (as added by push_command_queue())
      If i is 0, the command queue created with device initialization is
      used **/
  inline void set_command_queue(const int i) { _default_cq=i; }

  /// Block until all commands in the default stream have completed
  inline void sync() { sync(_default_cq); }

  /// Block until all commands in the specified stream have completed
  inline void sync(const int i) { ucl_sync(cq(i)); }

  /// Get the number of command queues currently available on device
  inline int num_queues()
    { return _cq.size(); }

  /// Add a command queue for device computations (with profiling enabled)
  inline void push_command_queue() {
    cl_int errorv;
    _cq.push_back(cl_command_queue());

#ifdef CL_VERSION_2_0
    cl_queue_properties props[] = {CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0};
    _cq.back()=clCreateCommandQueueWithProperties(_context, _cl_device, props, &errorv);
#else
    _cq.back()=clCreateCommandQueue(_context, _cl_device, CL_QUEUE_PROFILING_ENABLE, &errorv);
#endif
    if (errorv!=CL_SUCCESS) {
      std::cerr << "Could not create command queue on device: " << name()
                << std::endl;
      UCL_GERYON_EXIT;
    }
  }

  /// Remove a stream for device computations
  /** \note You cannot delete the default stream **/
  inline void pop_command_queue() {
    if (_cq.size()<2) return;
    CL_SAFE_CALL(clReleaseCommandQueue(_cq.back()));
    _cq.pop_back();
  }

  /// Get the current OpenCL device name
  inline std::string name() { return name(_device); }
  /// Get the OpenCL device name
  inline std::string name(const int i) {
    return std::string(_properties[i].name); }

  /// Get a string telling the type of the current device
  inline std::string device_type_name() { return device_type_name(_device); }
  /// Get a string telling the type of the device
  inline std::string device_type_name(const int i);

  /// Get current device type (UCL_CPU, UCL_GPU, UCL_ACCELERATOR, UCL_DEFAULT)
  inline enum UCL_DEVICE_TYPE device_type() { return device_type(_device); }
  /// Get device type (UCL_CPU, UCL_GPU, UCL_ACCELERATOR, UCL_DEFAULT)
  inline enum UCL_DEVICE_TYPE device_type(const int i);

  /// Returns true if host memory is efficiently addressable from device
  inline bool shared_memory() { return shared_memory(_device); }
  /// Returns true if host memory is efficiently addressable from device
  inline bool shared_memory(const int i)
    { return _shared_mem_device(_cl_devices[i]); }

  /// Returns preferred vector width
  inline int preferred_fp32_width() { return preferred_fp32_width(_device); }
  /// Returns preferred vector width
  inline int preferred_fp32_width(const int i)
    {return _properties[i].preferred_vector_width32;}
  /// Returns preferred vector width
  inline int preferred_fp64_width() { return preferred_fp64_width(_device); }
  /// Returns preferred vector width
  inline int preferred_fp64_width(const int i)
    {return _properties[i].preferred_vector_width64;}

  /// Returns true if double precision is support for the current device
  inline bool double_precision() { return double_precision(_device); }
  /// Returns true if double precision is support for the device
  inline bool double_precision(const int i)
    {return _properties[i].double_precision;}

  /// Get the number of compute units on the current device
  inline unsigned cus() { return cus(_device); }
  /// Get the number of compute units
  inline unsigned cus(const int i)
    { return _properties[i].compute_units; }

  /// Get the gigabytes of global memory in the current device
  inline double gigabytes() { return gigabytes(_device); }
  /// Get the gigabytes of global memory
  inline double gigabytes(const int i)
    { return static_cast<double>(_properties[i].global_mem)/1073741824; }

  /// Get the bytes of global memory in the current device
  inline size_t bytes() { return bytes(_device); }
  /// Get the bytes of global memory
  inline size_t bytes(const int i) { return _properties[i].global_mem; }

  /// Return the GPGPU revision number for current device
  //inline double revision() { return revision(_device); }
  /// Return the GPGPU revision number
  //inline double revision(const int i)
  //  { return //static_cast<double>(_properties[i].minor)/10+_properties[i].major;}

  /// Clock rate in GHz for current device
  inline double clock_rate() { return clock_rate(_device); }
  /// Clock rate in GHz
  inline double clock_rate(const int i) { return _properties[i].clock*1e-3;}

  /// Return the address alignment in bytes
  inline int alignment() { return alignment(_device); }
  /// Return the address alignment in bytes
  inline int alignment(const int i) { return _properties[i].alignment; }

  /// Return the timer resolution
  inline size_t timer_resolution() { return timer_resolution(_device); }
  /// Return the timer resolution
  inline size_t timer_resolution(const int i)
    { return _properties[i].timer_resolution; }

  /// Get the maximum number of threads per block
  inline size_t group_size() { return group_size(_device); }
  /// Get the maximum number of threads per block
  inline size_t group_size(const int i)
    { return _properties[i].work_group_size; }
  /// Get the maximum number of threads per block in dimension 'dim'
  inline size_t group_size_dim(const int dim)
    { return group_size_dim(_device, dim); }
  /// Get the maximum number of threads per block in dimension 'dim'
  inline size_t group_size_dim(const int i, const int dim)
    { return _properties[i].work_item_size[dim]; }

  /// Get the shared local memory size in bytes
  inline size_t slm_size() { return slm_size(_device); }
  /// Get the shared local memory size in bytes
  inline size_t slm_size(const int i)
    { return _properties[i].shared_mem; }

  /// Return the maximum memory pitch in bytes for current device
  inline size_t max_pitch() { return max_pitch(_device); }
  /// Return the maximum memory pitch in bytes
  inline size_t max_pitch(const int i) { return 0; }

  /// Returns false if accelerator cannot be shared by multiple processes
  /** If it cannot be determined, true is returned **/
  inline bool sharing_supported() { return sharing_supported(_device); }
  /// Returns false if accelerator cannot be shared by multiple processes
  /** If it cannot be determined, true is returned **/
  inline bool sharing_supported(const int i)
    { return true; }

  /// True if the device is a sub-device
  inline bool is_subdevice()
    { return is_subdevice(_device); }
  /// True if the device is a sub-device
  inline bool is_subdevice(const int i)
    { return _properties[i].is_subdevice; }
  /// True if splitting device into equal subdevices supported
  inline bool fission_equal()
    { return fission_equal(_device); }
  /// True if splitting device into equal subdevices supported
  inline bool fission_equal(const int i)
    { return _properties[i].partition_equal; }
  /// True if splitting device into subdevices by specified counts supported
  inline bool fission_by_counts()
    { return fission_by_counts(_device); }
  /// True if splitting device into subdevices by specified counts supported
  inline bool fission_by_counts(const int i)
    { return _properties[i].partition_counts; }
  /// True if splitting device into subdevices by affinity domains supported
  inline bool fission_by_affinity()
    { return fission_by_affinity(_device); }
  /// True if splitting device into subdevices by affinity domains supported
  inline bool fission_by_affinity(const int i)
    { return _properties[i].partition_affinity; }
  /// True if the device has subgroup support
  inline bool has_subgroup_support()
    { return has_subgroup_support(_device); }
  /// True if the device has subgroup support
  inline bool has_subgroup_support(const int i)
    { return _properties[i].has_subgroup_support; }
  /// True if the device supports shuffle intrinsics
  inline bool has_shuffle_support()
    { return has_shuffle_support(_device); }
  /// True if the device supports shuffle intrinsics
  inline bool has_shuffle_support(const int i)
    { return _properties[i].has_shuffle_support; }

  /// Maximum number of subdevices allowed from device fission
  inline int max_sub_devices()
    { return max_sub_devices(_device); }
  /// Maximum number of subdevices allowed from device fission
  inline int max_sub_devices(const int i)
    { return _properties[i].max_sub_devices; }
  /// OpenCL version supported by the device
  inline int cl_device_version()
    { return cl_device_version(_device); }
  /// OpenCL version supported by the device
  inline int cl_device_version(const int i)
    { return _properties[i].cl_device_version; }

  /// List all devices along with all properties
  inline void print_all(std::ostream &out);

  /// Return the OpenCL type for the device
  inline cl_device_id & cl_device() { return _cl_device; }

  /// Automatically set the platform by type, vendor, and/or CU count
  /** If first_device is positive, search restricted to platforms containing
    * this device IDs. If ndevices is positive, search is restricted
    * to platforms with at least that many devices  **/
  inline int auto_set_platform(const enum UCL_DEVICE_TYPE type=UCL_GPU,
                               const std::string vendor="",
                               const int ndevices=-1,
                               const int first_device=-1);

 private:
  int _num_platforms;          // Number of platforms
  int _platform;               // UCL_Device ID for current platform
  cl_platform_id _cl_platform; // OpenCL ID for current platform
  cl_platform_id _cl_platforms[20]; // OpenCL IDs for all platforms
  cl_context _context;              // Context used for accessing the device
  std::vector<cl_command_queue> _cq;// The default command queue for this device
  int _device;                            // UCL_Device ID for current device
  cl_device_id _cl_device;                // OpenCL ID for current device
  std::vector<cl_device_id> _cl_devices;  // OpenCL IDs for all devices
  int _num_devices;                       // Number of devices
  std::vector<OCLProperties> _properties; // Properties for each device

  inline void add_properties(cl_device_id);
  inline int create_context();
  int _default_cq;
};

// Grabs the properties for all devices
UCL_Device::UCL_Device() {
  _device=-1;

  // --- Get Number of Platforms
  cl_uint nplatforms;
  cl_int errorv=clGetPlatformIDs(20,_cl_platforms,&nplatforms);

  if (errorv!=CL_SUCCESS) {
    _num_platforms=0;
    return;
  } else
    _num_platforms=static_cast<int>(nplatforms);
  set_platform(0);
}

UCL_Device::~UCL_Device() {
  clear();
}

void UCL_Device::clear() {
  _properties.clear();

  #ifdef GERYON_NUMA_FISSION
  #ifdef CL_VERSION_1_2
  for (int i=0; i<_cl_devices.size(); i++)
    CL_DESTRUCT_CALL(clReleaseDevice(_cl_devices[i]));
  #endif
  #endif

  _cl_devices.clear();
  if (_device>-1) {
    for (size_t i=0; i<_cq.size(); i++) {
      CL_DESTRUCT_CALL(clReleaseCommandQueue(_cq.back()));
      _cq.pop_back();
    }
    CL_DESTRUCT_CALL(clReleaseContext(_context));
  }
  _device=-1;
  _num_devices=0;
}

int UCL_Device::set_platform(int pid) {
  clear();
  cl_int errorv;

  _cl_device=0;
  _device=-1;
  _num_devices=0;
  _default_cq=0;

  #ifdef UCL_DEBUG
  assert(pid<num_platforms());
  #endif
  _platform=pid;
  _cl_platform=_cl_platforms[_platform];

  // --- Get Number of Devices
  cl_uint n;
  errorv=clGetDeviceIDs(_cl_platform,CL_DEVICE_TYPE_ALL,0,nullptr,&n);
  _num_devices=n;
  if (errorv!=CL_SUCCESS || _num_devices==0) {
    _num_devices=0;
    return UCL_ERROR;
  }
  cl_device_id *device_list = new cl_device_id[_num_devices];
  CL_SAFE_CALL(clGetDeviceIDs(_cl_platform,CL_DEVICE_TYPE_ALL,n,device_list,
                              &n));

  #ifndef GERYON_NUMA_FISSION
  // --- Store properties for each device
  for (int i=0; i<_num_devices; i++) {
    _cl_devices.push_back(device_list[i]);
    add_properties(device_list[i]);
  }
  #else
  // --- Create sub-devices for anything partitionable by NUMA and store props
  int num_unpart = _num_devices;
  _num_devices = 0;
  for (int i=0; i<num_unpart; i++) {
    cl_uint num_subdevices = 1;

    #ifdef CL_VERSION_1_2
    cl_device_affinity_domain adomain;
    CL_SAFE_CALL(clGetDeviceInfo(device_list[i],
                                 CL_DEVICE_PARTITION_AFFINITY_DOMAIN,
                                 sizeof(cl_device_affinity_domain),
                                 &adomain,NULL));

    cl_device_partition_property props[3];
    props[0]=CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN;
    props[1]=CL_DEVICE_AFFINITY_DOMAIN_NUMA;
    props[2]=0;
    if (adomain & CL_DEVICE_AFFINITY_DOMAIN_NUMA)
      CL_SAFE_CALL(clCreateSubDevices(device_list[i], props, 0, NULL,
                                      &num_subdevices));
    if (num_subdevices > 1) {
      cl_device_id *subdevice_list = new cl_device_id[num_subdevices];
      CL_SAFE_CALL(clCreateSubDevices(device_list[i], props, num_subdevices,
                                      subdevice_list, &num_subdevices));
      for (int j=0; j<num_subdevices; j++) {
        _cl_devices.push_back(device_list[i]);
        add_properties(device_list[i]);
        _num_devices++;
      }
      delete[] subdevice_list;
    } else {
      _cl_devices.push_back(device_list[i]);
      add_properties(device_list[i]);
      _num_devices++;
    }
    #endif
  } // for i
  #endif

  delete[] device_list;
  return UCL_SUCCESS;
}

int UCL_Device::create_context() {
  cl_int errorv;
  cl_context_properties props[3];
  props[0]=CL_CONTEXT_PLATFORM;
  props[1]=_platform;
  props[2]=0;
  _context=clCreateContext(0,1,&_cl_device,nullptr,nullptr,&errorv);
  if (errorv!=CL_SUCCESS) {
    #ifndef UCL_NO_EXIT
    std::cerr << "UCL Error: Could not access accelerator number " << _device
              << " for use.\n";
    UCL_GERYON_EXIT;
    #endif
    return UCL_ERROR;
  }
  push_command_queue();
  _default_cq=0;
  return UCL_SUCCESS;
}

void UCL_Device::add_properties(cl_device_id device_list) {
  OCLProperties op;
  char buffer[1024];
  cl_bool ans_bool;

  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_NAME,1024,buffer,nullptr));
  op.name=buffer;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_GLOBAL_MEM_SIZE,
                               sizeof(op.global_mem),&op.global_mem,nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_LOCAL_MEM_SIZE,
                               sizeof(op.shared_mem),&op.shared_mem,nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
                               sizeof(op.const_mem),&op.const_mem,nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_TYPE,
                               sizeof(op.device_type),&op.device_type,nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MAX_COMPUTE_UNITS,
                               sizeof(op.compute_units),&op.compute_units,
                               nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MAX_CLOCK_FREQUENCY,
                               sizeof(op.clock),&op.clock,nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MAX_WORK_GROUP_SIZE,
                               sizeof(op.work_group_size),&op.work_group_size,
                               nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MAX_WORK_ITEM_SIZES,
                               3*sizeof(op.work_item_size[0]),op.work_item_size,
                               nullptr));
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_MEM_BASE_ADDR_ALIGN,
                               sizeof(cl_uint),&op.alignment,nullptr));
  op.alignment/=8;

  cl_uint float_width;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,
                               sizeof(float_width),&float_width,nullptr));
  op.preferred_vector_width32=float_width;

  cl_uint double_width;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,
                               sizeof(double_width),&double_width,nullptr));
  op.preferred_vector_width64=double_width;

  // Determine if double precision is supported: All bits in the mask must be set.
  cl_device_fp_config double_mask = (CL_FP_FMA|CL_FP_ROUND_TO_NEAREST|CL_FP_ROUND_TO_ZERO|
                                     CL_FP_ROUND_TO_INF|CL_FP_INF_NAN|CL_FP_DENORM);
  cl_device_fp_config double_avail;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_DOUBLE_FP_CONFIG,
                               sizeof(double_avail),&double_avail,nullptr));
  if ((double_avail & double_mask) == double_mask)
    op.double_precision=true;
  else
    op.double_precision=false;

  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_PROFILING_TIMER_RESOLUTION,
                               sizeof(size_t),&op.timer_resolution,nullptr));


  op.ecc_support=false;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_ERROR_CORRECTION_SUPPORT,
                               sizeof(ans_bool),&ans_bool,nullptr));
  if (ans_bool==CL_TRUE)
    op.ecc_support=true;

  op.c_version="";
  op.is_subdevice=false;
  op.partition_equal=false;
  op.partition_counts=false;
  op.partition_affinity=false;
  op.max_sub_devices=1;
  op.cl_device_version=0;
  op.has_subgroup_support=false;
  op.has_shuffle_support=false;

  #ifdef CL_VERSION_1_2
  size_t return_bytes;
  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_OPENCL_C_VERSION,1024,
                               buffer,nullptr));
  op.c_version=buffer;

  cl_device_partition_property pinfo[4];
  CL_SAFE_CALL(clGetDeviceInfo(device_list, CL_DEVICE_PARTITION_TYPE,
                               4*sizeof(cl_device_partition_property),
                               &pinfo, &return_bytes));
  if (return_bytes == 0) op.is_subdevice=false;
  else if (pinfo[0]) op.is_subdevice=true;
  else op.is_subdevice=false;

  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_PARTITION_PROPERTIES,
                               4*sizeof(cl_device_partition_property),
                               pinfo,&return_bytes));
  int nprops=return_bytes/sizeof(cl_device_partition_property);
  for (int i=0; i<nprops; i++) {
    if (pinfo[i]==CL_DEVICE_PARTITION_EQUALLY)
      op.partition_equal=true;
    else if (pinfo[i]==CL_DEVICE_PARTITION_BY_COUNTS)
      op.partition_counts=true;
    else if (pinfo[i]==CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN)
      op.partition_affinity=true;
  }

  CL_SAFE_CALL(clGetDeviceInfo(device_list,
                               CL_DEVICE_PARTITION_MAX_SUB_DEVICES,
                               sizeof(cl_uint),&op.max_sub_devices,nullptr));

  CL_SAFE_CALL(clGetDeviceInfo(device_list,CL_DEVICE_VERSION,1024,buffer,nullptr));
  int cl_version_maj = buffer[7] - '0';
  int cl_version_min = buffer[9] - '0';
  op.cl_device_version = cl_version_maj * 100 + cl_version_min * 10;

  size_t ext_str_size_ret;
  CL_SAFE_CALL(clGetDeviceInfo(device_list, CL_DEVICE_EXTENSIONS, 0, nullptr,
                               &ext_str_size_ret));
  char buffer2[ext_str_size_ret];
  CL_SAFE_CALL(clGetDeviceInfo(device_list, CL_DEVICE_EXTENSIONS,
                               ext_str_size_ret, buffer2, nullptr));
  #if defined(CL_VERSION_2_1) || defined(CL_VERSION_3_0)
  if (op.cl_device_version >= 210) {
    if ((std::string(buffer2).find("cl_khr_subgroups") != std::string::npos) ||
        (std::string(buffer2).find("cl_intel_subgroups") != std::string::npos))
      op.has_subgroup_support=true;
    if (std::string(buffer2).find("cl_intel_subgroups") != std::string::npos)
      op.has_shuffle_support=true;
  }
  #endif
  if (std::string(buffer2).find("cl_nv_device_attribute_query") !=
      std::string::npos) {
    #ifndef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
    #define CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV 0x4000
    #endif
    #ifndef CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
    #define CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV 0x4001
    #endif
    cl_uint major, minor;
    CL_SAFE_CALL(clGetDeviceInfo(device_list,
                                 CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,
                                 sizeof(cl_uint), &major, nullptr));
    CL_SAFE_CALL(clGetDeviceInfo(device_list,
                                 CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,
                                 sizeof(cl_uint), &minor, nullptr));
    double arch = static_cast<double>(minor)/10+major;
    if (arch >= 3.0)
      op.has_shuffle_support=true;
  }
  #endif

  _properties.push_back(op);
}

std::string UCL_Device::platform_name() {
  char info[1024];

  CL_SAFE_CALL(clGetPlatformInfo(_cl_platform,CL_PLATFORM_VENDOR,1024,info,
                                 nullptr));
  std::string ans=std::string(info)+' ';

  CL_SAFE_CALL(clGetPlatformInfo(_cl_platform,CL_PLATFORM_NAME,1024,info,
                                 nullptr));
  ans+=std::string(info)+' ';

  CL_SAFE_CALL(clGetPlatformInfo(_cl_platform,CL_PLATFORM_VERSION,1024,info,
               nullptr));
  ans+=std::string(info);

  return ans;
}

// Get a string telling the type of the device
std::string UCL_Device::device_type_name(const int i) {
  if (_properties[i].device_type==CL_DEVICE_TYPE_CPU)
    return "CPU";
  else if (_properties[i].device_type==CL_DEVICE_TYPE_GPU)
    return "GPU";
  else if (_properties[i].device_type==CL_DEVICE_TYPE_ACCELERATOR)
    return "ACCELERATOR";
  else
    return "DEFAULT";
}

// Get a string telling the type of the device
enum UCL_DEVICE_TYPE UCL_Device::device_type(const int i) {
  if (_properties[i].device_type==CL_DEVICE_TYPE_CPU)
    return UCL_CPU;
  else if (_properties[i].device_type==CL_DEVICE_TYPE_GPU)
    return UCL_GPU;
  else if (_properties[i].device_type==CL_DEVICE_TYPE_ACCELERATOR)
    return UCL_ACCELERATOR;
  else
    return UCL_DEFAULT;
}

// Set the CUDA device to the specified device number
int UCL_Device::set(int num) {
  _device=num;
  _cl_device=_cl_devices[_device];
  return create_context();
}

// List all devices from all platforms along with all properties
void UCL_Device::print_all(std::ostream &out) {
  // --- loop through the platforms
  for (int n=0; n<_num_platforms; n++) {

    set_platform(n);

    out << "\nPlatform " << n << ":\n";

    if (num_devices() == 0)
      out << "There is no device supporting OpenCL\n";
    for (int i=0; i<num_devices(); ++i) {
      out << "\nDevice " << i << ": \"" << name(i).c_str() << "\"\n";
      out << "  Type of device:                                "
          << device_type_name(i).c_str() << std::endl;
      out << "  Supported OpenCL Version:                      "
          << _properties[i].cl_device_version / 100 << "."
          << _properties[i].cl_device_version % 100 << std::endl;
      out << "  Is a subdevice:                                ";
      if (is_subdevice(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Double precision support:                      ";
      if (double_precision(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Total amount of global memory:                 "
          << gigabytes(i) << " GB\n";
      out << "  Number of compute units/multiprocessors:       "
          << _properties[i].compute_units << std::endl;
      //out << "  Number of cores:                               "
      //    << cores(i) << std::endl;
      out << "  Total amount of constant memory:               "
          << _properties[i].const_mem << " bytes\n";
      out << "  Total amount of local/shared memory per block: "
          << _properties[i].shared_mem << " bytes\n";
      //out << "  Total number of registers available per block: "
      //    << _properties[i].regsPerBlock << std::endl;
      //out << "  Warp size:                                     "
      //    << _properties[i].warpSize << std::endl;
      out << "  Maximum group size (# of threads per block)    "
          << _properties[i].work_group_size << std::endl;
      out << "  Maximum item sizes (# threads for each dim)    "
          << _properties[i].work_item_size[0] << " x "
          << _properties[i].work_item_size[1] << " x "
          << _properties[i].work_item_size[2] << std::endl;
      //out << "  Maximum sizes of each dimension of a grid:     "
      //    << _properties[i].maxGridSize[0] << " x "
      //    << _properties[i].maxGridSize[1] << " x "
      //    << _properties[i].maxGridSize[2] << std::endl;
      //out << "  Maximum memory pitch:                          "
      //    << _properties[i].memPitch) << " bytes\n";
      //out << "  Texture alignment:                             "
      //    << _properties[i].textureAlignment << " bytes\n";
      out << "  Clock rate:                                    "
          << clock_rate(i) << " GHz\n";
      //out << "  Concurrent copy and execution:                 ";
      out << "  ECC support:                                   ";
      if (_properties[i].ecc_support)
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Device fission into equal partitions:          ";
      if (fission_equal(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Device fission by counts:                      ";
      if (fission_by_counts(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Device fission by affinity:                    ";
      if (fission_by_affinity(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Maximum subdevices from fission:               "
          << max_sub_devices(i) << std::endl;
      out << "  Shared memory system:                          ";
      if (shared_memory(i))
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Subgroup support:                              ";
      if (_properties[i].has_subgroup_support)
        out << "Yes\n";
      else
        out << "No\n";
      out << "  Shuffle support:                               ";
      if (_properties[i].has_shuffle_support)
        out << "Yes\n";
      else
        out << "No\n";
    }
  }
}

int UCL_Device::auto_set_platform(const enum UCL_DEVICE_TYPE type,
                                  const std::string vendor,
                                  const int ndevices,
                                  const int first_device) {
  if (_num_platforms < 2) return set_platform(0);

  int last_device = -1;
  if (first_device > -1) {
    if (ndevices)
      last_device = first_device + ndevices - 1;
    else
      last_device = first_device;
  }

  bool vendor_match=false;
  bool type_match=false;
  int max_cus=0;
  int best_platform=0;

  std::string vendor_upper=vendor;
  for (int i=0; i<vendor.length(); i++)
    if (vendor_upper[i]<='z' && vendor_upper[i]>='a')
      vendor_upper[i]=toupper(vendor_upper[i]);

  for (int n=0; n<_num_platforms; n++) {
    set_platform(n);
    if (last_device > -1 && last_device >= num_devices()) continue;
    if (ndevices > num_devices()) continue;

    int first_id=0;
    int last_id=num_devices()-1;
    if (last_device > -1) {
      first_id=first_device;
      last_id=last_device;
    }

    if (vendor_upper!="") {
      std::string pname = platform_name();
      for (int i=0; i<pname.length(); i++)
        if (pname[i]<='z' && pname[i]>='a')
          pname[i]=toupper(pname[i]);

      if (pname.find(vendor_upper)!=std::string::npos) {
        if (vendor_match == false) {
          best_platform=n;
          max_cus=0;
          vendor_match=true;
        }
      } else if (vendor_match)
        continue;
    }

    if (type != UCL_DEFAULT) {
      bool ptype_matched=false;
      for (int d=first_id; d<=last_id; d++) {
        if (type==device_type(d)) {
          if (type_match == false) {
            best_platform=n;
            max_cus=0;
            type_match=true;
            ptype_matched=true;
          }
        }
      }
      if (type_match==true && ptype_matched==false)
        continue;
    }

    for (int d=first_id; d<=last_id; d++) {
      if (cus(d) > max_cus) {
        best_platform=n;
        max_cus=cus(d);
      }
    }
  }
  return set_platform(best_platform);
}

} // namespace ucl_opencl

#endif
