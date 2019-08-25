/***************************************************************************
                                ocl_kernel.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with OpenCL kernels

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Sun Feb 7 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef OCL_KERNEL
#define OCL_KERNEL

#include "ocl_device.h"
#include <fstream>

namespace ucl_opencl {

class UCL_Texture;
template <class numtyp> class UCL_D_Vec;
template <class numtyp> class UCL_D_Mat;
template <class hosttype, class devtype> class UCL_Vector;
template <class hosttype, class devtype> class UCL_Matrix;
#define UCL_MAX_KERNEL_ARGS 256

/// Class storing 1 or more kernel functions from a single string or file
class UCL_Program {
 public:
  inline UCL_Program() : _init_done(false) {}
  inline UCL_Program(UCL_Device &device) : _init_done(false) { init(device); }
  inline UCL_Program(UCL_Device &device, const void *program,
                     const char *flags="", std::string *log=NULL) :
      _init_done(false) {
    init(device);
    load_string(program,flags,log);
  }

  inline ~UCL_Program() { clear(); }

  /// Initialize the program with a device
  inline void init(UCL_Device &device) {
    clear();
    _device=device.cl_device();
    _context=device.context();
    _cq=device.cq();
    CL_SAFE_CALL(clRetainContext(_context));
    CL_SAFE_CALL(clRetainCommandQueue(_cq));
    _init_done=true;
  }

  /// Clear any data associated with program
  /** \note Must call init() after each clear **/
  inline void clear() {
    if (_init_done) {
      CL_DESTRUCT_CALL(clReleaseProgram(_program));
      CL_DESTRUCT_CALL(clReleaseContext(_context));
      CL_DESTRUCT_CALL(clReleaseCommandQueue(_cq));
      _init_done=false;
    }
  }

  /// Load a program from a file and compile with flags
  inline int load(const char *filename, const char *flags="",
                  std::string *log=NULL) {
    std::ifstream in(filename);
    if (!in || in.is_open()==false) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not open kernel file: "
                << filename << std::endl;
      UCL_GERYON_EXIT;
      #endif
      return UCL_FILE_NOT_FOUND;
    }

    std::string program((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    in.close();
    return load_string(program.c_str(),flags,log);
  }

  /// Load a program from a string and compile with flags
  inline int load_string(const void *program, const char *flags="",
                         std::string *log=NULL) {
    cl_int error_flag;
    const char *prog=(const char *)program;
    _program=clCreateProgramWithSource(_context,1,&prog,NULL,&error_flag);
    CL_CHECK_ERR(error_flag);
    error_flag = clBuildProgram(_program,1,&_device,flags,NULL,NULL);
    if (error_flag!=-11)
      CL_CHECK_ERR(error_flag);
    cl_build_status build_status;
    CL_SAFE_CALL(clGetProgramBuildInfo(_program,_device,
                                       CL_PROGRAM_BUILD_STATUS,
                                       sizeof(cl_build_status),&build_status,
                                       NULL));

    if (build_status != CL_SUCCESS || log!=NULL) {
      size_t ms;
      CL_SAFE_CALL(clGetProgramBuildInfo(_program,_device,CL_PROGRAM_BUILD_LOG,0,
                                         NULL, &ms));
      char *build_log = new char[ms];
      CL_SAFE_CALL(clGetProgramBuildInfo(_program,_device,CL_PROGRAM_BUILD_LOG,ms,
                                         build_log, NULL));

      if (log!=NULL)
        *log=std::string(build_log);

      if (build_status != CL_SUCCESS) {
        #ifndef UCL_NO_EXIT
        std::cerr << std::endl
                  << "----------------------------------------------------------\n"
                  << " UCL Error: Error compiling OpenCL Program ("
                  << build_status << ") ...\n"
                  << "----------------------------------------------------------\n";
        std::cerr << build_log << std::endl;
        #endif
        delete[] build_log;
        return UCL_COMPILE_ERROR;
      } else delete[] build_log;
    }

    return UCL_SUCCESS;
  }

  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return _cq; }
  /// Change the default command queue associated with matrix
  inline void cq(command_queue &cq_in) { _cq=cq_in; }

  friend class UCL_Kernel;
 private:
  bool _init_done;
  cl_program _program;
  cl_device_id _device;
  cl_context _context;
  cl_command_queue _cq;
};

/// Class for dealing with OpenCL kernels
class UCL_Kernel {
 public:
  UCL_Kernel() : _dimensions(1), _function_set(false), _num_args(0)
    {  _block_size[0]=0; _num_blocks[0]=0; }

  inline UCL_Kernel(UCL_Program &program, const char *function) :
    _dimensions(1), _function_set(false), _num_args(0)
    {  _block_size[0]=0; _num_blocks[0]=0; set_function(program,function); }

  inline ~UCL_Kernel() { clear(); }

  /// Clear any function associated with the kernel
  inline void clear() {
    if (_function_set) {
      clReleaseKernel(_kernel);
      clReleaseProgram(_program);
      clReleaseCommandQueue(_cq);
      _function_set=false;
    }
  }

  /// Get the kernel function from a program
  /** \return UCL_ERROR_FLAG (UCL_SUCCESS, UCL_FILE_NOT_FOUND, UCL_ERROR) **/
  inline int set_function(UCL_Program &program, const char *function);

  /// Set the kernel argument.
  /** If not a device pointer, this must be repeated each time the argument
    * changes **/
  template <class dtype>
  inline void set_arg(const cl_uint index, const dtype * const arg) {
    CL_SAFE_CALL(clSetKernelArg(_kernel,index,sizeof(dtype),arg));
    if (index>_num_args) {
      _num_args=index;
      #ifdef UCL_DEBUG
      if (_num_args>_kernel_info_nargs) {
        std::cerr << "TOO MANY ARGUMENTS TO OPENCL FUNCTION: "
                  << _kernel_info_name << std::endl;
        assert(0==1);
      }
      #endif
    }
  }

  /// Set a geryon container as a kernel argument.
  template <class numtyp>
  inline void set_arg(const UCL_D_Vec<numtyp> * const arg)
    { set_arg(&arg->begin()); }

  /// Set a geryon container as a kernel argument.
  template <class numtyp>
  inline void set_arg(const UCL_D_Mat<numtyp> * const arg)
    { set_arg(&arg->begin()); }

  /// Set a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void set_arg(const UCL_Vector<hosttype, devtype> * const arg)
    { set_arg(&arg->device.begin()); }

  /// Set a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void set_arg(const UCL_Matrix<hosttype, devtype> * const arg)
    { set_arg(&arg->device.begin()); }

  /// Add a kernel argument.
  template <class dtype>
  inline void add_arg(const dtype * const arg) {
    CL_SAFE_CALL(clSetKernelArg(_kernel,_num_args,sizeof(dtype),arg));
    _num_args++;
    #ifdef UCL_DEBUG
    if (_num_args>_kernel_info_nargs) {
      std::cerr << "TOO MANY ARGUMENTS TO OPENCL FUNCTION: "
                << _kernel_info_name << std::endl;
      assert(0==1);
    }
    #endif
  }

  /// Add a geryon container as a kernel argument.
  template <class numtyp>
  inline void add_arg(const UCL_D_Vec<numtyp> * const arg)
    { add_arg(&arg->begin()); }

  /// Add a geryon container as a kernel argument.
  template <class numtyp>
  inline void add_arg(const UCL_D_Mat<numtyp> * const arg)
    { add_arg(&arg->begin()); }

  /// Add a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void add_arg(const UCL_Vector<hosttype, devtype> * const arg)
    { add_arg(&arg->device.begin()); }

  /// Add a geryon container as a kernel argument.
  template <class hosttype, class devtype>
  inline void add_arg(const UCL_Matrix<hosttype, devtype> * const arg)
    { add_arg(&arg->device.begin()); }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks, const size_t block_size) {
    _dimensions=1;
    _num_blocks[0]=num_blocks*block_size;
    _block_size[0]=block_size;
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue for the kernel is changed to cq **/
  inline void set_size(const size_t num_blocks, const size_t block_size,
                       command_queue &cq)
    { _cq=cq; set_size(num_blocks,block_size); }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y) {
    _dimensions=2;
    _num_blocks[0]=num_blocks_x*block_size_x;
    _block_size[0]=block_size_x;
    _num_blocks[1]=num_blocks_y*block_size_y;
    _block_size[1]=block_size_y;
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue for the kernel is changed to cq **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y,
                       command_queue &cq)
    {_cq=cq; set_size(num_blocks_x, num_blocks_y, block_size_x, block_size_y);}

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x,
                       const size_t block_size_y, const size_t block_size_z) {
    _dimensions=3;
    const size_t num_blocks_z=1;
    _num_blocks[0]=num_blocks_x*block_size_x;
    _block_size[0]=block_size_x;
    _num_blocks[1]=num_blocks_y*block_size_y;
    _block_size[1]=block_size_y;
    _num_blocks[2]=num_blocks_z*block_size_z;
    _block_size[2]=block_size_z;
  }

  /// Set the number of thread blocks and the number of threads in each block
  /** \note This should be called before any arguments have been added
      \note The default command queue is used for the kernel execution **/
  inline void set_size(const size_t num_blocks_x, const size_t num_blocks_y,
                       const size_t block_size_x, const size_t block_size_y,
                       const size_t block_size_z, command_queue &cq) {
    _cq=cq;
    set_size(num_blocks_x, num_blocks_y, block_size_x, block_size_y,
             block_size_z);
  }

  /// Run the kernel in the default command queue
  inline void run();

  /// Clear any arguments associated with the kernel
  inline void clear_args() { _num_args=0; }

  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return _cq; }
  /// Change the default command queue associated with matrix
  inline void cq(command_queue &cq_in) { _cq=cq_in; }
  #include "ucl_arg_kludge.h"

 private:
  cl_kernel _kernel;
  cl_program _program;
  cl_uint _dimensions;
  size_t _block_size[3];
  size_t _num_blocks[3];
  bool _function_set;

  cl_command_queue _cq;        // The default command queue for this kernel
  unsigned _num_args;

  #ifdef UCL_DEBUG
  std::string _kernel_info_name;
  unsigned _kernel_info_nargs;
  //std::string _kernel_info_args[256];
  #endif
};

inline int UCL_Kernel::set_function(UCL_Program &program, const char *function) {
  clear();
  _function_set=true;
  _cq=program._cq;
  CL_SAFE_CALL(clRetainCommandQueue(_cq));
  _program=program._program;
  CL_SAFE_CALL(clRetainProgram(_program));
  cl_int error_flag;
  _kernel=clCreateKernel(program._program,function,&error_flag);

  if (error_flag!=CL_SUCCESS) {
    #ifndef UCL_NO_EXIT
    std::cerr << "UCL Error: Could not find function: " << function
              << " in program.\n";
    UCL_GERYON_EXIT;
    #endif
    return UCL_FUNCTION_NOT_FOUND;
  }

  #ifdef UCL_DEBUG
  _kernel_info_name=function;
  cl_uint nargs;
  CL_SAFE_CALL(clGetKernelInfo(_kernel,CL_KERNEL_NUM_ARGS,sizeof(cl_uint),
                               &nargs,NULL));
  _kernel_info_nargs=nargs;
  #ifdef NOT_TEST_CL_VERSION_1_2
  char tname[256];
  size_t ret;
  for (cl_uint i=0; i<nargs; i++) {
    CL_SAFE_CALL(clGetKernelArgInfo(_kernel,i,CL_KERNEL_ARG_TYPE_NAME,256,
                                    tname,&ret));
    _kernel_info_args[i]=tname;
  }
  #endif
  #endif

  return UCL_SUCCESS;
}

void UCL_Kernel::run() {
  CL_SAFE_CALL(clEnqueueNDRangeKernel(_cq,_kernel,_dimensions,NULL,
                                      _num_blocks,_block_size,0,NULL,NULL));
}

} // namespace

#endif

